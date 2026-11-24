#!/usr/bin/env python3
"""Exact arithmetic for Note 062's degree-two comparison.

This checks the elementary-series identity and triangular extraction in a
finite range. It does not evaluate higher-genus obstruction integrals or
verify the geometric/transversality proof of the all-genus theorem.
Only the Python standard library is needed.
"""

from fractions import Fraction as Q
from math import factorial


MAX_GENUS = 24
ORDER = 2 * MAX_GENUS - 2


def multiply(left, right):
    result = [Q(0) for _ in range(ORDER + 1)]
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right[: ORDER + 1 - i]):
                if b:
                    result[i + j] += a * b
    return result


def power(series, exponent):
    result = [Q(1)] + [Q(0)] * ORDER
    for _ in range(exponent):
        result = multiply(result, series)
    return result


def sine_series(scale=1):
    """Return 2 sin(scale*t/2), using its Taylor definition."""
    result = [Q(0) for _ in range(ORDER + 1)]
    for j in range((ORDER + 1) // 2):
        n = 2 * j + 1
        result[n] = Q((-1) ** j * scale**n, 2 ** (n - 1) * factorial(n))
    return result


def main():
    s = sine_series()
    s_two = sine_series(2)
    s4, s6 = power(s, 4), power(s, 6)

    # Both partitions of 2 have hooks (1, 2). Subtract the disconnected
    # square to obtain the connected elementary degree-two coefficient.
    two_hooks = multiply(power(s, 2), power(s_two, 2))
    elementary = [2 * x - y / 2 for x, y in zip(two_hooks, s4)]
    assert elementary == [Q(15, 2) * x - 2 * y for x, y in zip(s4, s6)]

    mixed = [x - 8 * y + 2 * z for x, y, z in zip(elementary, s4, s6)]
    assert mixed == [-x / 2 for x in s4]

    # A separate cosine expansion gives -s(t)^4/2 = -3 + 4 cos(t)-cos(2t).
    cosine = [Q(0) for _ in range(ORDER + 1)]
    for j in range(1, ORDER // 2 + 1):
        cosine[2 * j] = Q((-1) ** j * (4 - 2 ** (2 * j)), factorial(2 * j))
    assert mixed == cosine

    # These three entries were obtained geometrically before Note 062.
    independently_known = {3: Q(-1, 2), 4: Q(1, 12), 5: Q(-1, 160)}
    for genus, value in independently_known.items():
        assert mixed[2 * genus - 2] == value

    # Recover the degree-two elementary exponents without inserting them
    # into the extraction: first subtract the degree-one core's degree-two
    # multiple contribution, then eliminate in increasing support genus.
    residual = [a - b for a, b in zip(mixed, elementary)]
    exponents = {}
    for genus in range(3, MAX_GENUS + 1):
        n = residual[2 * genus - 2]
        exponents[genus] = n
        basis = power(s, 2 * genus - 2)
        residual = [a - n * b for a, b in zip(residual, basis)]
    assert not any(residual)
    assert exponents[3] == -8 and exponents[4] == 2
    assert all(exponents[h] == 0 for h in range(5, MAX_GENUS + 1))

    for genus in range(3, 9):
        closed = Q((-1) ** genus * (2 ** (2 * genus - 2) - 4), factorial(2 * genus - 2))
        assert mixed[2 * genus - 2] == closed
        print(f"L_2,{genus} = {closed}")
    print(f"Exact series and exponent checks passed through source genus {MAX_GENUS}.")
    print("This is an arithmetic certificate, not a check of geometric transversality.")


if __name__ == "__main__":
    main()
