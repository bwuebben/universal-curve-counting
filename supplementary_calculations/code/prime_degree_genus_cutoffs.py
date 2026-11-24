#!/usr/bin/env python3
"""Exact cubic arithmetic accompanying Note 063.

Construct the connected elementary coefficient from partitions and hook
lengths, then check the polynomial with two undetermined integer exponents.
No candidate values for these exponents are assumed. This does not verify
projected-jet transversality, cluster transport, or geometric genus bounds.
Only the Python standard library is required.
"""

from fractions import Fraction as Q
from math import factorial


MAX_GENUS = 16
ORDER = 2 * MAX_GENUS - 2


def constant(value):
    return [Q(value)] + [Q(0)] * ORDER


def add(*terms):
    return [sum(values, Q(0)) for values in zip(*terms)]


def scale(series, value):
    return [Q(value) * entry for entry in series]


def multiply(left, right):
    result = constant(0)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right[: ORDER + 1 - i]):
                if b:
                    result[i + j] += a * b
    return result


def power(series, exponent):
    result = constant(1)
    for _ in range(exponent):
        result = multiply(result, series)
    return result


def sine(scale_factor):
    result = constant(0)
    for j in range((ORDER + 1) // 2):
        n = 2 * j + 1
        result[n] = Q((-1) ** j * scale_factor**n, 2 ** (n - 1) * factorial(n))
    return result


def partitions(size, largest=None):
    if size == 0:
        yield ()
        return
    if largest is None:
        largest = size
    for head in range(min(size, largest), 0, -1):
        for tail in partitions(size - head, head):
            yield (head,) + tail


def hooks(partition):
    for row, length in enumerate(partition):
        for column in range(length):
            below = sum(other > column for other in partition[row + 1 :])
            yield length - column + below


def disconnected(degree):
    result = constant(0)
    for partition in partitions(degree):
        weight = constant(1)
        for hook in hooks(partition):
            weight = multiply(weight, power(sine(hook), 2))
        result = add(result, weight)
    return result


def main():
    z = power(sine(1), 2)
    zp = {exponent: power(z, exponent) for exponent in range(MAX_GENUS)}
    z1, z2, z3 = (disconnected(degree) for degree in (1, 2, 3))
    elementary = add(z3, scale(multiply(z1, z2), -1), scale(power(z1, 3), Q(1, 3)))
    expected = add(scale(zp[3], Q(220, 3)), scale(zp[4], -70), scale(zp[5], 21), scale(zp[6], -2))
    assert elementary == expected

    base = add(elementary, scale(zp[3], -64), scale(zp[4], 64))
    assert base[6] == Q(28, 3)
    assert base[8] == Q(-25, 3)
    assert base[10] == Q(4189, 180)
    assert base[12] == Q(-898, 81)
    assert zp[5][10] == 1 and zp[5][12] == Q(-5, 12)
    assert zp[6][10] == 0 and zp[6][12] == 1

    # Several different pairs must remain arithmetically admissible. These
    # checks illustrate why neither integrality nor the cutoff fixes A,B.
    for a, b in ((0, 0), (-20, 2), (17, -9)):
        mixed = add(base, scale(zp[5], a), scale(zp[6], b))
        assert mixed[10] == Q(4189, 180) + a
        assert mixed[12] == Q(-898, 81) - Q(5, 12) * a + b
        residual = add(mixed, scale(elementary, -1))
        extracted = {}
        for genus in range(4, MAX_GENUS + 1):
            coefficient = residual[2 * genus - 2]
            extracted[genus] = coefficient
            residual = add(residual, scale(zp[genus - 1], -coefficient))
        assert extracted[4] == -64 and extracted[5] == 64
        assert extracted[6] == a and extracted[7] == b
        assert all(extracted[h] == 0 for h in range(8, MAX_GENUS + 1))
        assert not any(residual)

    print("Connected cubic hook formula and polynomial agree exactly.")
    print("L_3,6 = 4189/180 + A; L_3,7 = -898/81 - 5*A/12 + B.")
    print(f"Formal extraction checked through source genus {MAX_GENUS} for three distinct pairs (A,B).")
    print("A and B remain undetermined; this is arithmetic, not a geometric proof.")


if __name__ == "__main__":
    main()
