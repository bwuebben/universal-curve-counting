"""Verify the polynomial residues of the balanced scalar cubic handle powers."""

from argparse import ArgumentParser
from math import lcm

from sympy import Poly, QQ, symbols


def unit_part(n):
    """Return the exponent of 3 and the integral unit in positive n."""
    valuation = 0
    while n % 3 == 0:
        valuation += 1
        n //= 3
    return valuation, n


def check_row(n, trace):
    denominator = lcm(*(int(c.q) for c in trace.all_coeffs()))
    k = n // 2
    valuation, unit = unit_part(n)
    exponent = max(0, k - 2 - valuation) if n % 2 else k
    assert denominator == 3**exponent, (n, denominator, exponent)
    if exponent == 0:
        return 0

    reduced = {}
    for (power,), coefficient in trace.terms():
        digit = int(coefficient * denominator) % 3
        if digit:
            reduced[power + 3*n] = digit
    sign = (-1)**k
    expected = {9*k: 2*sign % 3}
    if n % 2:
        digit = 2*unit*sign
        expected = {9*k+2: -digit % 3, 9*k+3: digit % 3}
    assert reduced == expected, (n, reduced, expected)
    return 1


def verify(max_power):
    w = symbols("w")

    def polynomial(expression):
        return Poly(expression, w, domain=QQ)

    factor = w**2 + 54*w - 243
    a = polynomial(w**2 + 50*w + 81)
    b = polynomial(4*(w-6)*factor)
    c = polynomial(4*(w-4)*(w-3)**2*factor)
    # Newton sums for the three eigenvalues after dividing each by w^3.
    t0, t1, t2 = polynomial(3), a, a*a - 2*b/3
    count = 0
    for n in range(1, max_power + 1):
        if n == 1:
            trace = t1
        elif n == 2:
            trace = t2
        else:
            t0, t1, t2 = t1, t2, a*t2 - b*t1/3 + c*t0/3
            trace = t2
        count += check_row(n, trace)
    return count


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--max-power", type=int, default=80)
    args = parser.parse_args()
    if args.max_power < 1:
        parser.error("--max-power must be positive")
    nonintegral = verify(args.max_power)
    print(f"Checked powers 1 through {args.max_power}, including "
          f"the full leading residues of {nonintegral} nonintegral rows.")
