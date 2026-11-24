"""Extract balanced exponents from the full partition polynomial over Q[z].

Run with uv run --no-project --with sympy==1.14.0 python checks/oriented_holonomy_genera.py.
"""
from functools import cache
from fractions import Fraction
from itertools import product

from sympy import Poly, QQ, symbols

z = symbols('z')
ONE = Poly(1, z, domain=QQ)
ZERO = 0 * ONE


def partitions(size, bound):
    if size == 0:
        yield ()
    for first in range(min(size, bound), 0, -1):
        for rest in partitions(size - first, first):
            yield (first, *rest)


@cache
def reciprocal_power(power):
    """x^power + x^(-power), expressed in z = 2-x-x^(-1)."""
    if power == 0:
        return 2*ONE
    if power == 1:
        return Poly(2-z, z, domain=QQ)
    return reciprocal_power(1)*reciprocal_power(power-1) - reciprocal_power(power-2)


@cache
def disconnected(degree, b, c):
    result = ZERO
    for shape in partitions(degree, degree):
        boxes = {(i, j) for i, length in enumerate(shape) for j in range(length)}
        content = sum(j-i for i, j in boxes)
        # Average a partition and its transpose before summing all partitions.
        term = reciprocal_power(abs(c*content)).mul_ground(Fraction((-1)**(degree*c), 2))
        for i, j in boxes:
            hook = shape[i]-j + sum(j < length for length in shape[i+1:])
            term *= (2*ONE-reciprocal_power(hook))**b
        result += term
    return result


@cache
def connected(degree, b, c):
    # q Z' = Z q(log Z)'.
    correction = sum((i*connected(i,b,c)*disconnected(degree-i,b,c)
                      for i in range(1, degree)), ZERO)
    return disconnected(degree,b,c)-correction.mul_ground(Fraction(1,degree))


def exponents(b, c, cutoff):
    rows = {}
    for m in range(1, cutoff+1):
        row = connected(m,b,c)
        for d, lower in rows.items():
            if m % d == 0:
                for (power,), coeff in lower.terms():
                    if coeff:
                        row -= coeff * connected(m//d,power,0)
        rows[m] = row
    return rows


if __name__ == '__main__':
    for b,c in product((1,2),(0,1,2)):
        rows = exponents(b,c,6)
        assert rows[1] == Poly((-1)**c*z**b, z, domain=QQ)
        for m,row in rows.items():
            if c == 0 and m > 1:
                assert row.is_zero
                continue
            expected = m*b+(b+c)*m*(m-1)//2
            assert row.degree() == expected
            assert row.LC() == (-1)**(m*(b+c)+expected)
            assert all(coeff.q == 1 for coeff in row.all_coeffs())
        if (b,c) == (1,1):
            assert rows[2] == Poly(15*z**2 - 8*z**3 + z**4,z,domain=QQ)
        print(f'PASS: g={b+1}, c={c}, complete extraction in degrees 1 through 6')
