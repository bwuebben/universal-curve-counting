"""Exact elementary extraction and independent low-excess tests for Note 077.

Run with uv run --no-project --with sympy==1.14.0 python checks/nonorientable_spin_partition.py.
"""

from fractions import Fraction
from functools import cache
from math import comb

from etale_elementary_exponents import elementary_exponents
from first_excess_spin_hurwitz import first_exponents
from oriented_holonomy_genera import ONE, ZERO, connected, partitions, reciprocal_power
from prime_degree_genus_cutoffs import hooks, multiply, power, sine


@cache
def twisted_disconnected(degree, b):
    result = ZERO
    for shape in partitions(degree, degree):
        transpose = tuple(sum(length > j for length in shape) for j in range(shape[0]))
        if shape == transpose:
            weight = ONE
            for hook in hooks(shape):
                weight *= (2 * ONE - reciprocal_power(hook))**b
            result += weight
    return result


def extract(b, epsilon, limit):
    logs, rows = {}, {}
    for m in range(1, limit + 1):
        disconnected = epsilon**m * twisted_disconnected(m, b)
        lower = sum((i * logs[i] * epsilon**(m-i) * twisted_disconnected(m-i, b)
                     for i in range(1, m)), ZERO)
        logs[m] = disconnected - lower.mul_ground(Fraction(1, m))
        row = logs[m]
        for d, previous in rows.items():
            if m % d == 0:
                for (j,), value in previous.terms():
                    row -= value * connected(m//d, j, 0)
        rows[m] = row
    return logs, rows


def check(genus, epsilon, limit):
    b = genus - 1
    logs, rows = extract(b, epsilon, limit)
    etale = elementary_exponents(genus, limit, True, epsilon)
    first = first_exponents(genus, limit, True, epsilon)
    for m, row in rows.items():
        assert all(n.q == 1 for n in row.all_coeffs())
        assert row.nth(m*b) == etale[m] and row.nth(m*b+1) == first[m]
        assert all(row.nth(j) == 0 for j in range(m*b))
        if m > 1:
            highest = b * comb(m + 1, 2)
            assert row.degree() == highest
            assert row.LC() == -2 * epsilon * (-1)**(m*b + highest)
    if (genus, epsilon) == (2, 1):
        assert rows[3].as_dict() == {(3,): -64, (4,): 64, (5,): -20, (6,): 2}
        # Evaluate the cubic logarithm two ways, including a direct hook Taylor series.
        z_series = power(sine(1), 2)
        hook = multiply(power(sine(3), 2), power(z_series, 2))
        direct = [a + d/3 for a, d in zip(hook, power(z_series, 3))]
        by_polynomial = [Fraction(0) for _ in direct]
        for (j,), n in logs[3].terms():
            by_polynomial = [a + Fraction(n) * term
                             for a, term in zip(by_polynomial, power(z_series, j))]
        assert by_polynomial == direct
        assert direct[6:13:2] == [Fraction(28, 3), Fraction(-25, 3),
                                   Fraction(589, 180), Fraction(-61, 81)]
    print(f"g={genus}, sign={epsilon:+}: exact extraction, highest genus, "
          f"and both low-excess checks through degree {limit}")


if __name__ == '__main__':
    for genus in (2, 3, 4):
        for epsilon in (1, -1):
            check(genus, epsilon, 8)
