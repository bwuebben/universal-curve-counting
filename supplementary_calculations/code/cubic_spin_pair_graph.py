#!/usr/bin/env python3
"""Check the 45 Szegő quadratic differentials of Note 078 using integer arithmetic."""

from itertools import combinations
from math import prod


def quadratic(first, second, branch_values):
    """Coefficients of P'R−PR', with P and R the two common pairs."""
    if len(first & second) == 1:
        second = branch_values - second
    a, b = sorted(first & second)
    e, f = sorted(branch_values - (first | second))
    linear_p, constant_p = -(a + b), a * b
    linear_r, constant_r = -(e + f), e * f
    leading = linear_r - linear_p
    linear = 2 * (constant_r - constant_p)
    constant = linear_p * constant_r - linear_r * constant_p
    resultant = prod(i - j for i in (a, b) for j in (e, f))
    assert linear**2 - 4 * leading * constant == 4 * resultant != 0
    return leading, linear, constant


def label(triple):
    return "".join(map(str, sorted(triple)))


def main():
    branch_values = frozenset(range(1, 7))
    spins = [frozenset((1, i, j)) for i, j in combinations(range(2, 7), 2)]
    good = set()
    exceptional = {}
    infinity = set()
    for first, second in combinations(spins, 2):
        edge = (label(first), label(second))
        leading, linear, constant = quadratic(first, second, branch_values)
        branch_zeros = sorted(
            x for x in branch_values if (leading * x + linear) * x + constant == 0
        )
        if leading == 0:
            assert linear != 0  # Simple root at the unbranched value infinity.
            infinity.add(edge)
        if branch_zeros:
            exceptional[edge] = branch_zeros
        else:
            good.add(edge)

    assert exceptional == {("126", "135"): [4], ("135", "156"): [3]}
    assert len(infinity) == 7
    assert len(good) == 43
    assert all(("123", label(spin)) in good for spin in spins[1:])
    assert ("124", "125") in good

    print("43 of the 45 spin pairs have four simple zeros.")
    for (first, second), roots in exceptional.items():
        print(f"{first}–{second} has a zero at branch value {roots[0]} and at infinity.")
        assert (first, second) in infinity
    print("Seven projective quadratics have a simple zero at the unbranched infinity.")
    print("The good pairs include the star at 123 and the triangle 123–124–125.")


if __name__ == "__main__":
    main()
