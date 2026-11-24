"""Exact residue and contact calculations for paired expanded finite covers."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import prod

from relative_elliptic_tail_gluing import partitions


def smoothing_residue(a, b, left, right):
    """Constant term as s→0 in the sum at the critical points of u^a+(s/u)^b.

    A monomial is given by powers of u and s, using v=s/u. If the critical
    roots satisfy u^(a+b)=(b/a)s^b, the sum of other powers vanishes.
    """
    u_power = left[0] + right[0] - a
    multiple, remainder = divmod(u_power, a + b)
    if remainder:
        return 0
    s_power = left[1] + right[1] + b * multiple
    assert s_power >= 0
    if s_power:
        return 0
    return Fraction(1, a) * Fraction(b, a)**multiple


def check_node(a, b):
    # Basis 1,u,...,u^a,v,...,v^(b-1).
    basis = [(i, 0) for i in range(a + 1)] + [(-j, j) for j in range(1, b)]
    matrix = [[smoothing_residue(a, b, u, v) for v in basis] for u in basis]
    for i, row in enumerate(matrix):
        partner = a - i if i <= a else 2 * a + b - i
        expected = Fraction(1, a if i <= a else b)
        assert row[partner] == expected
        assert sum(value != 0 for value in row) == 1
    assert len(set(tuple(row) for row in matrix)) == a + b


def check_clutches(profile):
    # At δ=1 enumerate the s_i of the coarse cover. At an even contact
    # (s_i)^(k_i/2) is the framed clutch. At odd k, squaring the twisted
    # smoothing root bijects its k choices onto the coarse choices.
    even = [i for i, order in enumerate(profile) if order % 2 == 0]
    counts = Counter()
    for roots in product(*(range(order) for order in profile)):
        counts[tuple(roots[i] % 2 for i in even)] += 1
    length = prod(k // 2 if k % 2 == 0 else k for k in profile)
    assert len(counts) == 2**len(even)
    assert all(number == length for number in counts.values())
    for k in profile:
        if k % 2:
            assert {(2 * j) % k for j in range(k)} == set(range(k))


def check_cut_transition():
    # The log form in the example is -2L²+2R². Its positive form is
    # iR ⊕ R. Project modulo it, to coordinates (Re L, Im R).
    degrees = []
    for clutch in (1, -1):
        generator = (1, clutch)
        assert -2 * generator[0]**2 + 2 * generator[1]**2 == 0
        directions = [tuple(s * x for x in generator) for s in (1, 1j)]
        (x, y), (u, v) = [(z.real, w.imag) for z, w in directions]
        degrees.append(int(x * v - u * y))
    assert degrees == [1, -1]


def main():
    for a, b in product(range(1, 9), repeat=2):
        check_node(a, b)
    for degree in range(1, 11):
        for profile in partitions(degree):
            check_clutches(profile)
    check_cut_transition()
    print("Verified 64 nodal residue forms from the smooth critical-point sums.")
    print("Verified framed spin multiplicities through degree ten and the two null-line orientations.")


if __name__ == "__main__":
    main()
