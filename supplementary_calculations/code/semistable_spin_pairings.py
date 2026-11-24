#!/usr/bin/env python3
"""Finite checks for Note 076. The Fredholm and virtual arguments are geometric."""

from fractions import Fraction
from itertools import combinations_with_replacement, permutations, product
from math import gcd

from all_genus_double_covers import multiply, power, sine_series
from etale_spin_characters import compose, permutation_sign
from relative_first_excess_collisions import connected


def has_vertex_signs(edges):
    """Encode sign changes at vertices by bits."""
    return any(
        all(bits[i] ^ bits[j] == color for i, j, color in edges)
        for bits in product(range(2), repeat=4)
    )


def all_cycle_products_positive(edges):
    """Even-degree edge subsets form the binary cycle space, including bigons."""
    for mask in range(1, 1 << len(edges)):
        odd_vertices = parity = 0
        for position, (i, j, color) in enumerate(edges):
            if mask & (1 << position):
                odd_vertices ^= (1 << i) | (1 << j)
                parity ^= color
        if not odd_vertices and parity:
            return False
    return True


def graph_checks():
    types = tuple(product(range(2), range(2, 4), range(2), range(2)))
    tested = 0
    for length in range(1, 5):
        for graph in combinations_with_replacement(types, length):
            for inactive in (set(), {2}, {2, 3}):
                even = [(i, j, color) for i, j, color, odd in graph
                        if not odd and j not in inactive]
                assert has_vertex_signs(even) == all_cycle_products_positive(even)
                tested += 1
    assert not has_vertex_signs([(0, 2, 0), (0, 2, 1)])
    return tested


def contact_checks():
    for k in range(1, 33):
        f = 2 // gcd(2, k)
        c = k // gcd(2, k)
        assert Fraction(c, f) == Fraction(k, 2)
        for r in (c, 2*c, 5*c):
            source_order = 2*r // k
            assert source_order == f*r // c
            assert (source_order*k) == 2*r
            assert (1 if k % 2 == 0 else 2) == f
            assert source_order // f == r // c
            assert Fraction(2*c, f) == k
    for k in (1, 2, 3):
        r = k // gcd(2, k)
        assert r * Fraction(k*k, r) == k*k

    for j, odd, trivial in product(range(6), range(0, 12, 2), (False, True)):
        if trivial and odd:
            continue
        sign_genus = j - 1 + odd//2 + int(trivial)
        if not trivial and sign_genus < 0:
            continue  # A connected double torsor cannot have negative Hodge rank.
        assert (sign_genus == 0) == (
            (j, odd, trivial) in {(0, 0, True), (0, 2, False), (1, 0, False)}
        )


def double_cover_checks():
    for genus in range(2, 8):
        chars = range(1, 1 << (2*genus))
        repeated = Fraction(len(chars), 2)
        mixed = Fraction(sum((-1)**(a & 1) for a in chars), 2)
        root_count = 1 << (2*genus-2)
        assert mixed - repeated == -2*root_count
        assert Fraction(-1, 2) * root_count * 4 == mixed - repeated

    # Evaluate (6) independently using the exponent difference (-8,2).
    s = sine_series()
    s4, s6 = power(s, 4), power(s, 6)
    s2prod = multiply(power(s, 2), power(sine_series(2), 2))
    difference = [-2*z for z in s2prod]
    assert difference == [-8*x + 2*y for x, y in zip(s4, s6)]
    assert difference[4:7:2] == [-8, Fraction(10, 3)]
    assert [-2*v for v in difference[4:7:2]] == [16, Fraction(-20, 3)]


def cubic_minimum_genus():
    """Independently count smooth genus-two covers and the twice-punctured elliptic cut."""
    group = list(permutations(range(3)))
    identity = group[0]
    inverse = {a: tuple(a.index(i) for i in identity) for a in group}
    comm = {(a, b): compose(compose(a, b), compose(inverse[a], inverse[b]))
            for a, b in product(group, repeat=2)}
    repeated = mixed = 0
    for a, b, c, d in product(group, repeat=4):
        if compose(comm[a, b], comm[c, d]) == identity and connected(3, (a, b, c, d)):
            repeated += 1
            mixed += permutation_sign(a)
    transpositions = [a for a in group if permutation_sign(a) == -1]
    cut = 0
    for a, b, p, q in product(group, group, transpositions, transpositions):
        if compose(comm[a, b], compose(p, q)) == identity and connected(3, (a, b, p, q)):
            cut += 1
    masses = tuple(Fraction(x, len(group)) for x in (mixed, repeated, cut))
    assert masses == (Fraction(28, 3), Fraction(220, 3), 16)
    assert mixed - repeated == -4 * cut
    return masses


if __name__ == '__main__':
    samples = graph_checks()
    contact_checks()
    double_cover_checks()
    print(f"Passed {samples} colored graph and torsor-pruning tests.")
    print("Checked contact weights through 32, Hodge ranks, and character sums in genera 2–7.")
    print("Degree-two square integrals 16, -20/3; series checked through source genus 24.")
    mixed, repeated, cut = cubic_minimum_genus()
    print(f"Cubic minimum-genus check: {mixed} − {repeated} = −4 × {cut}.")
