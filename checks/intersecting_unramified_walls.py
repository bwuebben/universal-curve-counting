#!/usr/bin/env python3
"""Finite orbit and determinant checks for intersecting unramified walls.

This is independent of the other project checks and uses only the standard
library. Sphere Euler characteristics come from actual equivariant cells:
two points for a sign character, and the barycentric subdivision of a polygon
for the standard dihedral plane. Products are actual Cartesian products of
these cell sets. Every subgroup, not just subgroup conjugacy representatives,
is checked against the alternating determinant signs at all parameter corners.

This verifies finite representation theory and the minimum-degree deck
normalization. It does not prove analytic realization, compactness, or a wall
formula. Positive-corner determinant signs are inputs, checked as both +/-1.
"""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations, product


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    return tuple(a.index(i) for i in range(len(a)))


class Group:
    def __init__(self, name, elements):
        self.name = name
        self.elements = tuple(sorted(elements))
        self.identity = tuple(range(len(self.elements[0])))
        self.number = {g: i for i, g in enumerate(self.elements)}
        assert self.identity in self.elements
        assert all(compose(a, b) in self.number
                   for a, b in product(self.elements, repeat=2))
        trivial = frozenset((self.identity,))
        seen, pending = {trivial}, [trivial]
        while pending:
            subgroup = pending.pop()
            for g in self.elements:
                if g not in subgroup:
                    candidate = self.generated((*subgroup, g))
                    if candidate not in seen:
                        seen.add(candidate)
                        pending.append(candidate)
        self.subgroups = tuple(sorted(seen, key=self.subgroup_key))

    def generated(self, generators):
        reached, pending = {self.identity}, [self.identity]
        while pending:
            g = pending.pop()
            for h in generators:
                candidate = compose(h, g)
                if candidate not in reached:
                    reached.add(candidate)
                    pending.append(candidate)
        return frozenset(reached)

    def subgroup_key(self, subgroup):
        return len(subgroup), tuple(sorted(self.number[g] for g in subgroup))

    def conjugate(self, subgroup, g):
        return frozenset(compose(compose(g, h), inverse(g)) for h in subgroup)

    def canonical(self, subgroup):
        return min((self.conjugate(subgroup, g) for g in self.elements),
                   key=self.subgroup_key)

    def deck_order(self, subgroup):
        normalizer = sum(self.conjugate(subgroup, g) == subgroup
                         for g in self.elements)
        return normalizer // len(subgroup)


@dataclass
class GSet:
    size: int
    action: dict

    def fixed(self, subgroup):
        return sum(all(self.action[g][i] == i for g in subgroup)
                   for i in range(self.size))


def make_gset(group, action):
    size = len(action[group.identity])
    assert all(sorted(a) == list(range(size)) for a in action.values())
    assert all(action[compose(g, h)] == compose(action[g], action[h])
               for g, h in product(group.elements, repeat=2))
    return GSet(size, action)


def cartesian(group, left, right):
    return make_gset(group, {
        g: tuple(left.action[g][i] * right.size + right.action[g][j]
                 for i in range(left.size) for j in range(right.size))
        for g in group.elements
    })


def orbit_decomposition(group, gset):
    remaining, answer = set(range(gset.size)), Counter()
    while remaining:
        point = min(remaining)
        orbit = {gset.action[g][point] for g in group.elements}
        stabilizer = frozenset(g for g in group.elements
                               if gset.action[g][point] == point)
        assert len(orbit) * len(stabilizer) == len(group.elements)
        answer[group.canonical(stabilizer)] += 1
        remaining.difference_update(orbit)
    return answer


def multiply_cells(group, factors):
    singleton = make_gset(group, {g: (0,) for g in group.elements})
    terms = [(1, singleton)]
    for cells in factors:
        terms = [(a * b, cartesian(group, x, y))
                 for a, x in terms for b, y in cells]
    return terms


def decompose_cells(group, cells):
    answer = Counter()
    for coefficient, gset in cells:
        for subgroup, multiplicity in orbit_decomposition(group, gset).items():
            answer[subgroup] += coefficient * multiplicity
    return {h: value for h, value in answer.items() if value}


def coset_mark(group, subgroup, stabilizer):
    # Count fixed cosets from their representatives, each repeated |H| times.
    numerator = sum(all(compose(compose(inverse(g), k), g) in stabilizer
                        for k in subgroup) for g in group.elements)
    assert numerator % len(stabilizer) == 0
    return numerator // len(stabilizer)


def matrix_product(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def matrix_rank(rows, columns):
    work = [list(map(Fraction, row)) for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(work))
                      if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [x / scale for x in work[rank]]
        for i in range(len(work)):
            if i != rank:
                scale = work[i][column]
                work[i] = [x - scale * y for x, y in zip(work[i], work[rank])]
        rank += 1
    return rank


def fixed_dimension(representation, subgroup):
    dimension = len(next(iter(representation.values())))
    rows = [[representation[g][i][j] - int(i == j)
             for j in range(dimension)]
            for g in subgroup for i in range(dimension)]
    return dimension - matrix_rank(rows, dimension)


def sign_factor(group, character):
    representation = {g: ((character(g),),) for g in group.elements}
    points = make_gset(group, {
        g: (0, 1) if character(g) == 1 else (1, 0)
        for g in group.elements
    })
    return representation, [(1, points)]


def polygon_factor(group, n, representation):
    # Group elements act on the n vertices by j -> a + sign*j.
    data = {}
    for g in group.elements:
        a = g[0]
        step = (g[1] - g[0]) % n
        assert step in (1, n - 1)
        sign = 1 if step == 1 else -1
        assert all(g[j] == (a + sign * j) % n for j in range(n))
        data[g] = a, sign
    vertices = make_gset(group, {
        g: tuple((2 * a + sign * j) % (2 * n) for j in range(2 * n))
        for g, (a, sign) in data.items()
    })
    edges = make_gset(group, {
        g: tuple((2 * a + sign * j - int(sign == -1)) % (2 * n)
                 for j in range(2 * n))
        for g, (a, sign) in data.items()
    })
    # A stabilizer of an open half-edge must fix both endpoints. Without
    # barycentric subdivision an edge setwise stabilizer would not suffice.
    for g, vertex_action in vertices.action.items():
        for j in range(2 * n):
            if edges.action[g][j] == j:
                assert vertex_action[j] == j
                assert vertex_action[(j + 1) % (2 * n)] == (j + 1) % (2 * n)
    return representation, [(1, vertices), (-1, edges)]


def corner_sum(dimensions, positive_sign=1):
    total = 0
    for negative in product((0, 1), repeat=len(dimensions)):
        corner_sign = positive_sign * (-1) ** sum(
            bit * dimension for bit, dimension in zip(negative, dimensions))
        total += (-1) ** sum(negative) * corner_sign
    return total


def verify_case(label, group, factors, expected, expected_degree, expected_jumps):
    for representation, cells in factors:
        assert all(representation[compose(g, h)] ==
                   matrix_product(representation[g], representation[h])
                   for g, h in product(group.elements, repeat=2))
        for subgroup in group.subgroups:
            direct = sum(c * x.fixed(subgroup) for c, x in cells)
            assert direct == 1 - (-1) ** fixed_dimension(representation, subgroup)
    cells = multiply_cells(group, [x[1] for x in factors])
    actual = decompose_cells(group, cells)
    expected = {group.canonical(h): c for h, c in expected.items() if c}
    assert actual == expected, (label, actual, expected)
    marks = {}
    for subgroup in group.subgroups:
        dimensions = [fixed_dimension(r, subgroup) for r, _ in factors]
        direct = sum(c * x.fixed(subgroup) for c, x in cells)
        decomposed = sum(c * coset_mark(group, subgroup, h)
                         for h, c in actual.items())
        assert direct == decomposed == corner_sum(dimensions)
        for positive_sign in (-1, 1):
            assert corner_sum(dimensions, positive_sign) == positive_sign * direct
        marks[subgroup] = direct
    nonzero = [h for h, mark in marks.items() if mark]
    degree = min((len(group.elements) // len(h) for h in nonzero), default=None)
    assert degree == expected_degree
    jumps = {}
    if degree is not None:
        assert all(len(group.elements) // len(h) >= degree for h in actual)
        for h in nonzero:
            if len(group.elements) // len(h) != degree:
                continue
            representative = group.canonical(h)
            deck = group.deck_order(h)
            jump = Fraction(marks[h], deck)
            assert jump == actual[representative]
            dimensions = [fixed_dimension(r, h) for r, _ in factors]
            for positive_sign in (-1, 1):
                weighted_corners = Fraction(corner_sum(dimensions, positive_sign), deck)
                assert weighted_corners == positive_sign * actual[representative]
            jumps[representative] = (jump, deck)
    assert sorted(jump for jump, _ in jumps.values()) == sorted(expected_jumps)
    summary = sorted((len(group.elements) // len(h), c,
                      group.deck_order(h)) for h, c in actual.items())
    print(f"{label}: {len(group.subgroups)} subgroups, "
          f"{2 ** len(factors)} corners; orbit (degree, coefficient, deck) {summary}.")
    if degree is None:
        print("  Mixed element and every mixed mark vanish.")
    else:
        print(f"  First affected degree {degree}; "
              f"positive-corner-normalized jumps {sorted(jumps.values())}.")
    return marks


def elementary_abelian(rank):
    order = 1 << rank
    return Group(f"C2^{rank}",
                 {tuple(i ^ mask for i in range(order)) for mask in range(order)})


def check_characters(label, ambient_rank, masks, rank):
    group = elementary_abelian(ambient_rank)
    characters = [lambda g, mask=mask: (-1) ** ((g[0] & mask).bit_count())
                  for mask in masks]
    factors = [sign_factor(group, character) for character in characters]
    common_kernel = frozenset(g for g in group.elements
                              if all(character(g) == 1 for character in characters))
    degree = len(group.elements) // len(common_kernel)
    assert degree == 2 ** rank
    multiplicity = 2 ** (len(masks) - rank)
    verify_case(label, group, factors, {common_kernel: multiplicity},
                degree, [multiplicity])


def main():
    check_characters("Independent two characters", 2, [1, 2], 2)
    check_characters("Repeated/dependent two characters", 2, [1, 1], 1)
    check_characters("Three distinct characters of rank two", 2, [1, 2, 3], 2)
    check_characters("Four distinct characters of rank three", 3, [1, 2, 4, 7], 3)

    symmetric = Group("S3", set(permutations(range(3))))
    sign = lambda g: (-1) ** sum(g[i] > g[j] for i in range(3) for j in range(i + 1, 3))
    standard = {
        g: tuple(tuple(int(g[j] == i) - int(g[2] == i) for j in range(2))
                 for i in range(2)) for g in symmetric.elements
    }
    verify_case("S3 sign times standard (vanishing control)", symmetric,
                [sign_factor(symmetric, sign),
                 polygon_factor(symmetric, 3, standard)], {}, None, [])

    dihedral = Group("D4 (order 8)",
                     {tuple((a + sign * j) % 4 for j in range(4))
                      for a in range(4) for sign in (1, -1)})
    rotation = (1, 2, 3, 0)
    reflection = (0, 3, 2, 1)
    h_s = dihedral.generated((reflection,))
    h_rs = dihedral.generated((compose(rotation, reflection),))
    assert dihedral.canonical(h_s) != dihedral.canonical(h_rs)
    trivial = frozenset((dihedral.identity,))
    r_matrix = ((0, -1), (1, 0))
    s_matrix = ((1, 0), (0, -1))
    identity_matrix = ((1, 0), (0, 1))
    matrices = {}
    for g in dihedral.elements:
        matrix = identity_matrix
        for _ in range(g[0]):
            matrix = matrix_product(matrix, r_matrix)
        if (g[1] - g[0]) % 4 == 3:
            matrix = matrix_product(matrix, s_matrix)
        matrices[g] = matrix
    alpha = lambda g: (-1) ** g[0]  # alpha(r)=-1, alpha(s)=+1.
    marks = verify_case("D4 alpha times standard", dihedral,
                        [sign_factor(dihedral, alpha),
                         polygon_factor(dihedral, 4, matrices)],
                        {h_s: 2, trivial: -1}, 4, [2])
    assert marks[h_s] == 4
    assert marks[h_rs] == 0
    assert marks[trivial] == 0
    print("  Reflection classes distinguished: marks at <s>, <rs>, e are 4, 0, 0.")
    print("All finite checks passed. No analytic wall realization is asserted.")


if __name__ == "__main__":
    main()
