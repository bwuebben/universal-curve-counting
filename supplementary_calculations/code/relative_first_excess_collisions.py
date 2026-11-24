"""Compute first-excess cut/glue masses from labeled cover monodromies."""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, prod

from etale_spin_characters import compose
from first_excess_spin_hurwitz import spin_series


def cycles(element):
    todo = set(range(len(element)))
    parts = []
    while todo:
        point = min(todo)
        part = []
        while point in todo:
            part.append(point)
            todo.remove(point)
            point = element[point]
        parts.append(tuple(part))
    return parts


def connected(degree, generators, pair=()):
    orbit = {0}
    while True:
        enlarged = orbit.union(*(set(gen[i] for i in orbit) for gen in generators))
        if any(s in orbit for s in pair):
            enlarged.update(pair)
        if enlarged == orbit:
            return len(orbit) == degree
        orbit = enlarged


def genus_two_cut(degree):
    group = list(permutations(range(degree)))
    identity = group[0]
    inverse = {g: tuple(g.index(i) for i in identity) for g in group}
    triples = [g for g in group if sum(a != b for a, b in zip(g, identity)) == 3]
    pairs = list(combinations(identity, 2))
    # h identifies the paired fibers. Its phase on each cycle is a smoothing root.
    phases = defaultdict(lambda: defaultdict(list))
    for p in group:
        for h in group:
            q = compose(h, compose(inverse[p], inverse[h]))
            matching = tuple(frozenset(h[i] for i in part) for part in cycles(p))
            phases[p, q][matching].append(h)

    cut_degrees = Counter()
    masses = Counter()
    for a in group:
        for b in group:
            comm = compose(a, compose(b, compose(inverse[a], inverse[b])))
            for p in group:
                orders = tuple(sorted(map(len, cycles(p))))
                lhs = compose(comm, p)
                for branch in [identity, *triples]:
                    # The cut surface has relation [a,b] p q branch = 1.
                    q = compose(inverse[lhs], inverse[branch])
                    for roots in phases[p, q].values():
                        assert len(roots) == prod(orders)
                        # Changing h within one matching composes with rotations of
                        # p-cycles. These preserve every orbit of the seam group.
                        generators = (a, b, p, q, branch, roots[0])
                        if branch == identity:
                            choices = sum(connected(degree, generators, pair) for pair in pairs)
                            kind = 'nodal'
                        else:
                            choices = int(connected(degree, generators))
                            kind = 'triple'
                        degree_of_choice = Fraction(choices, factorial(degree))
                        cut_degrees[kind, orders] += degree_of_choice
                        masses[kind, orders] += degree_of_choice * len(roots)
    return cut_degrees, masses


def main():
    for degree in range(2, 5):
        cut_degrees, masses = genus_two_cut(degree)
        for (kind, profile), count in cut_degrees.items():
            assert masses[kind, profile] == prod(profile) * count
        for mixed in (False, True):
            nodal, triple, _, _ = spin_series(2, degree, mixed)
            for kind, smooth in [('nodal', nodal), ('triple', triple)]:
                signed = sum(value * ((-1) ** (degree - len(profile)) if mixed else 1)
                             for (sector, profile), value in masses.items() if sector == kind)
                assert signed == smooth[degree]
        if degree == 2:
            assert masses['nodal', (1, 1)] == masses['nodal', (2,)] == 4
            # Each generic double-cover deck group has order two.
            assert all(2 * masses['nodal', profile] == 8 for profile in [(1, 1), (2,)])
        if degree == 3:
            profiles = [(1, 1, 1), (1, 2), (3,)]
            assert [cut_degrees['nodal', p] for p in profiles] == [50, 52, 27]
            assert [cut_degrees['triple', p] for p in profiles] == [18, 36, 15]
        print(f"Genus 2, degree {degree}: verified ordinary and signed nodal and triple masses.")

    # All 2^(2g) normalizations, including the trivial disconnected one.
    for genus in range(2, 8):
        profiles = Counter()
        for character in range(1 << (2 * genus)):
            seam = character & 1
            profiles[seam] += 1
        assert profiles[0] == profiles[1] == 2**(2 * genus - 1)
        per_root = -Fraction(genus - 1, 8)
        assert sum((-1)**seam * n * per_root for seam, n in profiles.items()) == 0
        assert sum(profiles.values()) * per_root == -(genus - 1) * 2**(2 * genus - 3)
    print("Cut masses recover both mixed and equal-spin smooth formulas; genus-two controls are 0 and -2.")


if __name__ == "__main__":
    main()
