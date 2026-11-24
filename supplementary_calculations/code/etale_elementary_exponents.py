"""Compute least-genus exponents and independently count cubic covers by monodromy."""

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import permutations, product
from math import factorial

from etale_spin_characters import compose, hook_count, partitions, permutation_sign


@cache
def connected_counts(genus, limit, twisted):
    """Log of the disconnected hook series with positive base sign."""
    disconnected = [1] + [hook_count(m, genus, twisted) for m in range(1, limit + 1)]
    connected = [Fraction(0)]
    for m in range(1, limit + 1):
        # Compare coefficients of B' = B (log B)'.
        previous = sum(k * connected[k] * disconnected[m - k] for k in range(1, m))
        connected.append(disconnected[m] - Fraction(previous, m))
    return connected


def elementary_exponents(genus, limit, twisted, epsilon):
    connected = connected_counts(genus, limit, twisted)
    exponents = [0]
    for m in range(1, limit + 1):
        contribution = epsilon**m * connected[m]
        for d in range(1, m):
            if m % d == 0:
                k = m // d
                source_genus = 1 + d * (genus - 1)
                contribution -= exponents[d] * connected_counts(source_genus, k, False)[k]
        assert contribution.denominator == 1
        exponents.append(int(contribution))
    return exponents


def cubic_orbits(genus):
    """Enumerate actual transitive tuples and their simultaneous conjugacy orbits."""
    group = tuple(permutations(range(3)))
    identity = (0, 1, 2)
    inverse = {a: tuple(a.index(i) for i in identity) for a in group}
    conjugate = {(a, b): compose(compose(a, b), inverse[a]) for a in group for b in group}
    commutator = {(a, b): compose(conjugate[a, b], inverse[b]) for a in group for b in group}
    sign_bit = {a: int(permutation_sign(a) < 0) for a in group}
    seen = set()
    covers = Counter()
    tuple_count = signed_tuple_count = 0
    for generators in product(group, repeat=2 * genus):
        relation = identity
        for a, b in zip(generators[::2], generators[1::2]):
            relation = compose(relation, commutator[a, b])
        if relation != identity:
            continue
        # An intransitive action on three sheets has a fixed point.
        if any(all(a[i] == i for a in generators) for i in identity):
            continue
        tuple_count += 1
        signed_tuple_count += (-1) ** sign_bit[generators[0]]
        if generators in seen:
            continue
        orbit = {tuple(conjugate[s, a] for a in generators) for s in group}
        seen.update(orbit)
        deck_order = sum(all(conjugate[s, a] == a for a in generators) for s in group)
        assert len(orbit) * deck_order == len(group)
        character = tuple(sign_bit[a] for a in generators)
        covers[character, deck_order] += 1
    assert len(seen) == tuple_count
    ordinary = sum(Fraction(n, deck) for (_, deck), n in covers.items())
    twisted = sum(Fraction(n * (-1)**chi[0], deck) for (chi, deck), n in covers.items())
    assert ordinary == Fraction(tuple_count, 6)
    assert twisted == Fraction(signed_tuple_count, 6)
    return covers, tuple_count, signed_tuple_count, ordinary, twisted


def check_cubics(genus, n3):
    covers, labelled, signed_labelled, ordinary, twisted = cubic_orbits(genus)
    zero = (0,) * (2 * genus)
    cyclic = (3 ** (2 * genus) - 1) // 2
    per_character = (3 ** (2 * genus - 2) - 1) // 2
    expected = Counter({(zero, 3): cyclic})
    for chi in product((0, 1), repeat=2 * genus):
        if chi != zero:
            expected[chi, 1] = per_character
    assert covers == expected
    odd = sum(n for (chi, _), n in covers.items() if chi[0])
    assert ordinary == connected_counts(genus, 3, False)[3]
    assert twisted == connected_counts(genus, 3, True)[3]
    assert twisted - ordinary == -2 * odd == n3
    print(f"S3 g={genus}: labelled={labelled}, signed={signed_labelled}, "
          f"ordinary={ordinary}, twisted={twisted}, cyclic={cyclic}, "
          f"per nonzero character={per_character}, odd on a1={odd}")


if __name__ == "__main__":
    limit = 7
    for genus in (2, 3):
        ordinary = elementary_exponents(genus, limit, False, 1)
        assert ordinary == [0, 1] + [0] * (limit - 1)
        for epsilon in (1, -1):
            signed = elementary_exponents(genus, limit, True, epsilon)
            print(f"g={genus}, nontrivial twist, epsilon={epsilon:+d}, n1..n{limit}={signed[1:]}")
            double_covers = 2 ** (2 * genus) - 1
            assert signed[2] == (-1 - epsilon * double_covers) // 2
            for prime in (3, 5, 7):
                magnitude = -epsilon * signed[prime]
                group_order = factorial(prime)
                partition_count = sum(1 for _ in partitions(prime))
                bound = 2 * group_order * (group_order * partition_count) ** (genus - 2)
                assert magnitude >= bound
                assert magnitude % (2 ** (2 * genus)) == 0
                if epsilon == 1:
                    number = magnitude // (2 ** (2 * genus))
                    print(f"  p={prime}: M_p={number}, lower bound on |n_p|={bound}")
            if epsilon == 1:
                check_cubics(genus, signed[3])
