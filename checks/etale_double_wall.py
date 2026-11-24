"""Exact partition extraction and genus-two quartic monodromy checks.

Run with uv run --no-project --with sympy==1.14.0 python checks/etale_double_wall.py.
"""
from collections import Counter, defaultdict
from functools import cache
from fractions import Fraction
from itertools import permutations, product

from sympy import Poly, QQ

from etale_elementary_exponents import elementary_exponents
from etale_spin_characters import compose, permutation_sign
from oriented_holonomy_genera import ONE, connected, exponents, reciprocal_power, z


def check_wall(genus, upstairs_cutoff):
    b = 2 * (genus - 1)
    rows = exponents(b, b, upstairs_cutoff)
    expected = Poly(z**(2*b) * (4-z)**b, z, domain=QQ) * (reciprocal_power(b) - 2*ONE)
    assert rows[2] == expected
    assert connected(2, b, b).nth(2*b) == Fraction(2**(2*(b+1))-1, 2)
    assert expected.nth(2*b+1) == -b*b*4**b
    if genus == 2:
        assert expected.as_dict() == {(5,): -64, (6,): 48, (7,): -12, (8,): 1}
    for e, row in rows.items():
        assert row.degree() == b * e**2 and row.LC() == 1
        assert all(coefficient.q == 1 for coefficient in row.all_coeffs())
    print(f"g={genus}: partition extraction through projected degree {2*upstairs_cutoff}; "
          f"quartic least-genus connected jump {connected(2, b, b).nth(2*b)}")


def count_quartic_twists():
    """Enumerate genus-two surface-group representations on four letters up to relabelling."""
    group = list(permutations(range(4)))
    index = {element: i for i, element in enumerate(group)}
    identity = index[tuple(range(4))]
    table = [[index[compose(a, b)] for b in group] for a in group]
    inverse = [row.index(identity) for row in table]
    by_commutator = defaultdict(list)
    for a, b in product(range(24), repeat=2):
        commutator = table[table[table[a][b]][inverse[a]]][inverse[b]]
        by_commutator[commutator].append((a, b))

    @cache
    def generated_subgroup(generators):
        found = {identity}
        pending = [identity]
        while pending:
            x = pending.pop()
            for gen in generators:
                y = table[x][gen]
                if y not in found:
                    found.add(y)
                    pending.append(y)
        return frozenset(found)

    @cache
    def dihedral_rotations(subgroup):
        # The square symmetries have a unique cyclic subgroup of order four.
        r = next(x for x in subgroup
                 if table[x][x] != identity and table[table[x][x]][table[x][x]] == identity)
        return generated_subgroup((r,))

    labelled = Counter()
    for commutator, first_pairs in by_commutator.items():
        for (a, b), (c, d) in product(first_pairs, by_commutator[inverse[commutator]]):
            subgroup = generated_subgroup(tuple(sorted({a, b, c, d})))
            # The difference of the chosen spin structures pairs by evaluation on a.
            if len(subgroup) == 24 and permutation_sign(group[a]) == -1:
                labelled["symmetric"] += 1
            if len(subgroup) == 8 and a not in dihedral_rotations(subgroup):
                labelled["dihedral"] += 1
    # Centralizer orders are one and two. Each dihedral quotient has two
    # inequivalent four-letter actions (by reflections of the two classes).
    assert labelled["symmetric"] == 24 * 480
    assert labelled["dihedral"] == 12 * 2 * 96
    print("g=2: 480 sign-reversing S4 quartics and 96 pairs of D4 quartics")


def check_spin_exponent(genus):
    c = genus - 1
    negative = 2**(2*genus-1)
    primitive = 2**(2*c) * (3**(2*c)-1) * (2**(4*c)-1)
    predicted = negative*(negative-1)//2 - 2*primitive
    row = elementary_exponents(genus, 4, twisted=True, epsilon=1)
    assert row[2] == -negative and row[4] == predicted
    print(f"g={genus}: quartic spin exponent {predicted}, primitive S4 covers {primitive}")


if __name__ == "__main__":
    for genus in (2, 3):
        check_wall(genus, upstairs_cutoff=5)
    count_quartic_twists()
    for genus in (2, 3, 4, 5):
        check_spin_exponent(genus)
