"""Exact finite-cover and sign checks for transfer at an étale S3 cubic wall.

The surface generators are ordered a,b,c,d with [a,b][c,d]=1. A cover
is a transitive permutation action. An equivariant map to another action
is a lift over the base. Determinant signs are inputs to the calculation.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product

from etale_spin_characters import compose


def inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def commutator(a, b):
    return compose(compose(a, b), compose(inverse(a), inverse(b)))


def surface_relation(action):
    a, b, c, d = action
    return commutator(a, b) == inverse(commutator(c, d))


def check_subgroups(group):
    """Fixed sheets, dim of the fixed standard plane, and the Euler multiplier."""
    rows = Counter()
    for size in (1, 2, 3, 6):
        for tail in combinations(group[1:], size - 1):
            subgroup = {group[0], *tail}
            if any(compose(a, b) not in subgroup for a in subgroup for b in subgroup):
                continue
            cubic = sum(all(a[i] == i for a in subgroup) for i in range(3))
            regular = sum(all(compose(a, s) == s for a in subgroup) for s in group)
            orbits = {frozenset(a[i] for a in subgroup) for i in range(3)}
            rank = len(orbits) - 1  # subtract the constant line in R^3
            # An étale pullback of index zero changes determinant sign by (-1)^rank.
            assert 2 * cubic - regular == 1 - (-1) ** rank
            rows[(size, cubic, regular, rank)] += 1
    assert rows == {(1, 3, 6, 2): 1, (2, 1, 0, 1): 3,
                    (3, 0, 0, 0): 1, (6, 0, 0, 0): 1}
    # Generic marked evaluation V -> R^4: ranks 1,0 have codimension
    # (2-r)(4-r), greater than the two real parameters of a regular mark.
    assert [(2 - r) * (4 - r) for r in (1, 0)] == [3, 8]


def action_key(action, colors=()):
    """Canonicalize by labeling from each possible root in generator order."""
    keys = []
    for root in range(len(action[0])):
        sheets = [root]
        label = {root: 0}
        for sheet in sheets:
            for generator in action:
                target = generator[sheet]
                if target not in label:
                    label[target] = len(sheets)
                    sheets.append(target)
        if len(sheets) != len(action[0]):
            return None
        transitions = tuple(label[a[i]] for i in sheets for a in action)
        keys.append((transitions, tuple(colors[i] for i in sheets) if colors else ()))
    return min(keys)


def lift_count(source, target):
    """Enumerate all equivariant maps by propagating the image of sheet zero."""
    count = 0
    for initial in range(len(target[0])):
        images = {0: initial}
        todo = [0]
        for i in todo:
            for a, b in zip(source, target):
                if a[i] not in images:
                    images[a[i]] = b[images[i]]
                    todo.append(a[i])
        assert len(images) == len(source[0])
        if all(images[a[i]] == b[j] for i, j in images.items() for a, b in zip(source, target)):
            count += 1
    return count


def check_genus_two(group):
    swap, cycle = (1, 0, 2), (1, 2, 0)
    cubic = (swap, cycle, cycle, swap)
    closure = tuple(tuple(group.index(compose(a, s)) for s in group) for a in cubic)
    assert surface_relation(cubic) and surface_relation(closure)
    cubic_actions = {}
    for action in product(group, repeat=4):
        if surface_relation(action):
            key = action_key(action)
            if key is not None:
                cubic_actions[key] = action
    cubic_lifts = Counter((lift_count(u, cubic), lift_count(u, closure))
                          for u in cubic_actions.values())
    assert cubic_lifts == {(1, 0): 1, (0, 0): 99}
    # 40 cyclic cubics (nonzero mod-3 characters up to sign), and 60 S3 cubics.
    assert Counter(lift_count(u, u) for u in cubic_actions.values()) == {3: 40, 1: 60}

    # Every two-sheeted cover of Y is described over C by 12 fiber flips,
    # one per generator per sheet of Y. Sheet colors record the map to Y.
    def double_lifts(a):
        return [tuple(2 * a[i] + (bit ^ flips[i]) for i in range(3) for bit in (0, 1))
                for flips in product((0, 1), repeat=3)]

    upstairs = {}
    for action in product(*(double_lifts(a) for a in cubic)):
        if surface_relation(action):
            key = action_key(action, (0, 0, 1, 1, 2, 2))
            if key is not None:
                upstairs[key] = action
    assert len(upstairs) == 2 ** 8 - 1  # g(Y)=4; nontrivial characters of H_1(Y; F_2)
    downstairs = {action_key(action): action for action in upstairs.values()}
    assert len(downstairs) == len(upstairs) and action_key(closure) in downstairs
    by_lifts = Counter((lift_count(u, cubic), lift_count(u, closure), lift_count(u, u))
                       for u in downstairs.values())
    assert by_lifts == {(1, 0, 2): 254, (3, 6, 6): 1}
    return cubic_lifts, by_lifts


def check_signs(h, negatives):
    """Test both crossing orientations and both possible closure signs."""
    total = 2 ** (2 * h) - 1
    for n_minus, eta, closure_sign in product(negatives, (-1, 1), (-1, 1)):
        assert 0 <= n_minus <= total - 1
        signed_regular = eta * (total - 1 - 2 * n_minus)
        # Each of the total upstairs doubles has deck order two. Every
        # ordinary composite has deck order two and a rank-one crossing;
        # the closure has deck order six and an even-dimensional crossing.
        upper_double_count = Fraction(signed_regular + closure_sign, 2)
        transfer = 2 * upper_double_count - closure_sign
        assert transfer == 2 * Fraction(signed_regular, 2)
        # Subtract the leading double-cover coefficient of the elementary
        # genus-h curve times the cubic primitive jump 2 eta.
        primitive = transfer - (2 * eta) * Fraction(total, 2)
        assert primitive == -eta * (1 + 2 * n_minus)
        assert primitive.denominator == 1 and primitive.numerator % 2 == 1


if __name__ == "__main__":
    s3 = tuple(permutations(range(3)))
    check_subgroups(s3)
    cubics, sextics = check_genus_two(s3)
    check_signs(4, range(255))
    for g in (3, 4, 5):
        genus = 3 * g - 2
        ordinary = 2 ** (2 * genus) - 2
        check_signs(genus, (0, 1, ordinary // 2, ordinary - 1, ordinary))
    print("S3 fixed-sheet identity verified for all six subgroups.")
    print(f"Genus-two cubics by (Y lifts, Z lifts): {dict(sorted(cubics.items()))}")
    print(f"Composed sextics by (Y lifts, Z lifts, deck order): {dict(sorted(sextics.items()))}")
    print("Degree-six primitive change equals -eta (1 + 2 N_minus) for all genus-two sign totals.")
