"""Finite dihedral coset actions and the minimum-genus wall subtraction.

An element (a, s) acts on the n-gon by j -> a + s*j, with s = +/-1.
Kernel dimensions use the faithful real plane. Determinant signs are
inputs. The surface-group enumeration for the cubic is in etale_cubic_transfer.py.
"""

from fractions import Fraction
from itertools import product


def multiply(left, right, n):
    a, s = left
    b, t = right
    return ((a + s * b) % n, s * t)


def subgroups(n):
    """K intersects <r> in <r^step>; its reflections, if present, form one coset."""
    for step in range(1, n + 1):
        if n % step:
            continue
        rotations = frozenset((j, 1) for j in range(0, n, step))
        yield rotations
        for first in range(step):
            reflections = {(j, -1) for j in range(first, n, step)}
            yield rotations | reflections


def left_cosets(group, stabilizer, n):
    """A representative and the complete fiber of each point of G/H."""
    seen = set()
    cosets = []
    for g in group:
        if g not in seen:
            fiber = {multiply(g, h, n) for h in stabilizer}
            cosets.append((g, fiber))
            seen.update(fiber)
    assert len(cosets) * len(stabilizer) == len(group)
    return cosets


def fixed_points(cosets, subgroup, n):
    return sum(all(multiply(k, g, n) in fiber for k in subgroup) for g, fiber in cosets)


def fixed_plane_dimension(subgroup):
    # Every nonidentity rotation in the faithful plane fixes only zero.
    if any(a != 0 and s == 1 for a, s in subgroup):
        return 0
    # Two distinct reflections would generate a nonidentity rotation.
    return 2 if len(subgroup) == 1 else 1


def check_group(n):
    group = tuple(product(range(n), (1, -1)))
    identity = (0, 1)
    representatives = (0,) if n % 2 else (0, 1)
    reflections = [frozenset((identity, (a, -1))) for a in representatives]
    quotients = [left_cosets(group, h, n) for h in reflections]
    regular = left_cosets(group, {identity}, n)
    weights = (2,) if n % 2 else (1, 1)
    classified = list(subgroups(n))
    assert len(classified) == len(set(classified))
    rows = {}
    for subgroup in classified:
        lifts = tuple(fixed_points(quotient, subgroup, n) for quotient in quotients)
        closure_lifts = fixed_points(regular, subgroup, n)
        value = sum(weight * count for weight, count in zip(weights, lifts)) - closure_lifts
        assert value == 1 - (-1) ** fixed_plane_dimension(subgroup)
        rows[subgroup] = (*lifts, closure_lifts)
    if n == 4:
        # Columns are the two quartic reflection quotients and the octic closure.
        examples = [frozenset((identity,)), *reflections, frozenset((a, 1) for a in range(n))]
        sample = [rows[k] for k in examples]
        assert sample == [(4, 4, 8), (2, 0, 0), (0, 2, 0), (0, 0, 0)]
        print(f"n=4 fixed sheets for K=1, <s>, <rs>, <r>: {sample}.")
    return len(classified)


def check_odd_subtraction(n, base_genus):
    quotient_genus = 1 + n * (base_genus - 1)
    total = 2 ** (2 * quotient_genus) - 1
    ordinary = total - 1  # remove Z -> Y
    for negative, eta, closure in product((0, 1, ordinary // 2, ordinary - 1, ordinary),
                                            (-1, 1), (-1, 1)):
        # Positive-end signs of the ordinary doubles are eta times their twist signs.
        connected_on_y = Fraction(eta * (ordinary - 2 * negative) + closure, 2)
        connected_changes = {n: 2 * eta, 2 * n: 2 * connected_on_y - closure}
        primitive_changes = {}
        for d in range(1, 2 * n + 1):
            change = connected_changes.get(d, 0)
            for r, previous in primitive_changes.items():
                if d % r == 0 and previous:
                    # Through 2n the only proper multiple of a changing label is 2*n.
                    assert (r, d) == (n, 2 * n)
                    change -= previous * Fraction(total, 2)
            primitive_changes[d] = change
        assert [r for r, change in primitive_changes.items() if change] == [n, 2 * n]
        assert primitive_changes[2 * n] == -eta * (1 + 2 * negative)


if __name__ == "__main__":
    limit = 32
    count = sum(check_group(n) for n in range(3, limit + 1))
    for n in range(3, limit + 1, 2):
        for genus in (2, 3, 4):
            check_odd_subtraction(n, genus)
    print(f"Verified the fixed-sheet identity on all {count} subgroups for 3 <= n <= {limit}.")
    print("Odd n: minimum-genus jumps are 2 eta at n and -eta (1 + 2 N_minus) at 2n.")
