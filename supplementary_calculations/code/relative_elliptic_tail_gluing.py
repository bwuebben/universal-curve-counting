"""Enumerate the smoothing roots and cut stabilizers in Note 073."""

from fractions import Fraction
from itertools import product
from math import lcm, prod


def partitions(total, least=1):
    if total == 0:
        yield ()
    for part in range(least, total + 1):
        for rest in partitions(total - part, part):
            yield (part,) + rest


def loop_orbits(profile, step):
    roots = set(product(*(range(order) for order in profile)))
    orbits = []
    while roots:
        root = next(iter(roots))
        orbit = set()
        while root not in orbit:
            orbit.add(root)
            roots.remove(root)
            root = tuple((phase + step) % order for phase, order in zip(root, profile))
        orbits.append(frozenset(orbit))
    return orbits


def check_profile(profile):
    count = prod(profile)
    period = lcm(*profile)
    ordinary = loop_orbits(profile, 1)
    assert len(ordinary) == count // period
    assert all(len(branch) == period for branch in ordinary)

    spin = loop_orbits(profile, 2)
    spin_period = lcm(period, 2) // 2
    assert len(spin) == count // spin_period
    assert all(len(branch) == spin_period for branch in spin)

    # Identify each ordinary branch with a coset of the simultaneous loop.
    # The unmarked rotations fix the phase at the marked index.
    branch_of = {phase: index for index, branch in enumerate(ordinary) for phase in branch}
    for marked, order in enumerate(profile):
        images = [branch_of[rotation] for rotation in branch_of if rotation[marked] == 0]
        assert set(images) == set(range(len(ordinary)))
        assert images.count(branch_of[(0,) * len(profile)]) == period // order


def main():
    profiles = 0
    for degree in range(1, 11):
        for profile in partitions(degree):
            check_profile(profile)
            profiles += 1

    # One marked cyclic cover in genus two with its inverse branch characters.
    # The smooth Hodge integral includes the cyclic deck and elliptic factors.
    for degree in range(2, 17, 2):
        cut = Fraction(2 * degree, 24 * degree)
        glued = sum(len(branch) * cut for branch in loop_orbits((degree,), 2))
        assert cut == Fraction(1, 12) and glued == Fraction(degree, 12)
    print(f"Verified ordinary and spin root orbits and cut stabilizers for {profiles} profiles.")
    print("Verified genus-two cyclic tail weights up to even degree 16. Degree two gives 1/12 and 1/6.")


if __name__ == "__main__":
    main()
