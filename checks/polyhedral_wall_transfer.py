"""Check the fixed-sheet formula for the three polyhedral rotation groups.

The three special orbits on the unit sphere contribute positively and the
regular orbit negatively. This tests the transfer coefficients independently
of the analytic existence of a crossing.
"""

from itertools import permutations

from etale_spin_characters import compose, permutation_sign


def generated_subgroup(identity, generators):
    reached = {identity}
    pending = [identity]
    while pending:
        element = pending.pop()
        for generator in generators:
            other = compose(generator, element)
            if other not in reached:
                reached.add(other)
                pending.append(other)
    return frozenset(reached)


def all_subgroups(group, identity):
    trivial = frozenset({identity})
    seen = {trivial}
    pending = [(trivial, ())]
    while pending:
        subgroup, generators = pending.pop()
        for element in group - subgroup:
            enlarged = generators + (element,)
            candidate = generated_subgroup(identity, enlarged)
            if candidate not in seen:
                seen.add(candidate)
                pending.append((candidate, enlarged))
    return seen


def left_cosets(group, subgroup):
    remaining = set(group)
    cosets = []
    while remaining:
        representative = next(iter(remaining))
        coset = {compose(representative, h) for h in subgroup}
        cosets.append((representative, coset))
        remaining.difference_update(coset)
    return cosets


def fixed_sheets(subgroup, cosets):
    return sum(all(compose(k, representative) in coset for k in subgroup)
               for representative, coset in cosets)


def check_group(name, generators, degrees, expected_subgroups):
    identity = tuple(range(len(generators[0])))
    group = {p for p in permutations(identity)
             if name == "S4" or permutation_sign(p) == 1}
    cyclic = {generated_subgroup(identity, (g,)) for g in group}
    stabilizers = [generated_subgroup(identity, (g,)) for g in generators]
    orbits = [left_cosets(group, h) for h in (*stabilizers, {identity})]
    assert tuple(map(len, orbits)) == degrees

    subgroups = all_subgroups(group, identity)
    assert len(subgroups) == expected_subgroups
    for subgroup in subgroups:
        # Rotations about a common axis form a cyclic group. Two distinct
        # axes leave no common nonzero vector in the standard real space.
        dimension = 3 if len(subgroup) == 1 else int(subgroup in cyclic)
        marks = [fixed_sheets(subgroup, orbit) for orbit in orbits]
        assert sum(marks[:3]) - marks[3] == 1 - (-1) ** dimension

    if name == "S4":
        edge_half_turn = stabilizers[0]
        face_half_turn = generated_subgroup(identity, (compose(generators[2], generators[2]),))
        assert [fixed_sheets(edge_half_turn, orbit) for orbit in orbits] == [2, 0, 0, 0]
        assert [fixed_sheets(face_half_turn, orbit) for orbit in orbits] == [0, 0, 2, 0]
    print(f"{name}: checked all {len(subgroups)} subgroups; orbit degrees {degrees}.")


if __name__ == "__main__":
    # In A4 the vertex orbit and the face-center orbit have conjugate C3
    # stabilizers. In S4, acting on body diagonals, an edge half-turn is a
    # transposition and a face half-turn is the square of a four-cycle.
    check_group("A4", [(1, 0, 3, 2), (1, 2, 0, 3), (1, 2, 0, 3)], (6, 4, 4, 12), 10)
    check_group("S4", [(1, 0, 2, 3), (1, 2, 0, 3), (1, 2, 3, 0)], (12, 8, 6, 24), 30)
    check_group("A5", [(1, 0, 3, 2, 4), (1, 2, 0, 3, 4), (1, 2, 3, 4, 0)], (30, 20, 12, 60), 59)
