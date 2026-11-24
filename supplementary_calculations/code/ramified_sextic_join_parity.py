"""Discrete checks of Note 119. The analytic local degree is not an input to these tests."""

from fractions import Fraction
from itertools import permutations

from etale_spin_characters import compose, hook_count
from ramified_wall_leading_degree import direct_pointed_count, odd_part_is_square, pointed_counts


def permutation_sign(permutation):
    inversions = sum(a > b for i, a in enumerate(permutation) for b in permutation[i + 1:])
    return (-1) ** inversions


def check_residues(group, tau):
    identity = (0, 1, 2)
    inverse = {g: tuple(g.index(i) for i in identity) for g in group}
    # A simple pole changes sign when the local ramification coordinate is negated.
    # Sections of the induced sign representation are odd functions on the right cosets.
    selected_flag = {identity: 1, tau: -1}
    for g in group:
        alternating = sum(permutation_sign(k) * selected_flag.get(compose(inverse[k], g), 0)
                          for k in group)
        assert Fraction(alternating, len(group)) == Fraction(permutation_sign(g), 3)

    i, j = [sheet for sheet in identity if tau[sheet] != sheet]
    remaining = next(sheet for sheet in identity if tau[sheet] == sheet)
    root = [0, 0, 0]
    root[i], root[j] = 1, -1
    fixed_vector = [1, 1, 1]
    fixed_vector[remaining] = -2
    assert fixed_vector[i] == fixed_vector[j]
    for sheet, velocity in ((i, 1), (j, -1)):
        moving_pole = velocity * (root[sheet] - root[remaining])
        common_pole = fixed_vector[sheet] - fixed_vector[remaining]
        assert 3 * moving_pole == common_pole


def main():
    group = list(permutations(range(3)))
    for tau in group:
        if permutation_sign(tau) == -1:
            check_residues(group, tau)
    print("All three inertia choices give sign projection 1/3 and moving cubic branch ratio 1/3.")

    for genus in (2, 3, 4):
        degrees = range(1, 19)
        disconnected = [1] + [hook_count(m, genus, False) for m in degrees]
        pointed = pointed_counts(disconnected)
        for m in degrees:
            assert pointed[m] % 2 == odd_part_is_square(m)
        if genus == 2:
            for m in (1, 2, 3, 4):
                assert pointed[m] == direct_pointed_count(m, genus)
            for m in (1, 2, 4):
                d = 6 + m
                genus_at_join = 1 + d * (genus - 1) + 3 + 1
                print(f"Genus-two base, two branch values: d={d}, k={genus_at_join}, "
                      f"{pointed[m]} join chart(s) per branch value.")
    print("Pointed-cover parity agrees with the square/twice-square criterion through m=18 "
          "for base genera 2, 3 and 4.")


if __name__ == "__main__":
    main()
