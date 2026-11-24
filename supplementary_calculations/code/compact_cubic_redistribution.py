"""Exact arithmetic for the compact cubic crossing and full-subset resonance."""

from itertools import combinations
from math import comb

from sympy import Matrix, ZZ, discriminant, prod, symbols
from sympy.matrices.normalforms import smith_normal_form

from degree_seven_resonance_lattice import jump_vector, subset_sums


def two_zero_blocks(size, first_block):
    """Only the two complete blocks and their union have zero nonempty sum."""
    second_block = sorted(set(range(size)) - set(first_block))
    weights = [0] * size
    for block in (first_block, second_block):
        for i in block[:-1]:
            weights[i] = 8 << i
        weights[block[-1]] = -sum(weights[i] for i in block[:-1])
    return weights, second_block


def full_wall_pair(weights):
    # Stay close to this relative chamber while crossing only the total sum.
    vectors = []
    for displacement in (-1, 1):
        side = [4 * w for w in weights]
        side[-1] += displacement
        vectors.append(jump_vector(subset_sums(side), len(side)))
    return vectors


def independent_coordinates(vector):
    half = len(vector) // 2
    assert all(vector[j - 1] == -vector[-j - 1] for j in range(1, half))
    assert vector[half - 1] == 0
    return vector[:half - 1] + vector[-1:]


def check_full_subset(size):
    half = size // 2
    # On the full hyperplane this starting point has no vanishing proper sum.
    start = [1 << i for i in range(size - 1)]
    start.append(-sum(start))
    lower, upper = [independent_coordinates(v) for v in full_wall_pair(start)]
    basis = [lower, [b - a for a, b in zip(lower, upper)]]
    assert basis[-1] == [0] * (half - 1) + [2]
    for j in range(1, half):
        # Up to relabeling, this covers every complementary pair of these sizes.
        first = list(range(j))
        weights, second = two_zero_blocks(size, first)
        reduced = []
        for direction in (-1, 1):
            relative = weights.copy()
            relative[first[-1]] += direction
            relative[second[-1]] -= direction
            reduced.append(independent_coordinates(full_wall_pair(relative)[1]))
        difference = [b - a for a, b in zip(*reduced)]
        assert difference == [2 * (-1) ** j if k == j else 0 for k in range(1, half + 1)]
        basis.append(difference)
    smith = smith_normal_form(Matrix(basis), domain=ZZ)
    assert abs(prod(smith[i, i] for i in range(half))) == 2 ** (half - 1)
    return 2 ** (half - 1)


def main():
    z, w, q = symbols("z w q")
    assert discriminant(z ** 3 - w * z + q, z) == 4 * w ** 3 - 27 * q ** 2
    print("The three new branch values solve 4w³ = 27q² in the centered model.")
    for size in range(2, 15, 2):
        interior = jump_vector(subset_sums([1] * size), size)
        exterior = [-value for value in interior]
        assert exterior == [(-1) ** (j + 1) * comb(size, j) for j in range(1, size + 1)]
        assert sum(exterior) == 1 and sum(map(abs, exterior)) == 2 ** size - 1
        index = check_full_subset(size)
        print(f"b={size}: maximum absolute changes sum to {2 ** size - 1}, "
              f"full-subset lattice index {index}.")
    # On either side of the full resonance with two flags, the singleton signs cancel.
    assert {tuple(v) for v in full_wall_pair([1, -1])} == {(0, -1), (0, 1)}
    # A nonzero rational bivector is cyclic under coordinate scalings/permutations.
    # Check the Vandermonde step used in dimension six in the geometric proof.
    pairs = list(combinations(range(6), 2))
    weights = [2 ** i + 2 ** j for i, j in pairs]
    assert len(set(weights)) == len(pairs) == 15
    print("The 15 diagonal weights on the exterior square in dimension six are distinct.")


if __name__ == "__main__":
    main()
