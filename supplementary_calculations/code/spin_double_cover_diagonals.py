"""Check diagonal intersections, étale signs and spin exponents against hook coefficients."""

from fractions import Fraction
from math import comb


def spin_parity(label, genus):
    """Arf sign in coordinates (a_i,b_i) of a fixed symplectic basis."""
    mask = (1 << genus) - 1
    return (-1) ** ((label & mask) & (label >> genus)).bit_count()


def check_genus(genus):
    order = 1 << (2 * genus)
    signs = [spin_parity(label, genus) for label in range(order)]
    diagonal_square = -4 * (genus - 1)
    balanced_least = Fraction(order - 1, 2)
    balanced_next = -Fraction((genus - 1) * (5 * order // 2 - 1), 12)
    next_degree_one = -Fraction(genus - 1, 6)
    expected_magnitude = (genus - 1) * order // 8
    # Each translate is a possible difference of theta characteristics. Origins
    # of both parities check the orientation at the embedded curve separately.
    for first in [0, 1 + (1 << genus)]:
        first_row = [signs[first] * signs[first ^ x] for x in range(order)]
        for twist in range(order):
            second = first ^ twist
            epsilon = signs[first] * signs[second]
            def intersection_at_node(x):
                return first_row[x] * signs[second] * signs[second ^ x]

            total = sum(intersection_at_node(x) for x in range(order))
            finite = Fraction(diagonal_square * total, 2 * 4**2)
            etale = Fraction(total - intersection_at_node(0), 2)
            tails = -Fraction(genus - 1, 6) * etale
            least_exponent = etale - epsilon * balanced_least
            next_exponent = finite + tails - epsilon * balanced_next
            next_exponent -= least_exponent * next_degree_one
            if twist:
                assert (finite, finite + tails, next_exponent) == (
                    0, Fraction(genus - 1, 12), epsilon * expected_magnitude
                )
            else:
                assert (finite, finite + tails, next_exponent) == (
                    -expected_magnitude, balanced_next, 0
                )
    print(f"genus {genus}, {order} twists, absolute mixed exponent {expected_magnitude}")


def check_smith_factors(genus):
    """Diagonalize the character matrix using integral elementary operations."""
    dimension = 2 * genus
    size = 1 << dimension
    matrix = [[(-1) ** (row & col).bit_count() for col in range(size)]
              for row in range(size)]
    for axis in range(dimension):
        bit = 1 << axis
        for upper in range(size):
            if upper & bit:
                continue
            lower = upper | bit
            for col in range(size):
                matrix[lower][col] -= matrix[upper][col]
            for row in range(size):
                matrix[row][lower] -= matrix[row][upper]
    assert all(matrix[i][j] == 0 for i in range(size) for j in range(i))
    assert all(matrix[i][j] == 0 for i in range(size) for j in range(i + 1, size))
    actual = sorted(abs(matrix[i][i]) for i in range(size))
    expected = [2**j for j in range(dimension + 1) for _ in range(comb(dimension, j))]
    assert actual == expected
    print(f"integral diagonalization of the {size} by {size} character matrix verified")


if __name__ == "__main__":
    for genus in range(2, 6):
        check_genus(genus)
    for genus in (2, 3):
        check_smith_factors(genus)
