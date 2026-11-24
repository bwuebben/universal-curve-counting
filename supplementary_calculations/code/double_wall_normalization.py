"""Exact arithmetic for the wall formulas of Notes 079 and 080."""

from fractions import Fraction
from math import comb


def normalization(base_genus, source_genus):
    c = base_genus - 1
    j = source_genus - (2 * base_genus - 1)
    if not 1 <= j <= c:
        return 0
    return 2 * 4 ** (c - j) * (comb(c, 2 * j) - (-1) ** j * comb(c, j))


def spin_wall_count(base_genus, branch_pairs):
    c = base_genus - 1
    if branch_pairs > c:
        return 0
    return -2 * (-1) ** branch_pairs * comb(c, branch_pairs) * 4 ** (c - branch_pairs)


def main():
    for g in range(2, 33):
        c = g - 1
        minimum = 2 * g - 1
        # After removing the common x^(2c), the two sides have degree at
        # most c in x. Evaluate at c+1 distinct rational squares.
        for i in range(1, c + 2):
            s = Fraction(i, c + 3)
            x = s**2
            scalar = (4 - 2 * s) ** c + (4 + 2 * s) ** c
            balanced = 2 * (4 - x) ** c
            correction = sum(normalization(g, minimum + j) * x**j for j in range(c + 1))
            spin_difference = sum(spin_wall_count(g, j) * x**j for j in range(c + 1))
            assert correction == scalar - balanced
            assert spin_difference == -balanced
    for s in [Fraction(1, 7), Fraction(3, 5), Fraction(5, 4)]:
        assert 1 / (4 - 2 * s) + 1 / (4 + 2 * s) == 2 / (4 - s**2)
    for m in range(2, 9):
        for g in range(2, 9):
            for branch_pairs in range(1, 9):
                h = 1 + m * (g - 1) + branch_pairs
                moving_tail = -Fraction(m * (g - 1), 12)
                ramification_terms = -Fraction(2 * branch_pairs, 24)
                assert moving_tail + ramification_terms == -Fraction(2 * h - 2, 24)
    assert [spin_wall_count(4, j) for j in range(3)] == [-128, 96, -24]
    assert [normalization(3, h) for h in (5, 6, 7)] == [0, 24, -2]
    examples = []
    for degree in range(3, 28):
        genus = degree * (degree - 1) // 2
        sections = (degree + 1) * (degree + 4) // 2
        primitive_count = (-1) ** (sections - 1) * sections
        h = 3 * genus - 2
        assert all(normalization(g, h) == 0 for g in range(genus))
        discrepancy = primitive_count * normalization(genus, h)
        assert discrepancy == (-1) ** (degree + 1) * (degree + 1) * (degree + 4)
        if degree in (3, 4, 15, 27):
            examples.append((genus, h, discrepancy))
    assert examples == [(3, 7, 28), (6, 16, -40), (105, 313, 304), (351, 1051, 868)]
    print("Scalar, balanced and spin expansions verified for base genera 2–32.")
    print("One-tail decomposition verified against the primitive sine coefficient.")
    print("Compact examples (primitive genus, output genus, difference):", examples)


if __name__ == "__main__":
    main()
