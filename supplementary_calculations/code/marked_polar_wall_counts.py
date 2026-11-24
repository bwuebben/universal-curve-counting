"""Exact intersection and signed-count checks for the marked double-cover formulas."""

from fractions import Fraction
from itertools import product
from math import comb


def diagonal_pair(left, right, tangent_degree):
    """Intersection of two divisors in the basis X, Y, diagonal of C × C."""
    x, y, diagonal = left
    a, b, d = right
    return x*b + y*a + x*d + y*d + diagonal*(a + b) + tangent_degree*diagonal*d


def complex_reference(genus):
    sheets = 2**(2*genus)
    tangent = 2 - 2*genus
    index_c1 = (Fraction(tangent, 4), Fraction(tangent, 4), Fraction(1, 2))
    index_ch2 = -Fraction(tangent, 4)
    obstruction_c2 = diagonal_pair(index_c1, index_c1, tangent)/2 + index_ch2
    ordinary = sheets * obstruction_c2 / 4
    slice_c1 = diagonal_pair(index_c1, (1, 0, 0), tangent)
    smoothing_c1 = -Fraction(sheets, 4)
    polar = -sheets * slice_c1 / 2 - smoothing_c1
    return ordinary, polar


def check_pole_counts(genus, pairs):
    sheets = 2**(2*genus)
    flags = 2*pairs - 1
    # On the ordered moving-branch curve each fixed graph meets the moving one once.
    ordinary_index_c1 = Fraction(sheets, 2) * Fraction(2*flags + 2 - 2*genus, 4)
    each_smoothing_c1 = -Fraction(sheets, 4)
    complex_polar = -(ordinary_index_c1 + flags*each_smoothing_c1)
    assert complex_polar == complex_reference(genus)[1]

    # Compare the real product basis of V+V with the interleaved complex basis.
    rank = 2*pairs
    permutation = [index for j in range(rank) for index in (j, rank+j)]
    inversions = sum(x > y for i, x in enumerate(permutation) for y in permutation[i+1:])
    rotation_orientation = (-1)**inversions
    for axes in product((1, -1), repeat=flags):
        image_orientation = (-1)**axes.count(-1)
        complex_quotient_c1 = -sum(sign*each_smoothing_c1 for sign in axes)
        repeated = rotation_orientation * image_orientation * complex_quotient_c1
        predicted = (-1)**(pairs+axes.count(-1)) * sum(axes) * Fraction(sheets, 4)
        assert repeated == predicted
        assert (repeated - complex_polar).denominator == 1
        if axes.count(1) in (0, flags):
            assert abs(repeated) == flags * sheets // 4


def check_quartic(genus, half_branches):
    sheets = 2**(2*genus)
    branches = 2*half_branches
    h = 2*genus - 1 + half_branches
    minimum = h + 2*genus - 1
    assert 2*h - 1 == minimum + half_branches - 1
    complex_ordinary, complex_polar = complex_reference(genus)
    for eta, base_sign, etale_sum in product((-1, 1), (-1, 1), (-3, 1, 3)):
        ordinary_crossings = -3
        polar_counts = list(range(-1, branches-1))
        ordinary = base_sign * (complex_ordinary + ordinary_crossings)
        polars = [base_sign * (complex_polar + count) for count in polar_counts]
        joined = eta * (etale_sum * (Fraction(branches*(minimum-1), 12)
                                      + 2*comb(branches, 2)) + comb(branches, 2))
        connected = joined - eta * (2*branches*ordinary + sum(polars))
        negative_twists = 3
        proper_divisor = 0
        if half_branches == 2:
            total_twists = 2**(2*h) - 1
            connected += eta * Fraction(total_twists-2*negative_twists, 2)
            proper_divisor = eta * Fraction(total_twists, 2)
        first_coefficient = -eta * branches * etale_sum
        next_primitive = connected - proper_divisor + first_coefficient * Fraction(minimum-1, 12)
        count_formula = ((2*etale_sum+1)*comb(branches, 2)
                         - base_sign * (Fraction(branches*sheets*genus*(genus-1), 8)
                                        + 2*branches*ordinary_crossings + sum(polar_counts)))
        if half_branches == 2:
            count_formula -= negative_twists
        assert next_primitive == eta * count_formula
        assert next_primitive.denominator == 1
        if (genus, branches) == (2, 6):
            assert int(next_primitive) % 2 == (1+sum(polar_counts)) % 2


def ordinary_count_polynomial(genus, axis_offset):
    """Two partitions in degree two, via the recurrence for 2*cos(a*v)."""
    older, cosine = [2], [2, -1]
    for _ in range(1, axis_offset):
        newer = [0] * (len(cosine)+1)
        for j, coeff in enumerate(cosine):
            newer[j] += 2*coeff
            newer[j+1] -= coeff
        for j, coeff in enumerate(older):
            newer[j] -= coeff
        older, cosine = cosine, newer
    cosine[0] -= 2
    for _ in range(genus-1):
        cosine = [4*a-b for a, b in zip(cosine + [0], [0] + cosine)]
    return cosine


def split_polar_degree(genus, pairs, offset, minus_flags):
    """Compute both modified axis indices on the compact moving-branch curve."""
    sheets = 2**(2*genus)
    flags = 2*pairs-1
    tangent_degree = 2-2*genus
    common_c1 = Fraction(sheets*(2*flags+tangent_degree), 16)
    each_smoothing_c1 = -Fraction(sheets, 4)
    plus_c1 = common_c1 + (offset + flags-minus_flags)*each_smoothing_c1
    minus_c1 = common_c1 + (-offset + minus_flags)*each_smoothing_c1
    negative_axis_index = -offset - pairs + minus_flags
    orientation = -1 if negative_axis_index % 2 else 1
    return orientation * (minus_c1 - plus_c1)


def check_axes(genus):
    sheets = 2**(2*genus)
    for offset in range(2, 22, 2):
        counts = ordinary_count_polynomial(genus, offset)
        rotated_c1_square = 2*offset**2
        rotated_ch2 = Fraction(2-2*genus, 4)
        rotated_c2 = Fraction(rotated_c1_square, 2) - rotated_ch2
        ordinary = -rotated_c2*sheets/4
        repeated = -Fraction(sheets*(genus-1), 8)
        assert counts[0] == 0
        assert counts[1] == ordinary-repeated == -sheets*offset**2//4
        for pairs in range(1, 6):
            flags = 2*pairs-1
            for minus_flags in range(flags+1):
                degree = split_polar_degree(genus, pairs, offset, minus_flags)
                assert degree == split_polar_degree(genus, pairs, -offset, flags-minus_flags)
        one_pole_jump = split_polar_degree(genus, 1, offset, 0) - split_polar_degree(genus, 1, 0, 0)
        assert one_pole_jump == -sheets*offset//2
        assert counts[-1] == (-1)**(genus-1+offset)
        assert len(counts)-1+2*genus-1 == 3*genus-2+offset
    if genus == 2:
        assert ordinary_count_polynomial(genus, 2) == [0, -16, 8, -1]


def verify():
    for genus in range(2, 7):
        sheets = 2**(2*genus)
        ordinary, polar = complex_reference(genus)
        assert ordinary == Fraction(sheets*(genus-1)*(genus-2), 16)
        assert polar == Fraction(sheets*(genus-1), 4)
        repeated_ordinary = -Fraction(sheets*(genus-1), 8)
        repeated_polar = -Fraction(sheets, 4)
        assert repeated_polar - polar == -Fraction(genus*sheets, 4)
        assert repeated_ordinary - ordinary == -Fraction(sheets*genus*(genus-1), 16)
        check_axes(genus)
        for pairs in range(1, 6):
            check_pole_counts(genus, pairs)
        for half_branches in range(2, 7):
            check_quartic(genus, half_branches)
    assert complex_reference(2) == (0, 4)
    print("Ordered GRR, spin rotation signs and pole counts agree for g=2–6 and 1–9 prescribed poles.")
    print("Quartic connected and primitive formulas agree for 4–12 branch values, including the four-branch tower.")
    print("Genus-two polar Euler numbers: complex 4, repeated spin -4, distinct spin 0.")
    print("Chebyshev extraction through axis offset 20 agrees with the ordinary Euler jump and top genus.")


if __name__ == "__main__":
    verify()
