"""Compare the cubic wall sums with the independent low-genus endpoint calculations.

Run from the project root with
uv run --no-project --with sympy==1.14.0 python checks/cubic_wall_counts.py.
"""

from math import comb

from sympy import Poly, QQ

from etale_elementary_exponents import elementary_exponents
from first_excess_spin_hurwitz import first_exponents
from nonorientable_spin_partition import extract
from oriented_holonomy_genera import connected, z


def wall_sums(base_genus):
    c = base_genus - 1
    etale = 4**c * (1 - 9**c)
    two_branch = c * 4 ** (c - 1) * (9**c - 1 + 8 * 3 ** (2 * c - 1))
    # The first two coefficients receive contributions only from cubic walls.
    # At the next genus each two-branch double crossing has two joining nodes.
    second_coefficient = (
        -4**c * comb(2 * c, 2) * 3 ** (2 * c - 2)
        - c * 4 ** (c - 1) * (2 * c) * 3 ** (2 * c - 1)
    )
    if c >= 2:
        second_coefficient += comb(c, 2) * 4 ** (c - 2) * (1 - 9**c)
    double_two_branch = 2 * c * 4 ** (c - 1)
    four_branch = second_coefficient + double_two_branch
    return etale, two_branch, four_branch


def check_all_genera(c, double_difference, cubic_difference):
    double_counts = double_difference.exquo(Poly(z ** (2 * c), z, domain=QQ))
    cubic_change = cubic_difference.exquo(Poly(z ** (3 * c), z, domain=QQ))
    mixed = Poly(0, z, domain=QQ)
    for (excess,), count in double_counts.terms():
        for joins in range(1, 2 * excess + 1):
            term = count * (-1) ** joins * comb(2 * excess, joins)
            mixed += Poly(term * z ** (excess + joins), z, domain=QQ)

    smooth = (cubic_change - mixed).exquo_ground(2)
    predicted = Poly(
        (4 - z * (1 - z) ** 2) ** c - ((4 - z) * (3 - z) ** 2) ** c,
        z,
        domain=QQ,
    )
    assert smooth == predicted
    assert smooth.degree() == 3 * c - 1
    assert smooth.LC() == -8 * c * (-1) ** (c - 1)
    assert all(value.q == 1 for value in smooth.all_coeffs())
    return smooth


def check(base_genus):
    c = base_genus - 1
    mixed_log, _ = extract(c, 1, 3)
    difference = mixed_log[3] - connected(3, c, 0)
    etale, branched, four_branch = wall_sums(base_genus)
    assert difference.nth(3 * c) == 2 * etale
    assert difference.nth(3 * c + 1) == 2 * branched
    double_difference = mixed_log[2] - connected(2, c, 0)
    double_two_branch = double_difference.nth(2 * c + 1)
    assert double_two_branch == 2 * c * 4 ** (c - 1)
    assert difference.nth(3 * c + 2) == 2 * four_branch - 2 * double_two_branch
    smooth = check_all_genera(c, double_difference, difference)
    assert tuple(smooth.nth(j) for j in range(3)) == (etale, branched, four_branch)

    # Each degree-three proper-divisor term depends only on degree one.
    # The endpoint degree-one series agree, so their contributions cancel.
    first_change = first_exponents(base_genus, 3, True, 1)[3] - first_exponents(
        base_genus, 3, False, 1
    )[3]
    least_change = elementary_exponents(base_genus, 3, True, 1)[3] - elementary_exponents(
        base_genus, 3, False, 1
    )[3]
    assert (least_change, first_change) == (2 * etale, 2 * branched)
    if base_genus == 2:
        assert (etale, branched, four_branch) == (-32, 32, -8)
    return etale, branched, four_branch


if __name__ == "__main__":
    for g in range(2, 10):
        etale_sum, two_branch_sum, four_branch_sum = check(g)
        print(
            f"g={g}: signed cubic sums {etale_sum} (unramified), "
            f"{two_branch_sum} (two branch points), {four_branch_sum} (four branch points); "
            f"full polynomial verified through highest source genus {6 * (g - 1)}"
        )
