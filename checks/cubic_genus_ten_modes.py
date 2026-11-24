# /// script
# dependencies = ["sympy==1.14.0"]
# ///
"""Compute genus-ten odd modes and the leading weight at each source genus."""

from functools import cache

import sympy as sp

from local_curve_degree_three import build, in_w, w


@cache
def connected_residue(weight):
    traces = {}
    for degree in (2, 3):
        handle = build(degree, weight, 1)[-1]
        traces[degree] = sp.Poly(in_w(sp.cancel(sp.trace(handle**2))), w)
    return sp.Poly(traces[3].as_expr() - w**2 * traces[2].as_expr() + w**6 / 3, w)


def odd_coefficients(power):
    # At source genus at most ten only symmetric odd modes 1 and 3 can occur.
    def row(weight):
        z = sp.Integer(weight)
        return [2 * (z**j + z**-j) for j in (1, 3)]

    def difference(weight):
        return connected_residue(weight).nth(power) - connected_residue(-weight).nth(power)

    matrix = sp.Matrix([row(1), row(2)])
    modes = matrix.inv() * sp.Matrix([difference(1), difference(2)])
    assert (sp.Matrix([row(3)]) * modes)[0] == difference(3)
    return tuple(modes)


def leading_weight_coefficients():
    # At w^(6+j) the symmetric weight polynomial has degree at most j in z+1/z.
    # Its leading coefficient is the coefficient of either z^j or z^(-j).
    weights = (1, -1, 2, -2, 3, -3, 4, -4)
    variable = sp.Symbol("Y")
    rows = []
    for order in range(1, 7):
        samples = [
            (weight + sp.Rational(1, weight), connected_residue(weight).nth(6 + order))
            for weight in weights
        ]
        polynomial = sp.Poly(sp.interpolate(samples[:order + 1], variable), variable)
        assert all(polynomial.eval(y) == value for y, value in samples[order + 1:])
        rows.append(polynomial.nth(order))
    return rows


if __name__ == "__main__":
    rows = {power: odd_coefficients(power) for power in (7, 8, 9)}
    assert rows[7][1] == rows[8][1] == 0
    assert rows[9] == (sp.Rational(674, 3), sp.Rational(506, 3))
    assert 4 * sum(rows[9]) == sp.Rational(4720, 3)
    assert connected_residue(2).nth(9) == 549
    assert connected_residue(-2).nth(9) == -sp.Rational(19891, 6)
    print("Odd modes (1, 3) in powers of w:", rows)

    # w = 2-2cos(u). Recover individual connected stable-map integrals.
    source = {}
    for mode in range(2):
        genus8 = rows[7][mode]
        genus9 = rows[8][mode] - sp.Rational(7, 12) * genus8
        genus10 = rows[9][mode] - sp.Rational(2, 3) * genus9 - sp.Rational(161, 720) * genus8
        source[(1, 3)[mode]] = (genus8, genus9, genus10)
    print("Connected coefficients at source genera 8, 9, 10:", source)
    diagonal = leading_weight_coefficients()
    assert diagonal == [
        3612, 1228, sp.Rational(506, 3), sp.Rational(21, 2),
        sp.Rational(8, 27), sp.Rational(1, 324),
    ]
    print("Leading weights in source genera 8 through 13:", diagonal)
