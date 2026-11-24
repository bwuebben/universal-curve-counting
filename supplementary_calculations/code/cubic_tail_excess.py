#!/usr/bin/env python3
"""Exact local algebra and arithmetic for Note 066.

Run: uv run --no-project --with sympy==1.14.0 python checks/cubic_tail_excess.py

This verifies the displayed ramification-chart ideal and rational calibration.
It does not prove the compact support classification, the Fredholm comparison,
the stable-map versal chart, or the geometric Hodge integral.
"""

from fractions import Fraction as F

import sympy as sp


def main():
    x, v, t, alpha, beta = sp.symbols("x v t alpha beta")
    field = sp.QQ.frac_field(alpha, beta)
    ga, gb = v + alpha * t, v + beta * t
    assert sp.det(sp.Matrix([[sp.diff(ga, v), sp.diff(ga, t)],
                            [sp.diff(gb, v), sp.diff(gb, t)]])) == beta - alpha
    # Both divisors meet the triple-ramification divisor t=0 transversely.
    assert ga.subs(t, 0) == gb.subs(t, 0) == v
    # On the ramification-marked surface, t=-3*x**2. The quotient has
    # basis 1,x over the coefficient field: the local intersection length is 2.
    groebner = sp.groebner([ga.subs(t, -3*x*x), gb.subs(t, -3*x*x)],
                          v, x, order="lex", domain=field)
    assert [p.as_expr() for p in groebner.polys] == [v, x*x]
    assert groebner.reduce(x)[1] == x
    assert groebner.reduce(x*x)[1] == 0

    # Check that the already established explicit Szego example has four
    # simple points on the hyperelliptic curve, not only two base roots.
    q = 2*t*t - 14*t + 17
    branch = sp.prod(t-j for j in range(1, 7))
    assert sp.discriminant(q, t) == 60
    assert sp.gcd(q, branch) == 1
    covers = 2 * sp.degree(q, t)
    assert covers == 4

    length = 2
    normal_degree = -F(1, 24)  # Geometric input: c1(N)|E=-lambda, integral 1/24.
    per_cover = length * normal_degree
    mixed = int(covers) * per_cover
    assert per_cover == -F(1, 12)
    assert mixed == -F(1, 3)

    equal_main = F(161, 6)
    exponent_offset = F(43, 2)
    # I_ab = I_equal - 2*(S_a+S_b) + 4*T_ab.
    assert equal_main + 4*mixed == F(51, 2)
    assert equal_main + 4*mixed - exponent_offset == 4
    # The candidate exponent -20 would require S_a+S_b=12; no value
    # of either self-intersection is derived by this check.
    assert (4 - F(-20))/2 == 12

    print("PASS: transverse core ideals pull back to (v,x^2), of length 2.")
    print("PASS: explicit Szego example gives four simple cover points.")
    print(f"PASS: per-cover pairing {per_cover}; total mixed pairing {mixed}.")
    print("PASS: I_ab=51/2-2*(S_a+S_b); n_3,6=4-2*(S_a+S_b).")
    print("Not checked numerically: global support, orientations, family comparison,")
    print("Hodge degree, or either unknown self-intersection or cubic exponent.")


if __name__ == "__main__":
    main()
