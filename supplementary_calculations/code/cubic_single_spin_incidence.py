#!/usr/bin/env python3
"""Check the weighted singular-cubic incidence and first-excess calibration.

Run with sympy==1.14.0. Geometry and orientations are separate proof inputs.
The incidence degree 5 is NOT the stable-map self-intersection S_theta.
"""

from fractions import Fraction

import sympy as sp


def main():
    x, h, z = sp.symbols("xi h z")

    def reduce_base(expr):
        return sp.series(expr, h, 0, 2).removeO().expand()

    incidence = reduce_base((2*x+2*h)*(3*x+3*h)*(3*x+5*h))
    assert incidence == 18*x**3 + 66*h*x**2
    # Weighted projective bundle of lines: push(xi^4)=1/36,
    # xi^5=-h*xi^4, and degree(h)=1 on the genus-two base.
    def ambient_degree(expr):
        poly = sp.Poly(reduce_base(expr), x, h)
        unexpected = [(m, c) for m, c in poly.terms()
                      if m not in ((5, 0), (4, 1)) and c != 0]
        assert not unexpected, unexpected
        return (poly.coeff_monomial(x**4*h)
                - poly.coeff_monomial(x**5))/36

    xx = ambient_degree(incidence*x*x)
    xh = ambient_degree(incidence*x*h)
    assert xx == sp.Rational(4, 3)
    assert xh == sp.Rational(1, 2)

    # W = K + 2 K^2 - theta_p^(-3) K^2, with c1(K)=-xi.
    total = ((1-x*z)*(1-2*x*z)**2)/(1+(-2*x-3*h)*z)
    expansion = sp.series(total, z, 0, 3).removeO().expand()
    c1 = reduce_base(expansion.coeff(z, 1))
    c2 = reduce_base(expansion.coeff(z, 2))
    assert c1 == -3*x+3*h
    assert c2 == 2*x*x-3*x*h
    gamma = reduce_base(c2-x*c1+x*x)
    assert gamma == 6*x*x-6*x*h
    formal_degree = ambient_degree(incidence*gamma)
    assert formal_degree == 5

    # Once the GLOBAL determinantal-divisor Euler formula is proved,
    # the known finite genus-five integrals and pairwise intersection
    # determine each divisor square. This is a different surface from M3.
    equal = Fraction(-70)
    mixed = Fraction(-6)
    pair = Fraction(4)
    difference_square = (2*equal-2*mixed)/4
    pair_of_squares = difference_square+2*pair
    individual_square = pair_of_squares/2
    assert difference_square == -32
    assert pair_of_squares == -24
    assert individual_square == -12
    tail_test = -individual_square/Fraction(12)
    assert tail_test == 1

    print("PASS: incidence class 18*xi^3+66*h*xi^2.")
    print(f"PASS: incidence integrals xi^2={xx}, xi*h={xh}.")
    print(f"PASS: W has c1={c1}, c2={c2}.")
    print(f"PASS: formal incidence class gamma={gamma}; degree={formal_degree}.")
    print("PASS: global first-excess divisor calibration gives D_theta^2=-12.")
    print("PASS: specified ramification-tail test integral -D_theta^2/12=1.")
    print("Not proved by this script: incidence geometry, global Cartier/Euler")
    print("formula, stable reduction, or the genus-six self-intersection S_theta.")
    print("In particular, this check does not identify 5 or 5+1 with S_theta.")


if __name__ == "__main__":
    main()
