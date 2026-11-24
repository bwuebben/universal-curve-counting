#!/usr/bin/env python3
"""Exact algebra for Note 065, using SymPy 1.14.0.

Run: uv run --no-project --with sympy==1.14.0 python checks/cubic_szego_construction.py

These identities check the cubic construction and example. They do not
compute a compact Euler integral or a kernel-locus intersection number.
"""

from sympy import Matrix, Poly, discriminant, expand, factor, gcd, simplify, solve, sqrt, symbols


def main():
    a, b, r, x, y, z, k, ell = symbols("A B r X Y z k ell")
    cubic = b * x**3 - 3 * r * a * x**2 * y + 3 * r**2 * b * x * y**2 - r**3 * a * y**3
    modified = expand(z * cubic.subs({x: k + r * ell / z, y: ell / z}, simultaneous=True))
    expected = z*b*k**3 + 3*r*(b-a)*k**2*ell + 6*r**2*(b-a)*k*ell**2/z + 4*r**3*(b-a)*ell**3/z**2
    assert simplify(modified - expected) == 0
    d = a**2 - b**2
    assert factor(discriminant(cubic.subs(y, 1), x)) == -108 * r**6 * (a-b)**2 * (a+b)**2

    ambient = (-r*a, -b, -r**2*b, -r**3*a)
    changed = (r*(b-a), -z*b, -2*r**2*(b-a)/z, 4*r**3*(b-a)/z**2)
    for alpha, beta, delta, gamma in (ambient, changed):
        linear = 3 * (beta * delta - alpha**2)
        constant = 3*alpha*beta*delta - 2*alpha**3 - beta**2*gamma
        assert simplify(linear + 3*r**2*d) == 0
        assert simplify(constant - 2*r**3*a*d) == 0

    xx, tt = symbols("x t")
    p = lambda v: (v-1)*(v-2)
    rr = lambda v: (v-5)*(v-6)
    first = p(xx)*(xx-3)*rr(tt)*(tt-4) + p(tt)*(tt-3)*rr(xx)*(xx-4)
    second = p(xx)*(xx-4)*rr(tt)*(tt-3) + p(tt)*(tt-4)*rr(xx)*(xx-3)
    numerator = factor((first-second) / (4*(xx-tt)**2))
    assert expand(numerator - (2*xx*tt-7*xx-7*tt+17)) == 0
    matrix = Matrix([[2, -7], [-7, 17]])
    assert matrix.det() == -15
    diagonal = expand(numerator.subs(xx, tt))
    roots = solve(diagonal, tt)
    assert set(roots) == {(7-sqrt(15))/2, (7+sqrt(15))/2}
    branch_polynomial = Poly(1, tt)
    for value in range(1, 7):
        branch_polynomial *= Poly(tt-value, tt)
    assert gcd(Poly(diagonal, tt), branch_polynomial).degree() == 0
    assert Poly(numerator - (2*tt-7)*(xx-tt), xx, tt).as_expr() == diagonal
    print("Modified cubic, binary discriminant, and both actual algebra charts agree exactly.")
    print("Example bidifferential numerator: 2*x*t - 7*x - 7*t + 17; determinant -15.")
    print("Two distinct nonbranch canonical values: (7 ± sqrt(15))/2, giving four points on C.")
    print("No compact Euler integral or intersection multiplicity has been evaluated.")


if __name__ == "__main__":
    main()
