"""Symbolic local-chart, coefficient and crossing checks for Notes 102–110.

Run with uv run --no-project --with sympy==1.14.0 python checks/two_copy_quartic_smoothing.py.
"""

from itertools import product

import sympy as sp


def check_two_copy_chart():
    a, b, p, z = sp.symbols("a b p z", nonzero=True)
    s, q, x = sp.symbols("s q x")
    y = q / z
    target = x + z**2 + 2*p*z + y**2 + 2*s*y
    for sign in (-1, 1):
        square = (z + sign*y + p)**2
        assert sp.expand(target.subs({s: sign*p, x: p**2 + 2*sign*q}) - square) == 0

    rows = [a*p - b*s, x - p**2 - 2*b*q/a, x - s**2 - 2*a*q/b]
    balanced, = sp.solve(rows, (s, x, q), dict=True)
    assert sp.simplify(balanced[q] + a*p**2/(2*b)) == 0
    assert balanced[x] == 0
    mass = sp.factor((2*q*(b*p - a*s)).subs(balanced))
    assert sp.simplify(mass - a*(a**2-b**2)*p**3/b**2) == 0

    # In the framed signed permutation group count elements fixing the kernel line.
    for height, expected in [((2, 1), 2), ((1, 0), 4), ((1, 1), 4)]:
        fixed = 0
        for first, second, swap in product((-1, 1), (-1, 1), (False, True)):
            lhs, rhs = (height[::-1] if swap else height)
            fixed += first*lhs*height[1] == second*rhs*height[0]
        assert fixed == expected


def check_double_zero_residual():
    # Use independent formal conjugates; u is the angle on the unit circle.
    e, u, p, pb, lam, lb = sp.symbols("e u p pbar lambda lambdabar", nonzero=True)

    def sheet(w, shift, smoothing):
        branch = w * sp.sqrt((1-e*shift/w)**2 - e**2*smoothing/w**2)
        return sp.series(branch, e, 0, 4).removeO().expand()

    d = sheet(u, p, lam)
    db = sheet(1/u, pb, lb)
    height = sp.expand(d*(u + e*p))
    assert height.coeff(e, 2) == -p**2-lam/2
    assert height.coeff(e, 3) == -lam*p/u

    # The true residual in the odd coefficient is divided by d. Extract its
    # radial term by taking the constant Laurent mode in the angle.
    numerator = sp.expand(db*(1/u + e*pb) - height/u**4)
    residual = sp.series(numerator / d, e, 0, 4).removeO().expand()
    radial = [residual.coeff(e, n).coeff(u, 0) for n in range(4)]
    assert radial == [0, 0, 0, -lb*pb]
    assert radial[3].subs(lb, -2*pb**2) == 2*pb**3
    # Conjugation on a unit circle sends pbar to the reciprocal angle.
    assert sp.expand(radial[3].subs({lb: -2/u**2, pb: 1/u})) == 2/u**3


def check_primitive_subtraction():
    t, s, eta, r, c, f2, tail = sp.symbols("T S eta R c F2 tail")
    minimum = eta*(t - 2*s)
    known_next = eta*(-tail*t + s*(2*tail + 2) + 1)
    proper_divisor = eta*(f2 - tail*c)
    primitive_tail = -tail*(minimum - eta*c)
    residual = sp.expand(r + known_next - proper_divisor - primitive_tail)
    assert residual == r + 2*s*eta + eta - eta*f2
    h = 4  # lower genus when g(C) = 2
    assert residual.subs(f2, -(h-1)*2**(2*h-3)) == r + eta*(2*s) + 97*eta


def check_displaced_chambers():
    lower_sign, ordinary, pole1, pole2, tower = sp.symbols("eta F E1 E2 U")

    # Increasing displacement reverses the upper time over C. The ordinary
    # and both pole integrals decrease by the same signed lower crossing.
    marked_before = -lower_sign*(4*ordinary + pole1 + pole2)
    marked_after = -lower_sign*(4*(ordinary-lower_sign)
                               + pole1-lower_sign + pole2-lower_sign)
    tower_after = tower - 2*3*lower_sign
    difference = sp.expand(marked_after - marked_before + lower_sign*(tower_after - tower))
    assert difference == 0

    genus = 2
    for offset in range(2, 18, 2):
        x = sp.Symbol("x")
        polynomial = sp.Poly(2*(4-x)**(genus-1)*(sp.chebyshevt(offset, 1-x/2)-1), x)
        assert abs(polynomial.LC()) == 1
        source_genus = 2*genus - 1 + polynomial.degree()
        branch_count = 2*(source_genus - 1) - 2*(2*genus - 2)
        assert source_genus == 4+offset and branch_count == 2+2*offset
        odd_direct_image_degree = sp.Rational(2*source_genus-4-branch_count, 2)
        assert odd_direct_image_degree + 2*(1-genus) == -1
        assert 3*branch_count == 6+6*offset
    assert 2**8*(4-1)*(4-2)//16 == 96
    assert 2**4*genus*(genus-1)//4 == 8


if __name__ == "__main__":
    check_two_copy_chart()
    check_double_zero_residual()
    check_primitive_subtraction()
    check_displaced_chambers()
    print("Passed: local charts, stabilizers, residues, primitive subtraction, chamber cancellation and branch counts.")
