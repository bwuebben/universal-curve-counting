"""Check the cusp length, oriented intersections, and the genus-five handle polynomial."""

from math import comb, gcd

from sympy import Poly, Rational, cos, groebner, symbols


def scalar_minus_balanced(w):
    """Paper 6, arithmetic section, fourth handle power (core genus five)."""
    n = 4
    factor = w**2 + 54*w - 243
    a = w**2 + 50*w + 81
    b = 4*(w - 6)*factor
    c = 4*(w - 4)*(w - 3)**2*factor
    traces = [3, a, a**2 - 2*b/3]
    for i in range(3, n + 1):
        traces.append(a*traces[i-1] - b*traces[i-2]/3 + c*traces[i-3]/3)
    doubles = 2*sum(
        4**(n-i)*(comb(n, 2*i) - (-1)**i*comb(n, i))*w**(2*n+i)
        for i in range(1, n+1)
    )
    u = w**3*(4-w)*(3-w)**2
    v = w**3*(3-w)**2
    return Poly(w**12*traces[n] - 2*u**n - v**n - w**n*doubles, w)


def check():
    d, r, w, t = symbols("d r w t")
    assert list(groebner([2*d-3*r, d**2-3*r**2], d, r, domain="QQ")) == [d-3*r/2, r**2]
    lines = 8*81*comb(5, 3)
    planes = 2**5
    cusps = 2*5-2 + 2*6
    weights = (Rational(lines*3, 9), Rational(planes*cusps*3*2, 9), -Rational(planes, 9))
    assert weights == (2160, Rational(1280, 3), -Rational(32, 9))
    correction = Poly(sum(value*w**power for value, power in zip(weights, (15, 17, 18))), w)
    assert gcd(9, planes) == 1
    # w(t) = t^2 times this unit. Only three further Taylor orders enter.
    unit = ((2-2*cos(t))/t**2).series(t, 0, 8).removeO()
    for genus, expected in [(18, Rational(6275, 3)), (19, -Rational(17967, 14))]:
        value = 0
        for (power,), coef in correction.terms():
            order = 2*(genus-1-power)
            if order >= 0:
                expansion = (unit**power).series(t, 0, order+1).removeO()
                value += coef*expansion.coeff(t, order)
        assert value == expected, (genus, value)
    handle = scalar_minus_balanced(w)
    assert all(c.is_Integer for c in (handle + correction).all_coeffs())
    residues = {p: int(9*c) % 9 for (p,), c in handle.terms() if (9*c) % 9}
    assert residues == {17: 3, 18: 5}, residues
    print(f"{lines} lines, {planes} pencils, {3*cusps} cusp hits per pencil, local length 2.")
    print(f"Isolated minus scalar coefficients in genera 16, 18, 19: {weights}.")
    print(f"Ordinary GW corrections in genera 18, 19: {Rational(6275,3)}, {-Rational(17967,14)}.")
    print(f"9 times the independent handle row modulo 9: {residues}. All fractional terms cancel.")


if __name__ == "__main__":
    check()
