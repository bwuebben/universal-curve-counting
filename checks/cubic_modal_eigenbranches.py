# /// script
# dependencies = ["sympy==1.14.0"]
# ///
"""Compute the genus-ten Laurent coefficients by lifting Frobenius characters.

Run with uv run checks/cubic_modal_eigenbranches.py.
Weights are t1=-l^2, t2=1, and v=(r-1)/(r+1).
The calculation takes place in Q(l)[v]/(v^7).
"""

from sympy import QQ, Rational, expand
from sympy.polys.rings import ring
from sympy.polys.ring_series import rs_series_inversion, rs_trunc

field = QQ.frac_field("l")
ring_v, v = ring("v", field)
l = ring_v(field.gens[0])
one = ring_v.one
k = -l**2
t1, t2 = k, one


def trunc(p):
    return rs_trunc(p, v, 7)


def inv(p):
    return rs_series_inversion(trunc(p), v, 7)


def lift(seed, polynomial, derivative):
    root = seed
    # Each Newton step doubles the known precision, starting at one.
    for _ in range(3):
        root = trunc(root - trunc(polynomial(root)) * inv(derivative(root)))
    assert not trunc(polynomial(root))
    return root


a = -(t1 + t2) * v
b = trunc(-12 * (t1 + t2) * v * inv(3 + v**2))
cap_factor = trunc(4 * v**2 * inv(3 + v**2))


def cubic(x):
    return (x - b) * (x**2 - a * x + 3 * k) + 6 * k * x


def cubic_derivative(x):
    return 3 * x**2 - 2 * (a + b) * x + a * b + 9 * k


def cubic_cap(t, x):
    y = (x**2 - a * x + 3 * k) / 3
    return trunc(t**3 + t**2 * x * v + t * y * cap_factor)


def cubic_square_trace():
    total = ring_v.zero
    for seed in (ring_v.zero, 3 * l, -3 * l):
        x = lift(seed, cubic, cubic_derivative)
        # The inverse Frobenius norm of the primitive idempotent.
        numerator = -6 * k**3 * cubic_derivative(x)
        denominator = x**2 - (a + b) * x + a * b + 6 * k
        denominator *= cubic_cap(t1, x) * cubic_cap(t2, x)
        normalized = trunc(numerator * inv(denominator))
        total += trunc(normalized**2)
    return total


def quadratic_square_trace():
    total = ring_v.zero
    for seed in (l, -l):
        x = lift(seed, lambda x: x**2 - a * x + k, lambda x: 2 * x - a)
        numerator = 2 * k**2 * (2 * x - a)
        denominator = (x - a) * (t1**2 + t1 * x * v) * (t2**2 + t2 * x * v)
        normalized = trunc(numerator * inv(denominator))
        total += trunc(normalized**2)
    return total


def w3_modes(series):
    assert all(not series.get((i,)) for i in (1, 3, 5))
    f2, f4, f6 = (series.get((i,), field.zero) for i in (2, 4, 6))
    w3 = expand(((-f2 + 2 * f4 - f6) / 64).as_expr())
    l_symbol = field.gens[0].as_expr()
    out = {}
    if w3 == 0:
        return {}
    for term in w3.as_ordered_terms():
        coefficient, power = term.as_coeff_exponent(l_symbol)
        assert power.is_Integer and power % 2 == 0
        j = int(power) // 2
        out[j] = coefficient * (-1 if j % 2 else 1)
    assert all(out[j] == out[-j] for j in out)
    return dict(sorted(out.items()))


if __name__ == "__main__":
    raw = cubic_square_trace()
    quadratic = quadratic_square_trace()
    # [q]Z=w^2, so the full q^3 logarithm is w^6 times this normalized expression.
    connected = raw - quadratic + one / 3
    for name, series in (("degree-three", raw), ("degree-two", quadratic), ("connected", connected)):
        modes = w3_modes(series)
        print(name, modes)
    assert not w3_modes(quadratic)
    assert modes[1] == Rational(674, 3) and modes[3] == Rational(506, 3)
    scalar = sum(modes.values())
    antidiagonal = sum((-1 if j % 2 else 1) * c for j, c in modes.items())
    assert scalar == Rational(292, 3) and antidiagonal == -1476
    assert scalar - antidiagonal == Rational(4720, 3)
    print("scalar minus antidiagonal =", scalar - antidiagonal)
