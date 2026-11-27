"""Check the two-slice formula with p,s and all four relative step weights free.

Run with SageMath: sage -python checks/quadratic_interpolation.py
"""

from functools import cache
from itertools import combinations

from sage.all import PolynomialRing, QQ


field = PolynomialRing(QQ, names=("p", "s", "d", "b", "e")).fraction_field()
p, s, d, b, e = field.gens()
positive = (field.one(), d)
negative = (b, e)
adams_ratio = (s + p) * (1 - p) / ((s - p) * (1 + p))


def kernel(q, volume, ratio):
    return (1 - ratio) * (1 - q * ratio / volume) / (
        (1 - q * ratio) * (1 - ratio / volume)
    )


def transport(root, x):
    return (x - root) / (x - s * root)


def one_pair(a, b, x):
    return (1 - transport(a, x)) * (1 - s * transport(b, x))


def two_selections(pair, sign, x):
    u, v = pair
    q, volume, shift = (p, s, 1) if sign == 1 else (1 / p, 1 / s, s)
    forward = kernel(q, volume, v / u)
    backward = kernel(q, volume, u / v)
    ru, rv = transport(u, x), transport(v, x)
    mixed = ru * backward + rv * forward if sign == 1 else ru * forward + rv * backward
    square = (forward + backward) * (1 + shift**2 * ru * rv) - 2 * shift * mixed
    elementary = volume * (u - v) * (v - u) / ((volume * u - v) * (volume * v - u))
    adams = (forward + backward - 2 * elementary) * (1 - shift**2 * ru * rv)
    return square, adams


@cache
def one_slice(x):
    plus, adams_plus = two_selections(positive, 1, x)
    minus, adams_minus = two_selections(negative, -1, x)
    return (plus * minus + adams_ratio * adams_plus * adams_minus) / 2


def two_slices(x, y):
    value = transport(1, x) * transport(d, x) * one_slice(y)
    value += transport(b, y) * transport(e, y) * one_slice(x)
    for a, other_a in (positive, positive[::-1]):
        for beta, other_beta in (negative, negative[::-1]):
            weight = kernel(p, s, other_a / a) * kernel(1 / p, 1 / s, other_beta / beta)
            value += (
                weight * transport(beta, y) * transport(other_a, x)
                * one_pair(a, beta, x) * one_pair(other_a, other_beta, y)
            )
    return value


if __name__ == "__main__":
    for x, y in combinations(range(-5, 0), 2):
        assert two_slices(x, y) == two_slices(y, x), (x, y)
    print("All ten off-diagonal differences vanish over Q(p,s,d,b,e).")
