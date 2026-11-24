"""Exact input checks for the low-divisibility calculation; these do not prove its geometric deductions.

Run: uv run --no-project --with sympy==1.14.0 python checks/partial_integrality.py
"""
import sympy as S


def partitions(n, largest=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n, largest or n), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def hooks(partition):
    for i, row_length in enumerate(partition):
        for j in range(row_length):
            yield row_length - j + sum(length > j for length in partition[i + 1:])


w, z, q, unknown_scalar_cubic = S.symbols("w z q unknown_scalar_cubic")
# w = s^2. The Chebyshev identity computes (2 sin(k u/2))^2.
hook_square = {
    k: S.expand(w * S.chebyshevu(k - 1, S.sqrt(1 - w / 4)) ** 2)
    for k in range(1, 4)
}
elementary = {
    d: S.expand(sum(S.prod(hook_square[h] for h in hooks(partition))
                    for partition in partitions(d)))
    for d in range(1, 4)
}
assert elementary[1] == w
assert elementary[2] == 8 * w**2 - 2 * w**3
assert elementary[3] == 81 * w**3 - 72 * w**4 + 21 * w**5 - 2 * w**6

scalar_series = 1 + w*q + 8*w**2*q**2 + unknown_scalar_cubic*q**3
elementary_series = 1 + sum(value*q**d for d, value in elementary.items())
log_difference = S.series(S.log(scalar_series) - S.log(elementary_series),
                         q, 0, 4).removeO().expand()
assert log_difference.coeff(q, 2) == 2*w**3
assert S.expand(log_difference.coeff(q, 3)
                - (unknown_scalar_cubic - 81*w**3
                   + 70*w**4 - 21*w**5 + 2*w**6)) == 0

# Exact local genus-zero degree-two PT expression, including its sign.
scalar_g0_d2 = 2 / (w**2 * (4 - w))
pt_g0_d2 = S.factor(scalar_g0_d2.subs(w, (1 + z)**2 / z))
assert S.cancel(pt_g0_d2 + 2*z**3 / ((1 + z)**4 * (1 - z)**2)) == 0
assert S.denom(pt_g0_d2).subs(z, 0) == 1

print("Hook enumeration at degree three:", elementary[3])
print("Exact connected cubic difference:", log_difference.coeff(q, 3))
print("Genus-zero degree-two PT expression:", pt_g0_d2)
print("PASS: exact inputs and subtraction in the low-divisibility calculation; no scalar cubic value assumed.")
