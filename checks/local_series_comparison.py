"""Exact degree-two BP matrix comparison; no truncation in the genus variable.

Run: uv run --no-project --with sympy==1.14.0 python checks/local_series_comparison.py
Matrices: Bryan--Pandharipande arXiv:math/0411037v3, Section 8, p. 40.
s = 2 sin(u/2), c = cos(u/2), c^2 = 1 - s^2/4.
"""
import sympy as S

s, c = S.symbols('s c', nonzero=True)
t1, t2 = S.symbols('t1 t2', nonzero=True)
G = S.Matrix([[4*(t1*t2)**2, -2*(t1*t2)**2*(t1+t2)*s/(2*c)],
              [-2*(t1*t2)*(t1+t2)*s/(2*c),
               4*(t1*t2)**2+2*t1*t2*(t1+t2)**2*s**2/(4*c**2)]])
A = S.Matrix([[t1**2/s**2, -t1**2*t2/(2*s*c)],
              [-t1/(2*s*c), t1*(t1+t2)/(4*c**2)+t1**2/s**2]])
Abar = A.xreplace({t1:t2, t2:t1})

def degree_two(weights):
    raw=S.factor(S.trace(G.subs(weights)*A.subs(weights).inv()*Abar.subs(weights).inv()))
    return S.factor(raw.subs(c**2, 1-s**2/4))

scalar=degree_two({t1:1,t2:1})
anti=degree_two({t1:1,t2:-1})
assert S.simplify(scalar-8*s**4)==0
assert S.simplify(anti-(8*s**4-2*s**6))==0
assert S.simplify(scalar-anti-2*s**6)==0
print('Scalar genus-two degree-two:',scalar)
print('Antidiagonal genus-two degree-two:',anti)
print('Exact difference:',S.factor(scalar-anti))

for g in range(2,21):
    delta=S.expand(s**(4*g-4)*((4-2*s)**(g-1)+(4+2*s)**(g-1)-2*(4-s**2)**(g-1)))
    correction=sum(2*4**(g-1-j)*(S.binomial(g-1,2*j)-(-1)**j*S.binomial(g-1,j))*s**(4*g-4+2*j)
                   for j in range(1,g))
    assert S.expand(delta-correction)==0
print('General binomial identity replayed at g=2,...,20; proof is binomial expansion.')
