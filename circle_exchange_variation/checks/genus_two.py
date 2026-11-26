"""Genus-two double, Picard-valued first variation in the two circle chambers."""
from itertools import permutations
from math import factorial

import sympy as sp


def verify_ribbon_prediction():
    t, delta, theta = sp.symbols('t delta theta')
    kappa = 1-delta
    # Chern roots of the two derived cohomologies have power sums -theta, -4theta.
    first_order = 0
    for p in (t, kappa/t):
        bare = kappa**-2 * (p**2-kappa)/(p**2-1)
        pole = sum(mult/((w-1)*(kappa/w-1)) for w, mult in ((p, 1), (p**2, 4)))
        first_order += sp.diff(bare, delta).subs(delta, 0) - theta*pole.subs(delta, 0)
    expected = 3 + 2*theta*(t/(1-t)**2 + 4*t**2/(1-t**2)**2)
    assert sp.cancel(first_order-expected) == 0
    coefficients = sp.series(expected, t, 0, 21).removeO().expand()
    for k in range(1, 21):
        assert coefficients.coeff(t, k) == 2*k*(1+2*(k%2 == 0))*theta


def verify_rotated_weights():
    u, v, tangent, functions = sp.symbols('u v tangent functions')
    kappa = u*v
    minimum = -(u+v+u**2+v/u)
    # Compare the actual framed complex with the weights used at each fixed component.
    global_functions = -1
    def tangent_character(root, tip, dual_root, dual_tip, ratio, dual_ratio):
        return sp.expand(root+tip/u + kappa*(dual_root+u*dual_tip)
                         - (2*global_functions+ratio/u+u*dual_ratio)*(1-u)*(1-v))
    tip_only = tangent_character(-1, -1+tangent, -1, -1-functions,
                                 -1+tangent, -1-functions)
    common = tangent_character(-1+tangent, -1+tangent, -1-functions, -1-functions, -1, -1)
    assert sp.expand(tip_only-minimum-tangent*(1+v/u-v)-functions*u*(1-u-v)) == 0
    assert sp.expand(common-minimum-tangent*(1+1/u)+kappa*functions*(1+u)) == 0
    for length in range(1, 21):
        for weights in ((length-1, -length-1, length-1, -length-1),
                        (-1, length-1, -length-1, -1)):
            neg2, neg1, pos1, pos2 = weights
            assert length + pos1 + 2*pos2 == -3
            assert length - neg1 - 2*neg2 == 3
    canonical, q = sp.symbols('canonical q')
    row = q*(1-q)*(1-q**-2)**2*(1-q**2)*(1-q*canonical)/(1-q**2*canonical)
    total = row+row.subs(q, 1/q)
    y = q+1/q
    assert sp.factor(total.subs(canonical, 1)-2*(y-2)*(y**2-4)) == 0
    assert sp.factor(sp.diff(total, canonical).subs(canonical, 1)) == 0


def cycles(permutation):
    unused = set(range(len(permutation)))
    answer = []
    while unused:
        length, pt = 0, next(iter(unused))
        while pt in unused:
            unused.remove(pt)
            pt = permutation[pt]
            length += 1
        answer.append(length)
    return answer


def verify_adams_push():
    # Two even characters of weights 1,2 and one odd of weight 3. The test sum is
    # an independent character-level check of the integral symmetric-functor proof.
    a = sp.Integer(2)
    for d in range(1, 8):
        total = 0
        for perm in permutations(range(d)):
            lengths = cycles(perm)
            det_one_minus_standard = d if len(lengths) == 1 else 0
            total += det_one_minus_standard * sp.prod(a**h+a**(2*h)-a**(3*h) for h in lengths)
        assert total/factorial(d) == a**d+a**(2*d)-a**(3*d)
    j = sp.symbols('j')
    # In the K-group of the symmetric square, (J-1)^3 vanishes.
    assert sp.factor(j*((j-1)*(2+j+1/j)-4*(j-1)) - (j-1)**3) == 0


if __name__ == '__main__':
    verify_ribbon_prediction()
    verify_rotated_weights()
    verify_adams_push()
    print('Verified the genus-two ribbon and fixed-divisor weights, first variations')
    print('through exchanged grading 20, and Adams characters through length seven.')
