"""Exact fixed-point, character, and plethystic checks for the two-circle comparison."""
from collections import Counter
from math import comb, factorial

import sympy as sp


def partitions(size, bound=None):
    if size == 0:
        yield ()
        return
    for first in range(min(size, bound if bound is not None else size), 0, -1):
        for rest in partitions(size - first, first):
            yield (first,) + rest


def plane_data(shape):
    boxes = [(i, j) for i, row in enumerate(shape) for j in range(row)]
    columns = [sum(row > j for row in shape) for j in range(shape[0])]
    hooks = [shape[i] - j + columns[j] - i - 1 for i, j in boxes]
    return hooks, sum(i for i, j in boxes), sum(j for i, j in boxes)


def multiply(left, right):
    result = Counter()
    for i, a in left.items():
        for j, b in right.items():
            result[i+j] += a*b
    return Counter({i: c for i, c in result.items() if c})


def check_split_levels():
    q = sp.Rational(3, 2)
    for b in (-1, 0, 1, 2):
        for degree in range(1, 8):
            for shape in partitions(degree):
                hooks, rows, columns = plane_data(shape)
                assert sum(hooks) == degree + rows + columns
                for delta in (-2, -1, 0, 1, 2):
                    first, second = b+delta, b-delta
                    original = ((-sp.S.One)**(degree*second)
                                * q**(-degree*b-rows*first-columns*second)
                                * sp.prod((1-q**h)**(2*b) for h in hooks))
                    rotated = ((-sp.S.One)**(degree*delta) * q**(delta*(columns-rows))
                               * sp.prod(((1-q**h)*(1-q**(-h)))**b for h in hooks))
                    assert original == rotated


def hook_character(degree, b):
    answer = Counter()
    for shape in partitions(degree):
        term = Counter({0: 1})
        for h in plane_data(shape)[0]:
            for _ in range(b):
                term = multiply(term, {h: -1, -h: -1, 0: 2})
        answer.update(term)
    return answer


def graded_power(rank, count):
    if rank < 0:
        return (-1)**count * comb(-rank, count) if count <= -rank else 0
    return comb(rank+count-1, count)


def check_connected(cutoff):
    for b in (0, 1, 2):
        # Add complete symmetric powers of each signed character as Euler factors.
        product = [Counter({0: 1})] + [Counter() for _ in range(cutoff)]
        for degree in range(1, cutoff+1):
            target = hook_character(degree, b)
            connected = target.copy()
            connected.subtract(product[degree])
            connected = {i: c for i, c in connected.items() if c}
            assert all(c == connected.get(-i) for i, c in connected.items())
            if b == 0:
                assert connected == {0: 1}
            elif degree > 1:
                extremal = b*degree*(degree+1)//2
                assert min(connected) == -extremal
                assert connected[-extremal] == 2 * (-1)**(degree*b)
            for weight, rank in connected.items():
                previous = [entry.copy() for entry in product]
                for n in range(degree, cutoff+1):
                    for count in range(1, n//degree+1):
                        factor = graded_power(rank, count)
                        for power, multiplicity in previous[n-count*degree].items():
                            product[n][power+count*weight] += multiplicity*factor
            assert +product[degree] == +target and -product[degree] == -target


def check_collision():
    q = sp.symbols('q', nonzero=True)
    y = q+1/q
    e = 2-y
    e2 = e.subs(q, q**2)
    plus = 2*e*e2 - (e**2+e2)/2
    minus = y*e*e2 - (e**2-e2)/2
    assert sp.expand(plus - (-2*e+8*e**2-2*e**3)) == 0
    assert sp.expand(plus-minus-e2*(e**2-1)) == 0
    # Independent sheet trace via the two blowup charts of the ordered Hilbert square.
    ordered = e*(y**2-4)**2*sum(1/((1-w**2)*(1+1/w)) for w in (q, 1/q))
    assert sp.cancel(ordered-e2-(plus-minus)) == 0
    v = sp.symbols('v')
    strict = 2*(y-2)*(y-v)**2
    proposed_wall = 2*(y-2)*(2*v*y+v**2-12)
    assert sp.expand(plus-strict).coeff(q, -2) == 4*v
    assert sp.expand(plus-strict-proposed_wall-2*(y-2)*(9-2*v**2)) == 0


def tableau_degrees(shape, variable):
    """Enumerate standard tableaux by adding the next number to a row."""
    terms = Counter()

    def extend(lengths, previous_row, filled, major):
        if lengths == shape:
            terms[major] += 1
            return
        for row, bound in enumerate(shape):
            if lengths[row] == bound or (row and lengths[row] == lengths[row-1]):
                continue
            next_lengths = list(lengths)
            next_lengths[row] += 1
            descent = filled if row > previous_row else 0
            extend(tuple(next_lengths), row, filled + 1, major + descent)

    extend((0,) * len(shape), 0, 0, 0)
    return sp.Poly.from_dict({(k,): a for k, a in terms.items()}, variable, domain=sp.ZZ)


def check_pairing_map():
    q = sp.Symbol("q")
    for m in range(1, 8):
        degree_bound = m * (m-1) // 2
        shapes = list(partitions(m))

        def reverse_weights(polynomial):
            terms = {(degree_bound-k,): a for (k,), a in polynomial.terms()}
            return sp.Poly.from_dict(terms, q, domain=sp.ZZ)

        def one_minus(power):
            return sp.Poly(1 - q**power, q, domain=sp.ZZ)

        numerator = sp.prod(one_minus(j) for j in range(1, m+1))
        shifted_fiber = sp.Poly(0, q, domain=sp.ZZ)
        shifted_average = sp.Poly(0, q, domain=sp.ZZ)
        for shape in shapes:
            fake_degree = tableau_degrees(shape, q)
            shifted_fiber += fake_degree * reverse_weights(fake_degree)
            # Independently trace the Koszul zero-fiber algebra on each
            # permutation class, then average the product of its two traces.
            trace = numerator.exquo(sp.prod(one_minus(j) for j in shape))
            stabilizer = sp.prod(j**count * factorial(count) for j, count in Counter(shape).items())
            class_size = factorial(m) // stabilizer
            shifted_average += class_size * trace * reverse_weights(trace)
        assert factorial(m) * shifted_fiber == shifted_average
        assert shifted_fiber.eval(1) == factorial(m)
        if m == 3:
            assert shifted_fiber.as_expr() == q**4 + 4*q**3 + q**2
        for b in range(4):
            state = hook_character(m, b)
            low = min(state)
            shifted_state = sp.Poly.from_dict({(j-low,): v for j, v in state.items()}, q, domain=sp.ZZ)
            assert shifted_fiber.gcd(shifted_state).as_expr() == (2 if m == 2 else 1)


def check_smooth_branches():
    z, a = sp.symbols('z a')
    x = z+2+1/z
    for sign in (1, -1):
        jacobian = (1+z)**4*(1+sign*z)**4/z**4
        predicted = x**4-4*(1-sign)*x**3+8*(1-sign)*x**2
        assert sp.cancel(jacobian-predicted) == 0
        branch = (1+z)**2*(1+sign*z)**4/z**3
        assert sp.cancel(jacobian*z/(1+z)**2-branch) == 0
        assert sp.expand(-2*branch).coeff(z, -2) == (-12 if sign == 1 else 4)
    # Evaluate in the exact cyclotomic ring via the two-dimensional integral basis.
    def reduce_cubic(poly):
        return sp.rem(poly, a**2+a+1, a).expand()
    for g in (3, 4):
        base = (a-a**2)*(1+z)**(2*g-1)*(1+a*z)**(2*g+2)*(1+a**2*z)**(2*g-1)
        factored = (a-a**2)*(1+z**3)**(2*g-1)*(1+a*z)**3
        assert reduce_cubic(sp.Poly(base-factored, z).as_expr()) == 0
    # External two-line spin Euler at independently variable volume.
    Q, K = sp.symbols('Q K', nonzero=True)
    assert sp.expand((Q-1/Q)*(1/(K*Q)-K*Q)
                     - (1-Q**2)*(1-(K*Q)**2)/(-K*Q**2)) == 0


if __name__ == '__main__':
    check_split_levels()
    check_connected(6)
    check_collision()
    check_smooth_branches()
    check_pairing_map()
    print('Split levels through degree 7, connected Euler products through degree 6,')
    print('ordered and unordered collisions, smooth characters and pairing maps through degree 7 verified.')
