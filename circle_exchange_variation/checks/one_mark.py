"""Check the one-mark Macdonald calculation. Run with sage -python checks/one_mark.py."""

from math import prod

from sage.all import PolynomialRing, QQ

from fixed_vertex import partitions, polarization, vertex


def diagram(shape):
    boxes = [(i, j) for j, width in enumerate(shape) for i in range(width)]
    removable = [(i, j) for i, j in boxes
                 if (i + 1, j) not in boxes and (i, j + 1) not in boxes]
    addable = [(width, j) for j, width in enumerate(shape)
               if j == 0 or width < shape[j - 1]] + [(0, len(shape))]
    return boxes, removable, addable


def corner_terms(shape, u, p, c):
    _, corners, addable = diagram(shape)
    v, s = c/u, p/c
    weight = lambda h: u**h[0] * v**h[1]
    R = lambda z: (1 - s*z)/(1 - z)
    return {(i, j): R(u)*R(v)
            * prod(R(weight((i, j))/weight(h)) for h in corners if h != (i, j))
            / prod(R(c*weight((i, j))/weight(h)) for h in addable)
            for i, j in corners}


def exchanged_terms(shape, u, p, c):
    boxes, corners, _ = diagram(shape)
    s = p/c
    result = {}
    for corner in corners:
        heights = [int(h == corner) for h in boxes]
        term = 1/s
        for (i, j, k), multiplicity in polarization(boxes, heights).items():
            z = u**i * p**j * s**k
            term *= ((1 - s*z)/(1 - z))**multiplicity
        result[corner] = term * (1 - p)/(1 - c)
    return result


def wick_terms(shape, u, p, c):
    """Use the residues of the two currents before finite differencing."""
    boxes, _, _ = diagram(shape)
    steps = range(-len(shape), shape[0])
    lower = {i: max(0, -i) + sum(k-j == i for k, j in boxes)
             for i in range(steps.start, steps.stop+1)}
    z = {i: u**i * c**lower[i] for i in lower}
    minus = [i for i in steps if lower[i+1] == lower[i]]
    plus = [i for i in steps if lower[i+1] < lower[i]]
    s = p/c

    def inside(i, alpha):
        if alpha < i:
            return 0
        h = s*z[alpha]
        return ((s-1) * prod((1-h/(s*z[j]))/(1-h/z[j]) for j in minus if j < i)
                * prod((1-h/z[j])/(1-h/(s*z[j])) for j in plus if j >= i and j != alpha))

    def outside(i, beta):
        if beta >= i:
            return 0
        w = z[beta]
        return (-(s-1)*s**(i-1)
                * prod((1-z[j]/(s*w))/(1-z[j]/w) for j in minus if j < i and j != beta)
                * prod((1-z[j]/(c*w))/(1-z[j]/(p*w)) for j in plus if j >= i))

    result = {}
    for beta in minus:
        for alpha in plus:
            if alpha <= beta:
                continue
            h = s*z[alpha]/z[beta]
            kernel = (1-h)*(1-h/(p*s))/((1-h/p)*(1-h/s))
            value = -s*kernel/(s-1)**2 * sum(
                (inside(i, alpha)-inside(i+1, alpha)) *
                (outside(i, beta)-s*outside(i-1, beta))
                for i in range(beta+1, alpha+1))
            if alpha == beta+1:
                result[alpha+lower[alpha]-1, lower[alpha]-1] = value
            else:
                assert value == 0, (shape, beta, alpha, value)
    return result


def telescope():
    field = PolynomialRing(QQ, names=("s", "a", "u", "t")).fraction_field()
    s, a, u, t = field.gens()
    forward = (a-t)/(a-s*t)
    backward = (a-s*u)/(a-u)
    assert (1-forward)*(backward-s) == (s-1)*t/(t-u)*(forward*backward-1)
    # At the birth and death use the cross-multiplied identity, without a pole.
    assert (1-forward.subs({a:u})) == ((s-1)*t/(t-u))*forward.subs({a:u})
    assert backward.subs({a:t})-s == -(s-1)*t/(t-u)


def genus_two_double():
    """Compare the full q^-2 t coefficient with the independent double formula."""
    field = PolynomialRing(QQ, names=("u", "K", "c")).fraction_field()
    u, K, c = field.gens()
    # q^3 times the two empty contributions, with q=u.
    common = (1-u)*(1-u/c)
    row = c*common*(1-u**2)*(1-u**2/c)
    column = c**2*common*(1-u**2/c)*(1-u**2/c**2)
    marked = (1-c)*(row*sum(corner_terms((2,),u,K,c).values())
                    + column*sum(corner_terms((1,1),u,K,c).values()))
    result = marked.derivative(u).subs({u:0})
    assert result == (c-1)*(1+K+c**2/K+c**2)
    # c=(1-delta)^-1, K=1+epsilon, epsilon^2=0.
    first = result.derivative(c).subs({c:1})
    second = first + result.derivative(c,2).subs({c:1})/2
    assert (first.subs({K:1}),first.derivative(K).subs({K:1})) == (4,0)
    assert (second.subs({K:1}),second.derivative(K).subs({K:1})) == (8,-2)


if __name__ == "__main__":
    telescope()
    genus_two_double()
    shapes = [s for m in range(1, 8) for s in partitions(m, m)]
    shapes += [(3, 3, 3), (4, 3, 2, 1), (5, 4, 3), (5, 5, 4, 2)]
    for p, c, u in [(QQ(2), QQ(3), QQ(5)/7), (QQ(5)/3, QQ(7)/4, QQ(2)/5)]:
        for shape in shapes:
            assert wick_terms(shape,u,p,c) == corner_terms(shape,u,p,c) == exchanged_terms(shape,u,p,c), shape
    # These comparisons retain both independent parameters of the rational functions
    # in u and p, with c a fixed generic rational; t-degree is one at all u-gradings.
    field = PolynomialRing(QQ, names=("u", "p")).fraction_field()
    u, p = field.gens()
    c = QQ(3)
    for shape in [(1,), (3,1), (2,2), (3,2,1), (3,3,2)]:
        assert wick_terms(shape,u,p,c) == corner_terms(shape,u,p,c) == exchanged_terms(shape,u,p,c), shape
    # Independent original finite-height computation including the empty hook product.
    p, c = QQ(2), QQ(3)
    for shape in [(4,2), (3,2,1), (3,3,2), (3,3,3)]:
        original = vertex(shape,p,c,6)
        rotated = vertex(shape,p,p/c,6)
        assert all(original[n][1] == rotated[1][n] for n in range(7)), shape
    print(f"Cleared telescoping, its endpoints, and the genus-two nonlinear check verified. "
          f"Tested all ordered residue pairs and "
          f"each corner on {len(shapes)} diagrams at two exact triples, five diagrams as rational "
          "functions of u,p, and four original height expansions through Euler excess six.")
