"""Exact genus-two quartic numerical residues in the two fiber weights.

Run with sage -python checks/quartic_weight_one_residues.py.
Bryan--Pandharipande, Theorems 4.1 and 6.6 and Lemma 6.3, give the
partition-basis calculation in prime_five_genus_two. In balanced level
(1,1), its single handle trace is the disconnected genus-two series.
Paper 6, Section 2, gives inversion symmetry and the Laurent degree bound
six for the fixed-source two-weight quartic. All evaluations below take
place in rational functions, before expanding in source genus.
"""

from functools import cache

from sage.all import PowerSeriesRing, QQ, matrix, vector

from prime_five_genus_two import K, W, connected, disconnected, hook_degree, in_w, w


@cache
def quartic(weight):
    terms = [K.one()] + [disconnected(d, weight, 1) for d in range(1, 5)]
    if weight == -1:
        assert terms[1:] == [hook_degree(d) for d in range(1, 5)]
    return in_w(connected(terms)[4])


def odd_row(weight):
    z = QQ(weight)
    return [2 * (z**j + z**(-j)) for j in (1, 3, 5)]


def main():
    weights = (1, 2, 3)
    # Subtract z -> -z. A symmetric pair of odd modes has the same
    # coefficient, and together the two appearances give the factor two.
    values = vector(W, [quartic(z) - quartic(-z) for z in weights])
    modes = matrix(QQ, [odd_row(z) for z in weights]).inverse() * values
    assert sum(a*b for a, b in zip(odd_row(4), modes)) == quartic(4) - quartic(-4)
    assert modes[0] == w**5 * W([1715, -QQ(1215)/4, QQ(2381)/12,
                                 -QQ(2071)/48, QQ(16)/3, -QQ(9)/32])
    assert modes[1] == w**7 * W([QQ(949)/12, -QQ(1157)/48, QQ(57)/16, -QQ(119)/576])
    assert modes[2] == w**9 * W([QQ(5)/48, -QQ(7)/576])

    # Section 8 of Bryan--Pandharipande in genus four, degree two.
    # The scalar trace has the sum of (4 +/- 2 sqrt(w))^3, and the
    # anti-diagonal has two identical hook products. Subtract half the
    # square of the degree-one term in both logarithms.
    scalar_double = w**6 * (128 + 96*w - QQ(1)/2)
    balanced_double = w**6 * (2*(4-w)**3 - QQ(1)/2)
    residual = modes[0] - scalar_double/2
    assert residual == w**5 * W([1715, -QQ(735)/2, QQ(1805)/12,
                                 -QQ(2071)/48, QQ(16)/3, -QQ(9)/32])
    scalar_minus_balanced = quartic(1) - quartic(-1)
    assert scalar_minus_balanced == 4*sum(modes)
    assert quartic(1) == w**4 * W([QQ(5275)/4, 4090, 1155, 42])
    assert quartic(-1) == w**4 * W([QQ(5275)/4, -2770, 2370, -1068, 269, -36, 2])
    elementary = scalar_minus_balanced - 2*balanced_double
    assert elementary == w**5 * W([6860, -1470, 1302, -317, 40, -2])
    inverse = -scalar_minus_balanced + 2*scalar_double
    assert inverse == w**5 * W([-6860, 1470, -918, 269, -36, 2])

    u = PowerSeriesRing(QQ, 'u', default_prec=20).gen()
    substitution = 2*(1-u.cos())
    source = modes[0](substitution)
    residual_source = residual(substitution)
    assert [source[2*h-2] for h in (6, 7, 8)] == [1715, -QQ(3055)/3, QQ(11837)/24]
    assert [residual_source[2*h-2] for h in (6, 7, 8)] == [1715, -QQ(12985)/12, QQ(5725)/12]
    print("Positive weight-one quartic connected residue:", modes[0])
    print("After subtracting half of the genus-four scalar double:", residual)
    print("Raw connected coefficients in source genera 6, 7, 8:",
          [source[2*h-2] for h in (6, 7, 8)])
    print("Numerical quartic exponents in the balanced basis, genus 6 through 11:",
          [elementary[i] for i in range(5, 11)])


if __name__ == '__main__':
    main()
