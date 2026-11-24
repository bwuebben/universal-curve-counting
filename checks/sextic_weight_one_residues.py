"""Exact genus-two sextic Laurent residues and their scalar proper-root terms.

Run with sage -python checks/sextic_weight_one_residues.py.
Bryan--Pandharipande's partition algebra computes the disconnected numbers.
Papers 2 and 6 specify the weight signs, genus variable and proper divisions.
The universal interpretation of the numerical subtraction uses the degree-six comparison.
"""
from concurrent.futures import ProcessPoolExecutor
from math import comb

from sage.all import Partitions, PowerSeriesRing, QQ, matrix, prod, vector

from prime_five_genus_two import K, W, connected, disconnected, hook_degree, in_w, r, w


def genus_two_logs(z):
    terms = [K.one()]
    for degree in range(1, 7):
        # At the antidiagonal the degree-six transposition operator does
        # not generate the center. Use the independent hook formula there.
        residue = hook_degree(degree) if z == -1 else disconnected(degree, z, 1)
        if z == -1 and degree < 6:
            assert residue == disconnected(degree, z, 1)
        terms.append(residue)
    return [in_w(term).list() for term in connected(terms)[1:]]


def odd_row(weight, exponents):
    z = QQ(weight)
    return [2 * (z**j + z**(-j)) for j in exponents]


def odd_modes(differences, degree):
    # n_d = binomial(d,2) for a genus-two core. Inversion symmetry
    # pairs j with -j. Replacing z by -z kills every even weight.
    exponents = range(1, comb(degree, 2) + 1, 2)
    count = len(exponents)
    weights = range(1, count + 1)
    values = vector(W, [differences[z][degree - 1] for z in weights])
    modes = matrix(QQ, [odd_row(z, exponents) for z in weights]).inverse() * values
    for z in range(count + 1, 10):
        predicted = sum(a*b for a, b in zip(odd_row(z, exponents), modes))
        assert predicted == differences[z][degree - 1]
    return modes


def scalar_double(genus):
    n = genus - 1
    # (4 + 2 sqrt(w))^n + (4 - 2 sqrt(w))^n, followed by the log.
    trace = 2*sum(comb(n, 2*j)*4**(n - 2*j)*(4*w)**j for j in range(n//2 + 1))
    return w**(2*n) * (trace - QQ(1)/2)


def genus_four_scalar_triple():
    # Paper 6, Section 4. The roots of X^3 - A X^2 + B X/3 - C/3
    # are the normalized eigenvalues of the genus-two degree-three handle.
    A = w**2 + 50*w + 81
    auxiliary = w**2 + 54*w - 243
    B = 4*(w - 6)*auxiliary
    C = 4*(w - 4)*(w - 3)**2*auxiliary
    trace_cubed = A**3 - A*B + C
    linear = w**3
    quadratic = scalar_double(4) + linear**2/2
    return w**9*trace_cubed - linear*quadratic + linear**3/3


def balanced_log(genus, degree):
    # The antidiagonal disconnected residue is the hook sum with exponent
    # genus-1 on each squared sine. Keep all lower powers in the logarithm.
    terms = [K.one()]
    for d in range(1, degree + 1):
        total = K.zero()
        for shape in Partitions(d):
            hooks = [hook for row in shape.hook_lengths() for hook in row]
            total += prod(2-r**h-r**(-h) for h in hooks)**(genus-1)
        terms.append(total)
    return in_w(connected(terms)[degree])


def main():
    weights = [sign*z for z in range(1, 10) for sign in (1, -1)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        logs = dict(zip(weights, pool.map(genus_two_logs, weights)))
    differences = {z: [W(p)-W(m) for p, m in zip(logs[z], logs[-z])]
                   for z in range(1, 10)}
    assert odd_modes(differences, 2)[0] == w**3/2
    assert odd_modes(differences, 3)[0] == 30*w**4 - 5*w**5 + QQ(4)/9*w**6
    modes = odd_modes(differences, 6)
    for j, mode in zip(range(1, 16, 2), modes):
        assert mode.valuation() >= 6+j
    # The outermost mode is the geometric outermost-weight scalar generator.
    outer_denominator = 2*3**2*4**3*5**4*6**5
    assert modes[-1] == w**21/outer_denominator
    raw = modes[0]
    expected_raw = W([
        8651506, 15352931, QQ(102753031)/6, -QQ(3921065)/6,
        QQ(733684957)/180, -QQ(64114396)/27, QQ(716237737)/600,
        -QQ(10003200103)/21600, QQ(33773279609)/243000, -QQ(862549187)/27000,
        QQ(357544757)/64800, -QQ(21044614)/30375, QQ(214897)/3600,
        -QQ(142)/45, QQ(25)/324,
    ])
    assert raw == w**7*expected_raw
    double = 30*scalar_double(5) - 5*scalar_double(6) + QQ(4)/9*scalar_double(7)
    triple = genus_four_scalar_triple()/2
    assert double == w**8 * W([15345, 23040, -QQ(18555)/2, -25600, QQ(1322)/3,
                               QQ(40960)/3, QQ(10240)/3, QQ(512)/9])
    assert triple == w**9 * W([QQ(140870)/3, 418155, QQ(735903)/2, 69322,
                               QQ(7527)/2, 75, QQ(1)/2])
    residual = raw - double - triple
    # the degree-six calculation counts 288 central amplitude lines. The normal-parameter
    # graph has self-intersection -1 and the scalar test has weight 3.
    isolated_triple_correction = QQ(288)/3
    simple_hits = residual - isolated_triple_correction*w**12/2
    assert simple_hits[12] == -QQ(65999284)/27
    assert all((outer_denominator**2*c).denominator() == 1 for c in simple_hits)
    # A cubic factor at q^2 and three quadratic factors at q^3 enter
    # the elementary degree-six conversion. Their logarithms are taken
    # at antidiagonal weights in the intermediate curve, before subtraction.
    scalar_minus_balanced = differences[1][5]
    assert scalar_minus_balanced == 4*sum(modes)
    from_quadratic = 2*balanced_log(4, 3)
    from_cubic = 120*balanced_log(5, 2) - 20*balanced_log(6, 2) + 2*balanced_log(7, 2)
    elementary = scalar_minus_balanced - from_quadratic - from_cubic
    assert elementary == w**7 * W([
        34606024, 61350344, 125323022, -36417600, 62126672, -41880636,
        22703802, -9395788, 2975620, -715462, 128364, -16652, 1476, -80, 2,
    ])
    print('Coefficient rows in w, increasing powers 7 through 21')
    for title, row in [('raw', raw), ('scalar double', double),
                       ('scalar triple', triple), ('residual', residual)]:
        print(title, [row[i] for i in range(7, 22)])
    u = PowerSeriesRing(QQ, 'u', default_prec=24).gen()
    sine_square = 2*(1-u.cos())
    predicted_source = residual(sine_square)
    assert [predicted_source[2*g-2] for g in range(8, 13)] == [
        8651506, QQ(61745245)/6, QQ(2973727027)/360,
        -QQ(4662673193)/432, QQ(15680715194051)/1814400,
    ]
    print('Coefficient rows in source genera 8 through 12')
    for title, row in [('raw', raw), ('scalar double', double),
                       ('scalar triple', triple), ('residual', residual)]:
        series = row(sine_square)
        print(title, [series[2*g-2] for g in range(8, 13)])
    print('Scalar-to-elementary sextic exponents in output genera 8 through 22')
    print([elementary[i] for i in range(7, 22)])


if __name__ == '__main__':
    main()
