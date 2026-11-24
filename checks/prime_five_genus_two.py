"""Connected degree-five scalar correction over a genus-two core.

Run with sage -python checks/prime_five_genus_two.py. Computations use
exact rational functions in r=exp(iu). No interpolation in genus is used.
"""
from collections import Counter
from math import factorial, prod

from sage.all import QQ, Partitions, PolynomialRing, matrix, vector


R = PolynomialRing(QQ, "r")
K = R.fraction_field()
r = K.gen()
P = PolynomialRing(K, "X")
X = P.gen()
W = PolynomialRing(QQ, "w")
w = W.gen()


def partition_norm(part):
    return prod(size**count * factorial(count) for size, count in Counter(part).items())


def partition_key(counts):
    return tuple(sorted(counts.elements(), reverse=True))


def disconnected(degree, t1, t2):
    """Bryan–Pandharipande cap and balanced handle in the partition basis."""
    parts = [tuple(p) for p in reversed(list(Partitions(degree)))]
    index = {p: i for i, p in enumerate(parts)}
    size = len(parts)
    operator = matrix(K, size, size)
    for column, part in enumerate(parts):
        counts = Counter(part)
        operator[column, column] = -(t1 + t2) / QQ(2) * sum(
            a * count * (a * (r**a + 1) / (r**a - 1) - (r + 1) / (r - 1))
            for a, count in counts.items()
        )
        for a in range(1, degree):
            for b in range(1, degree - a + 1):
                if counts[a+b]:
                    split = counts.copy()
                    split[a+b] -= 1
                    split[a] += 1
                    split[b] += 1
                    dest = partition_key(split)
                    operator[index[dest], column] -= (
                        QQ(t1*t2*(a+b)*counts[a+b]*partition_norm(dest))
                        / (2*partition_norm(part))
                    )
                pairs = counts[a]*(counts[b] - (a == b))
                if pairs:
                    joined = counts.copy()
                    joined[a] -= 1
                    joined[b] -= 1
                    joined[a+b] += 1
                    dest = partition_key(joined)
                    operator[index[dest], column] += (
                        QQ(a*b*pairs*partition_norm(dest)) / (2*partition_norm(part))
                    )

    # The transposition operator generates the Frobenius algebra here.
    # Convert each partition class into its polynomial in that generator.
    unit = vector(K, [int(i == 0) for i in range(size)])
    columns = [unit]
    for _ in range(size-1):
        columns.append(operator * columns[-1])
    krylov = matrix(K, columns).transpose()
    inverse = krylov.inverse()  # Also certifies cyclicity for this degree and weights.
    characteristic = P(operator.charpoly())
    classes = [P(list(column)) for column in inverse.columns()]
    handle = P.zero()
    cap1, cap2 = P.zero(), P.zero()
    for part, element in zip(parts, classes):
        metric_inverse = (-1)**len(part) * partition_norm(part) * (t1*t2)**len(part)
        handle += (metric_inverse * element**2) % characteristic
        factor = (-1)**degree / prod(1-r**(-i) for i in part)
        cap1 += t1**len(part)*factor*element
        cap2 += t2**len(part)*factor*element
    balanced = (r**degree * handle * (cap1*cap2).inverse_mod(characteristic)) % characteristic
    powers = [K(size)]
    op_power = operator
    for _ in range(1, size):
        powers.append(op_power.trace())
        op_power *= operator
    return sum(coef*trace for coef, trace in zip(balanced.list(), powers))


def connected(disconnected_terms):
    # q F'(q) = F(q) q (log F)'(q).
    coefficients = [K.zero()]
    for d in range(1, len(disconnected_terms)):
        correction = sum(k*coefficients[k]*disconnected_terms[d-k] for k in range(1, d))
        coefficients.append(disconnected_terms[d]-correction/QQ(d))
    return coefficients


def hook_degree(d):
    result = K.zero()
    for shape in Partitions(d):
        hooks = [length for row in shape.hook_lengths() for length in row]
        result += prod(2-r**length-r**(-length) for length in hooks)
    return result


def in_w(value):
    assert value.parent() is K
    denominator = value.denominator()
    assert len(denominator.dict()) == 1
    shift = denominator.degree()
    coeff = value.numerator() / denominator.leading_coefficient()
    assert all(coeff[shift-i] == coeff[shift+i] for i in range(1,shift+1))
    assert coeff.degree() <= 2*shift
    previous, current = W(2), 2-w
    result = W(coeff[shift])
    for i in range(1,shift+1):
        result += coeff[shift+i]*current
        previous, current = current, (2-w)*current-previous
    assert result(2-r-1/r) == value
    return result


if __name__ == "__main__":
    scalar = [K.one()]
    anti = [K.one()]
    for d in range(1,6):
        scalar.append(disconnected(d, 1, 1))
        anti.append(hook_degree(d))
        assert disconnected(d, 1, -1) == anti[-1]
    scalar_log, anti_log = connected(scalar), connected(anti)
    differences = [in_w(scalar_log[d]-anti_log[d]) for d in (2,3,5)]
    assert differences[:2] == [2*w**3, 120*w**4-20*w**5+2*w**6]
    expected = w**6 * W([444088, 142232, 340864, -156542, 75396,
                            -23737, 5054, -696, 56, -2])
    assert differences[-1] == expected
    print("Degree-five difference in genus two =", differences[-1])
