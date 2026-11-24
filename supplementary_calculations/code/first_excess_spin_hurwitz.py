"""First-excess spin degrees from collision counts, with finite group cross-checks."""

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import permutations
from math import comb, factorial, prod

from etale_elementary_exponents import connected_counts, elementary_exponents
from etale_spin_characters import compose, partitions, permutation_sign


def divide_series(numerator, denominator):
    """Truncated ordinary power series division with unit constant denominator."""
    assert denominator[0] == 1
    quotient = []
    for m, coefficient in enumerate(numerator):
        quotient.append(coefficient - sum(denominator[m-j] * q for j, q in enumerate(quotient)))
    return quotient


@cache
def spin_series(genus, limit, mixed):
    """Connected nodal, triple, finite and total degrees with base sign +1."""
    denominator, nodal, triple, hooks = [1], [0], [0], [0]
    for m in range(1, limit + 1):
        b = triple_sum = hook_sum = 0
        for shape in partitions(m):
            transpose = tuple(sum(row > col for row in shape) for col in range(shape[0]))
            if mixed and shape != transpose:
                continue
            cells = [(i, j) for i, row in enumerate(shape) for j in range(row)]
            lengths = [shape[i] + transpose[j] - i - j - 1 for i, j in cells]
            content_square = sum((j - i)**2 for i, j in cells)
            hook_square = sum(h**2 for h in lengths)
            assert hook_square == m**2 + content_square
            weight = prod(lengths)**(2 * genus - 2)
            b += weight
            triple_sum += weight * (content_square - comb(m, 2))
            hook_sum += weight * hook_square
        denominator.append(b)
        nodal.append(comb(m, 2) * b)
        triple.append(triple_sum)
        hooks.append(hook_sum)
    nodal = divide_series(nodal, denominator)
    triple = divide_series(triple, denominator)
    hooks = divide_series(hooks, denominator)
    finite = [-Fraction(genus - 1, 12) * (3*a + b) for a, b in zip(nodal, triple)]
    total = [-Fraction(genus - 1, 12) * h for h in hooks]
    etale = connected_counts(genus, limit, mixed)
    for m in range(1, limit + 1):
        assert total[m] == finite[m] - Fraction(m * (genus - 1), 12) * etale[m]
    return nodal, triple, finite, total


def first_exponents(genus, limit, mixed, epsilon):
    minimal = elementary_exponents(genus, limit, mixed, epsilon)
    finite = spin_series(genus, limit, mixed)[2]
    exponents = [0]
    for m in range(1, limit + 1):
        exponent = epsilon**m * finite[m]
        for d in range(1, m):
            if m % d == 0:
                k = m // d
                elementary_finite = spin_series(1 + d*(genus-1), k, False)[2][k]
                exponent -= minimal[d] * elementary_finite
        assert exponent.denominator == 1
        exponents.append(int(exponent))
    return exponents


def raw_collision_counts(degree, genus, mixed):
    """Enumerate product-of-commutator distributions independently of characters."""
    group = tuple(permutations(range(degree)))
    identity = group[0]
    inverse = {a: tuple(a.index(i) for i in identity) for a in group}
    ordinary, first = Counter(), Counter()
    for a in group:
        for b in group:
            comm = compose(compose(a, b), compose(inverse[a], inverse[b]))
            ordinary[comm] += 1
            first[comm] += permutation_sign(a) if mixed else 1
    distribution = first
    for _ in range(1, genus):
        updated = Counter()
        for a, count in distribution.items():
            for b, other in ordinary.items():
                updated[compose(a, b)] += count * other
        distribution = updated
    # A permutation with exactly three moved sheets is a 3-cycle.
    triple = sum(count for a, count in distribution.items() if sum(i != j for i, j in enumerate(a)) == 3)
    return (Fraction(distribution[identity], factorial(degree)),
            Fraction(triple, factorial(degree)))


if __name__ == "__main__":
    for genus in (2, 3):
        for mixed in (False, True):
            mass, raw_triples = [1], [0]
            for m in range(1, 6):
                b, t = raw_collision_counts(m, genus, mixed)
                mass.append(b)
                raw_triples.append(t)
            nodes, triples, _, _ = spin_series(genus, 5, mixed)
            assert triples == divide_series(raw_triples, mass)
            assert nodes == divide_series([comb(m, 2)*b for m, b in enumerate(mass)], mass)
            if not mixed:
                assert first_exponents(genus, 8, False, 1) == [0]*9
                continue
            for epsilon in (1, -1):
                exponents = first_exponents(genus, 8, mixed, epsilon)
                assert exponents[2] == epsilon * (genus-1) * 2**(2*genus-3)
                assert exponents[3] == epsilon * (genus-1) * (11*6**(2*genus-2) - 3*2**(2*genus-2)) // 6
                print(f"g={genus}, epsilon={epsilon:+}, first-excess exponents 1..8 = {exponents[1:]}")
