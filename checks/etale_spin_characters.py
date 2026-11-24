"""Verify the étale spin character formula by exact commutator counts in S_m."""

from collections import Counter
from itertools import permutations
from math import factorial, prod


def partitions(total, upper=None):
    if total == 0:
        yield ()
        return
    for first in range(min(total, upper or total), 0, -1):
        for remaining in partitions(total - first, first):
            yield (first,) + remaining


def hook_count(degree, genus, twisted):
    total = 0
    for shape in partitions(degree):
        transpose = tuple(sum(row > j for row in shape) for j in range(shape[0]))
        if twisted and shape != transpose:
            continue
        hook_product = prod(
            row - j + transpose[j] - i - 1
            for i, row in enumerate(shape)
            for j in range(row)
        )
        total += hook_product ** (2 * genus - 2)
    return total


def compose(left, right):
    return tuple(left[i] for i in right)


def permutation_sign(element):
    inversions = sum(a > b for i, a in enumerate(element) for b in element[i + 1:])
    return (-1) ** inversions


def commutator_counts(degree, genus):
    elements = list(permutations(range(degree)))
    identity = elements[0]
    inverse = {a: tuple(a.index(i) for i in identity) for a in elements}
    ordinary, twisted = Counter(), Counter()
    for a in elements:
        for b in elements:
            commutator = compose(compose(a, b), compose(inverse[a], inverse[b]))
            ordinary[commutator] += 1
            twisted[commutator] += permutation_sign(a)
    counts = []
    for first in (ordinary, twisted):
        distribution = first
        for _ in range(1, genus):
            updated = Counter()
            for a, left in distribution.items():
                for b, right in ordinary.items():
                    updated[compose(a, b)] += left * right
            distribution = updated
        quotient, remainder = divmod(distribution[identity], factorial(degree))
        assert remainder == 0
        counts.append(quotient)
    return counts


if __name__ == "__main__":
    for genus in (2, 3):
        for degree in (1, 2, 3, 4):
            actual = commutator_counts(degree, genus)
            assert actual == [hook_count(degree, genus, twist) for twist in (False, True)]
            print(f"g={genus}, degree={degree}: ordinary {actual[0]}, twisted {actual[1]}")
