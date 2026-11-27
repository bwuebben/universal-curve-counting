"""Verify two-slice commutation in homogeneous degree three over the rationals.

Run sage -python checks/cubic_slices.py --all-characters --workers 4
from the manuscript directory for the full interpolation calculation.
The default runs one character pair and reports the seven row ranks.
"""

import argparse
from collections import Counter
from functools import cache
from itertools import combinations, permutations, product
from math import factorial, prod
from multiprocessing import get_context

from sage.all import QQ, matrix, prime_range


def contraction(q, V, t):
    return (1 - t) * (1 - q * t / V) / ((1 - q * t) * (1 - t / V))


def allocations(labels, parts):
    if not parts:
        yield ()
        return
    for group in combinations(labels, parts[0]):
        other = tuple(z for z in labels if z not in group)
        for tail in allocations(other, parts[1:]):
            yield (group,) + tail


class SliceWeights:
    def __init__(self, p, s):
        self.p, self.s, self.c = p, s, p / s
        self.ell = (1 - self.c) / (1 - p)

    @cache
    def between(self, left, right, positive):
        q, V = (self.p, self.s) if positive else (1 / self.p, 1 / self.s)
        return prod(contraction(q, V, b / a) for a in left for b in right)

    @cache
    def elementary(self, group, positive):
        V = self.s if positive else 1 / self.s
        return prod(
            V * (a - b) * (b - a) / ((V * a - b) * (V * b - a))
            for a, b in combinations(group, 2)
        )

    @cache
    def power_sum(self, group, positive):
        if len(group) == 1:
            return 1
        singles = sum(
            prod(self.between((a,), (b,), positive) for a, b in combinations(order, 2))
            for order in permutations(group)
        )
        if len(group) == 2:
            return singles - 2 * self.elementary(group, positive)
        assert len(group) == 3
        mixed = sum(
            self.elementary(pair, positive) * self.between(single, pair, positive)
            for single, pair in allocations(group, (1, 2))
        )
        return singles - 3 * mixed + 3 * self.elementary(group, positive)

    @cache
    def transport(self, labels, X):
        return prod((X - z) / (X - self.s * z) for z in labels)

    @cache
    def side(self, labels, parts, positive, X):
        total = 0
        for groups in allocations(labels, parts):
            internal = prod(self.power_sum(group, positive) for group in groups)
            for flags in product((0, 1), repeat=len(parts)):
                weight = internal
                ordered = sorted(
                    range(len(parts)),
                    key=lambda i: (flags[i] if positive else -flags[i], i),
                )
                for i, j in combinations(ordered, 2):
                    weight *= self.between(groups[i], groups[j], positive)
                for group, shift in zip(groups, flags):
                    if shift:
                        weight *= -self.transport(group, X) * (
                            1 if positive else self.s ** len(group)
                        )
                total += weight
        return total


PARTS = {0: ((),), 1: ((1,),), 2: ((2,), (1, 1)), 3: ((3,), (2, 1), (1, 1, 1))}


def two_slices(data, labels, X, Y, positive):
    row = []
    for r in range(4):
        splittings = [
            (tuple(first), tuple(j for j in labels if j not in first))
            for first in combinations(labels, r)
        ]
        for left_parts in PARTS[r]:
            for right_parts in PARTS[3 - r]:
                value = 0
                for first, second in splittings:
                    transport = (
                        data.transport(second, X)
                        if positive
                        else data.transport(first, Y)
                    )
                    value += (
                        transport
                        * data.between(first, second, positive)
                        * data.side(first, left_parts, positive, X)
                        * data.side(second, right_parts, positive, Y)
                    )
                row.append(value)
    return row


def pulse_scalars(data):
    def scalar(parts):
        return prod(
            ((1 - data.c**a) / (a * (1 - data.p**a) * data.ell**a)) ** m
            / QQ(factorial(m))
            for a, m in Counter(parts).items()
        )

    return [
        scalar(left) * scalar(right)
        for r in range(4)
        for left in PARTS[r]
        for right in PARTS[3 - r]
    ]


def certificate(p, s, *, verbose=False):
    data = SliceWeights(p, s)
    weights = pulse_scalars(data)
    # Degrees in each label are at most ten after the common denominator is cleared.
    # Thirteen nodes leave eleven choices distinct from the other two labels.
    nodes = tuple(map(QQ, prime_range(11, 60)))
    assert len(nodes) == 13
    Y = QQ(-1)
    for xi in range(2, 9):
        X = QQ(-xi)
        positive_rows = []
        negative_rows = []
        for labels in combinations(nodes, 3):
            P = two_slices(data, labels, X, Y, True)
            P_swapped = two_slices(data, labels, Y, X, True)
            N = two_slices(data, labels, X, Y, False)
            N_swapped = two_slices(data, labels, Y, X, False)
            positive_rows.append(P + [-v for v in P_swapped])
            negative_rows.append(
                [w * v for w, v in zip(weights + weights, N + N_swapped)]
            )
        positive = matrix(QQ, positive_rows)
        negative = matrix(QQ, negative_rows)
        basis = positive.row_space().basis_matrix()
        assert (basis * negative.transpose()).is_zero()
        if verbose:
            print(
                "X/Y =",
                xi,
                "rank =",
                basis.nrows(),
                "all six-label assignments certified.",
                flush=True,
            )
        data.side.cache_clear()
        data.transport.cache_clear()
    # These method caches include the instance. Release the previous character pair
    # before this worker begins another specialization.
    data.between.cache_clear()
    data.power_sum.cache_clear()
    data.elementary.cache_clear()


def one_specialization(pair):
    p, s = pair
    certificate(QQ(p), QQ(s))
    return pair


def all_characters(workers):
    inputs = [(p, s) for p in range(2, 20) for s in range(21, 51)]
    verified = set()
    with get_context("spawn").Pool(workers) as pool:
        for pair in pool.imap_unordered(one_specialization, inputs):
            verified.add(pair)
            if len(verified) % 30 == 0:
                print(
                    len(verified),
                    "of",
                    len(inputs),
                    "character pairs verified.",
                    flush=True,
                )
    assert verified == set(inputs)
    print(
        "Passed over Q(p,s,X,Y,A1,A2,A3,B1,B2,B3) by the degree bounds in the manuscript.",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--all-characters", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.all_characters:
        all_characters(args.workers)
    else:
        certificate(QQ(3), QQ(4), verbose=True)
