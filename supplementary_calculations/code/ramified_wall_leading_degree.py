"""Exact checks for the first genus at a ramified double wall in each degree.

Compute pointed connected unramified-cover counts from the partition formula
and independently from transitive permutation actions in small degree. Check
the parity of the pointed count and the minimum-genus alternatives. These are
combinatorial tests. The individual Fredholm determinant signs are input data
of the wall.
"""

from collections import defaultdict
from itertools import permutations
from math import factorial, isqrt

from etale_spin_characters import compose, hook_count


def pointed_counts(disconnected):
    """Coefficients of q A'(q)/A(q), using A(0)=1 and exact division of series."""
    assert disconnected[0] == 1
    pointed = [0]
    for degree in range(1, len(disconnected)):
        count = degree * disconnected[degree]
        count -= sum(pointed[j] * disconnected[degree - j] for j in range(1, degree))
        pointed.append(count)
    return pointed


def direct_pointed_count(degree, genus):
    """Enumerate transitive actions with product of the handle commutators one."""
    group = list(permutations(range(degree)))
    identity = tuple(range(degree))
    inverse = {g: tuple(g.index(i) for i in identity) for g in group}
    pairs = defaultdict(list)
    for a in group:
        for b in group:
            commutator = compose(compose(a, b), compose(inverse[a], inverse[b]))
            pairs[commutator].append((a, b))

    def transitive(generators):
        reached = {0}
        todo = [0]
        while todo:
            i = todo.pop()
            for generator in generators:
                j = generator[i]
                if j not in reached:
                    reached.add(j)
                    todo.append(j)
        return len(reached) == degree

    def enumerate_handles(left, product, generators):
        if left == 1:
            return sum(transitive(generators + pair) for pair in pairs.get(inverse[product], ()))
        count = 0
        for commutator, choices in pairs.items():
            next_product = compose(product, commutator)
            for pair in choices:
                count += enumerate_handles(left - 1, next_product, generators + pair)
        return count

    labeled = enumerate_handles(genus, identity, ())
    # Relabelings fixing sheet 0 act freely on every transitive action.
    pointed, remainder = divmod(labeled, factorial(degree - 1))
    assert remainder == 0
    return pointed


def odd_part_is_square(n):
    while n % 2 == 0:
        n //= 2
    return isqrt(n) ** 2 == n


def staircase_parities(limit):
    """Formal division in F_2 using triangular exponents for A(q)."""
    triangles = []
    step = 1
    while (triangle := step * (step + 1) // 2) <= limit:
        triangles.append(triangle)
        step += 1
    exponents = set(triangles)
    parity = [0]
    for degree in range(1, limit + 1):
        rhs = degree % 2 if degree in exponents else 0
        for triangle in triangles:
            if triangle >= degree:
                break
            rhs ^= parity[degree - triangle]
        parity.append(rhs)
    return parity, exponents


def minimum_excess(degree, branch_pairs):
    """Minimize over all degree splits before adding ramification of the factors."""
    candidates = []
    for sheet_sum in range(1, degree // 2 + 1):
        zero_degree = degree - 2 * sheet_sum
        zero_counts = range(1, zero_degree + 1) if zero_degree else (0,)
        for active in range(1, sheet_sum + 1):
            for zero in zero_counts:
                excess = active + zero - 1 + sheet_sum * branch_pairs
                candidates.append((excess, (sheet_sum, active, zero)))
    minimum = min(excess for excess, _ in candidates)
    return minimum, {shape for excess, shape in candidates if excess == minimum}


def main():
    max_degree = 27
    parity_limit = 1000
    parity, triangles = staircase_parities(parity_limit)
    assert all(parity[n] == odd_part_is_square(n) for n in range(1, parity_limit + 1))
    for genus in (2, 3, 4):
        disconnected = [1] + [hook_count(n, genus, False) for n in range(1, max_degree + 1)]
        pointed = pointed_counts(disconnected)
        for n in range(1, max_degree + 1):
            assert disconnected[n] % 2 == (n in triangles)
            assert pointed[n] % 2 == parity[n]
        direct_limit = {2: 4, 3: 3, 4: 2}[genus]
        direct = [direct_pointed_count(n, genus) for n in range(1, direct_limit + 1)]
        assert direct == pointed[1:direct_limit + 1]
        print(f"g={genus} direct transitive counts through m={direct_limit}: {direct}")

    exceptional = set()
    divisor_terms = set()
    for branch_pairs in range(1, 11):
        minima = {}
        for degree in range(2, 33):
            threshold, shapes = minimum_excess(degree, branch_pairs)
            minima[degree] = threshold
            if degree == 2:
                assert threshold == branch_pairs and shapes == {(1, 1, 0)}
                continue
            assert threshold == branch_pairs + 1
            assert (1, 1, 1) in shapes
            for extra_shape in shapes - {(1, 1, 1)}:
                exceptional.add((degree, branch_pairs, extra_shape))
            for lower_degree in range(2, degree):
                if degree % lower_degree == 0:
                    multiple = degree // lower_degree
                    if multiple * minima[lower_degree] <= threshold:
                        divisor_terms.add((degree, branch_pairs, lower_degree))
    assert exceptional == {(4, 1, (2, 1, 0))}
    assert divisor_terms == {(4, 1, 2)}
    print(f"Hook recurrence verified through m={max_degree} and staircase parity through m={parity_limit}.")
    print("Genus and proper-divisor checks for 3 <= d <= 32, 1 <= b/2 <= 10: only (d,b/2)=(4,1) ties.")


if __name__ == "__main__":
    main()
