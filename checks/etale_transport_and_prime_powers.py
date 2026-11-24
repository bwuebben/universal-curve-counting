"""Exact finite-cover tower and least-genus arithmetic checks for the transport and prime-power formulas.

The tower uses the action of D4 on the square. The series attached to a
stabilizer are test inputs for the formal transfer identity. The arithmetic
tests use the ordinary and sign-twisted surface-group hook series.
"""

from fractions import Fraction
from itertools import product
from math import gcd

from dihedral_wall_transfer import subgroups
from etale_spin_characters import hook_count


def multiply(left, right):
    """Product of series available through the same cutoff."""
    return [sum(left[i] * right[n - i] for i in range(n + 1))
            for n in range(len(left))]


def inverse(series):
    assert series[0] == 1
    result = [1]
    for n in range(1, len(series)):
        result.append(-sum(series[i] * result[n - i] for i in range(1, n + 1)))
    return result


def substitute(series, index, cutoff):
    assert len(series) == cutoff // index + 1
    return [series[d // index] if d % index == 0 else 0 for d in range(cutoff + 1)]


def translate(g, point):
    """Points are tuples of subsets of the square's four vertices."""
    offset, direction = g
    return tuple(frozenset((offset + direction * v) % 4 for v in subset)
                 for subset in point)


def orbit_towers(group, points):
    remaining = set(points)
    while remaining:
        point = next(iter(remaining))
        orbit = {translate(g, point) for g in group}
        stabilizer = frozenset(g for g in group if translate(g, point) == point)
        assert len(orbit) * len(stabilizer) == len(group)
        remaining.difference_update(orbit)
        yield len(orbit), stabilizer


def transport(virtual_set, family, group, cutoff):
    """Shift by each H-orbit's degree and evaluate its stabilizer's family."""
    result = [1] + [0] * cutoff
    for sign, points in virtual_set:
        for index, stabilizer in orbit_towers(group, points):
            if index > cutoff:
                continue  # this orbit contributes only its constant 1
            value = family(stabilizer, cutoff // index)
            if sign < 0:
                value = inverse(value)
            result = multiply(result, substitute(value, index, cutoff))
    return result


def sample_series(group, cutoff):
    # Fixed vertices distinguish the two classes of order-two reflections.
    vertices_fixed = sum(all((a + s * v) % 4 == v for a, s in group) for v in range(4))
    return [1] + [1 + (3 + len(group)) * d + vertices_fixed * d**2
                  for d in range(1, cutoff + 1)]


def check_towers():
    vertices = [(frozenset([i]),) for i in range(4)]
    edges = [(frozenset([i, (i + 1) % 4]),) for i in range(4)]
    diagonals = [(frozenset([i, i + 2]),) for i in range(2)]
    point = [(frozenset(range(4)),)]
    b, c = [(1, vertices), (-1, point)], [(1, edges), (-1, diagonals)]
    bc = [(sign_b * sign_c, [x + y for x, y in product(xs, ys)])
          for sign_b, xs in b for sign_c, ys in c]
    full = frozenset(product(range(4), (1, -1)))
    # At total cutoff M only covers of degree [G:K] <= M are available.
    # A second restriction has the smaller budget floor(M/[G:K_x]).
    flags = [v + edge for v, edge in product(vertices, edges) if v[0] <= edge[0]]
    assert len(flags) == 8  # a free transitive action
    cutoffs = (1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17)
    for total_cutoff in cutoffs:
        def visible_series(stabilizer, bound):
            degree = len(full) // len(stabilizer)
            assert degree <= total_cutoff and bound == total_cutoff // degree
            return sample_series(stabilizer, bound)

        def through_c(stabilizer, bound):
            return transport(c, visible_series, stabilizer, bound)

        for group in subgroups(4):
            degree = len(full) // len(group)
            if degree > total_cutoff:
                continue
            cutoff = total_cutoff // degree
            direct = transport(bc, visible_series, group, cutoff)
            assert direct == transport(b, through_c, group, cutoff), (group, total_cutoff)
            free = transport([(1, flags)], visible_series, group, cutoff)
            if total_cutoff < 8:
                assert free == [1] + [0] * cutoff

    # Independently decompose the simultaneous square action. Vertex/edge
    # pairs give two free orbits; vertex/diagonal pairs give two vertex orbits.
    trivial = frozenset([(0, 1)])
    vertex_stab = frozenset([(0, 1), (0, -1)])
    edge_stab = frozenset([(0, 1), (1, -1)])
    diag_stab = frozenset((a, s) for a in (0, 2) for s in (1, -1))
    assert sample_series(vertex_stab, 1) != sample_series(edge_stab, 1)
    for cutoff in (7, 8, 16, 17):
        def shifted(group, degree):
            return substitute(sample_series(group, cutoff // degree), degree, cutoff)

        free_factor = shifted(trivial, 8)
        vertex_factor = shifted(vertex_stab, 4)
        numerator = multiply(shifted(diag_stab, 2), multiply(free_factor, free_factor))
        denominator = multiply(shifted(edge_stab, 4), multiply(vertex_factor, vertex_factor))
        assert transport(bc, sample_series, full, cutoff) == multiply(numerator, inverse(denominator))
    print("D4 virtual towers agree on the visible subgroups for 12 cutoffs up to degree 17.")


def ghosts(series):
    """Integer coefficients of q U'(q) / U(q)."""
    assert series[0] == 1
    result = [0]
    for n in range(1, len(series)):
        result.append(n * series[n] - sum(result[j] * series[n - j] for j in range(1, n)))
    return result


def check_prime_powers():
    bound = 16
    primes = (2, 3, 5, 7, 11, 13)
    # r=b=2 gives an independent square. The next two rows have two
    # components in the product of the character covers.
    # At upstairs degree two every nonzero binary character gives one double
    # cover of deck weight 1/2, so these unsigned mixed counts also have a
    # direct count independent of the hook series.
    degree_two_counts = {
        (2, 2, 2): Fraction(1023, 2),
        (2, 3, 2): 1023,
        (2, 4, 3): 262143,
        (3, 2, 2): Fraction(262143, 2),
    }
    upstairs_genera = {2} | {1 + 2**b * (g - 1) for g, _, b in degree_two_counts}
    for genus in sorted(upstairs_genera):
        for twisted, sign in product((False, True), (1, -1)):
            series = [1] + [sign**n * hook_count(n, genus, twisted) for n in range(1, bound + 1)]
            g_values = ghosts(series)
            if genus == 2 and not twisted and sign == 1:
                assert g_values[1:5] == [1, 15, 220, 5275]
            if genus == 2 and twisted and sign == 1:
                assert g_values[1:4] == [1, -1, 28]
            if genus == 5 and not twisted and sign == 1:
                assert g_values[1:5] == [1, 1023, 10095844, 882442672507]
            exponents = [0]
            for n in range(1, bound + 1):
                remainder = g_values[n] - sum(d * exponents[d] for d in range(1, n) if n % d == 0)
                exponent, residue = divmod(remainder, n)
                assert residue == 0, (genus, twisted, sign, n)
                exponents.append(exponent)

            for p in primes:
                power = p
                while power <= bound:
                    for n in range(1, bound // power + 1):
                        if n % p:
                            difference = g_values[n * power] - g_values[n * (power // p)]
                            assert difference % power == 0
                    assert g_values[power] % p == sign % p
                    for (base_genus, r, b), expected_double in degree_two_counts.items():
                        joint_degree = 2**b
                        if 1 + joint_degree * (base_genus - 1) != genus:
                            continue
                        multiplier = 2**(r - b)
                        mixed_count = multiplier * Fraction(g_values[power], power)
                        if power == 2 and not twisted and sign == 1:
                            assert mixed_count == expected_double
                        assert mixed_count.denominator == power // gcd(power, multiplier)
                        assert 1 + joint_degree * power * (base_genus - 1) == 1 + power * (genus - 1)
                    power *= p
    print("Ordinary and sign-twisted hook series at upstairs genera 2, 5, 9 pass all prime-power")
    print("ghost congruences through 16, including the two dependent-character denominator reductions.")


if __name__ == "__main__":
    check_towers()
    check_prime_powers()
