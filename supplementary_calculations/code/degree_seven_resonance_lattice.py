"""Exact subset signs, integral lattices, and elliptic-tail coefficients in Note 125."""

from fractions import Fraction
from itertools import product
from math import comb, factorial

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def subset_sums(weights):
    sums = [0]
    for value in weights:
        sums.extend([term + value for term in sums])
    return sums


def jump_vector(sums, size):
    jump = [0] * size
    for mask, value in enumerate(sums[1:], start=1):
        assert value
        degree = mask.bit_count()
        jump[degree - 1] += (-1) ** degree * (1 if value > 0 else -1)
    return jump


def crossing_witness(size, mask):
    """Choose one point on the selected hyperplane and no other, by binary uniqueness."""
    distinguished = (mask & -mask).bit_length() - 1
    weights = [2 << i for i in range(size)]
    weights[distinguished] = 0
    weights[distinguished] = -sum(w for i, w in enumerate(weights) if mask & (1 << i))
    sides = []
    for offset in (-1, 1):
        side = weights.copy()
        side[distinguished] += offset
        sides.append(subset_sums(side))
    return sides


def sine_power(power, coefficient):
    """Coefficient of v^(2*coefficient) in [2(1-cos v)]^power."""
    x = [Fraction(0)] + [Fraction(2 * (-1) ** (i + 1), factorial(2 * i))
                         for i in range(1, coefficient + 1)]
    result = [Fraction(1)] + [Fraction(0)] * coefficient
    for _ in range(power):
        result = [sum(result[k - i] * x[i] for i in range(1, k + 1))
                  for k in range(coefficient + 1)]
    return result[coefficient]


def main():
    tail = sine_power(11, 12)
    assert tail == Fraction(-11, 12) and sine_power(12, 12) == 1
    pairs = {tuple(jump_vector(subset_sums(w), 2))
             for w in product(range(-3, 4), repeat=2) if w[0] * w[1] * sum(w)}
    assert pairs == {(-2, 1), (0, -1), (0, 1), (2, -1)}
    assert {a * tail + b for a, b in pairs} == {Fraction(17, 6), -1, 1, Fraction(-17, 6)}
    print("Genus twelve/thirteen primitive pairs:", sorted(pairs))
    print("Connected genus-thirteen changes:", sorted(str(a * tail + b) for a, b in pairs))

    for size in (2, 4, 6, 8, 10):
        # Generate the one positive-weight row and actual chamber differences.
        all_positive = jump_vector(subset_sums([1] * size), size)
        generators = [all_positive]
        masks = range(1, 1 << size) if size <= 8 else [(1 << j) - 1 for j in range(1, size + 1)]
        for mask in masks:
            left, right = crossing_witness(size, mask)
            flips = [k for k in range(1, len(left)) if left[k] * right[k] < 0]
            assert flips == [mask]
            before, after = jump_vector(left, size), jump_vector(right, size)
            difference = [b - a for a, b in zip(before, after)]
            j = mask.bit_count()
            assert difference == [2 * (-1) ** j if i == j else 0 for i in range(1, size + 1)]
            assert all((value - comb(size, i)) % 2 == 0 for i, value in enumerate(after, 1))
            if mask == (1 << j) - 1:
                generators.append(difference)
        form = smith_normal_form(Matrix(generators), domain=ZZ)
        diagonal = [abs(form[i, i]) for i in range(size)]
        assert diagonal == [1] + [2] * (size - 1)
        print(f"b={size}: {len(masks)} hyperplanes, lattice index {2 ** (size - 1)}.")

    # Binary-encoded polynomials: in characteristic two (U^l+V^l)/T^(3l)
    # equals (1+T^2)^l+1. Expand via powers selected by binary digits of l.
    for l in range(1, 65):
        p = 1
        for bit in range(l.bit_length()):
            if l & (1 << bit):
                p ^= p << (2 << bit)
        p = (p ^ 1) << (3 * l)
        by_subsets = sum((comb(2 * l, j) % 2) << (3 * l + j) for j in range(1, 2 * l + 1))
        assert p == by_subsets and p.bit_length() - 1 == 5 * l
    print("Cross-branch parity polynomials verified through 128 branch values.")


if __name__ == "__main__":
    main()
