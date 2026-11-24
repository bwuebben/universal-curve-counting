#!/usr/bin/env python3
"""Exact branch-polynomial and cubic-calibration checks for Note 064.

Uses the standard library and the project's existing exact hook-count code.
It does not evaluate the unknown main mixed-spin integral.
"""

from fractions import Fraction as Q
from etale_elementary_exponents import connected_counts
from first_excess_spin_hurwitz import spin_series
from prime_degree_genus_cutoffs import add, disconnected, multiply, power, scale


def determinant(matrix):
    matrix = [list(map(Q, row)) for row in matrix]
    result = Q(1)
    for column in range(len(matrix)):
        pivot = next((i for i in range(column, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            return Q(0)
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            result = -result
        value = matrix[column][column]
        result *= value
        for row in range(column + 1, len(matrix)):
            ratio = matrix[row][column] / value
            for k in range(column + 1, len(matrix)):
                matrix[row][k] -= ratio * matrix[column][k]
    return result


def discriminant(coefficients):
    """Polynomial coefficients are listed in ascending order."""
    n = len(coefficients) - 1
    derivative = [i * coefficients[i] for i in range(1, n + 1)]
    f, df = coefficients[::-1], derivative[::-1]
    matrix = []
    for offset in range(n - 1):
        matrix.append([0] * offset + f + [0] * (n - 2 - offset))
    for offset in range(n):
        matrix.append([0] * offset + df + [0] * (n - 1 - offset))
    return Q((-1) ** (n * (n - 1) // 2), coefficients[-1]) * determinant(matrix)


def evaluate(coefficients, x):
    result = Q(0)
    for coefficient in reversed(coefficients):
        result = result * x + coefficient
    return result


def sine_coefficient(exponent):
    # Direct multinomial expansion of (sin(t/2)/(t/2))^exponent at t^4.
    return exponent * Q(1, 1920) + Q(exponent * (exponent - 1), 2 * 24**2)


def main():
    examples = {
        "(1,1)": (lambda w: [1, -w, 1], [-4, 0, 1]),
        "(2,1)": (lambda w: [1, -w, 0, 1], [-27, 0, 0, 4]),
        "(3,1)": (lambda w: [1, -w, 0, 0, 1], [256, 0, 0, 0, -27]),
        "(2,2)": (lambda w: [1, 0, -w, 1, 1], [229, -144, -128, 4, 16]),
        "tripod": (lambda w: [-1, w + 2, -w - 1, 1], [-23, 6, -5, -2, 1]),
    }
    for label, (fiber, branch) in examples.items():
        degree = len(fiber(0)) - 1
        # The resultant gives a polynomial in w of degree at most 2n-2.
        # More than that many exact evaluations prove the displayed identity.
        for value in range(2 * degree):
            assert discriminant(fiber(value)) == evaluate(branch, value)
        branch_discriminant = discriminant(branch)
        assert branch_discriminant
        print(f"Bubble {label}: branch discriminant {branch_discriminant}")

    for genus in range(2, 10):
        for degree in range(2, 20):
            k = degree * (genus - 1)
            assert sine_coefficient(2 * k) == Q(k * (5 * k - 1), 1440)

    counts = {mixed: connected_counts(2, 3, mixed)[3] for mixed in (False, True)}
    finite = {mixed: spin_series(2, 3, mixed)[2][3] for mixed in (False, True)}
    assert counts == {False: Q(220, 3), True: Q(28, 3)}
    assert finite == {False: -70, True: -6}
    tail = {mixed: -Q(1, 4) * finite[mixed] for mixed in finite}
    etale = {mixed: sine_coefficient(6) * counts[mixed] for mixed in counts}
    assert tail == {False: Q(35, 2), True: Q(3, 2)}
    assert etale == {False: Q(77, 36), True: Q(49, 180)}

    z1, z2, z3 = (disconnected(degree) for degree in (1, 2, 3))
    elementary = add(z3, scale(multiply(z1, z2), -1), scale(power(z1, 3), Q(1, 3)))
    equal_total = elementary[10]
    assert equal_total == Q(1673, 36)
    equal_main = equal_total - tail[False] - etale[False]
    assert equal_main == Q(161, 6)
    mixed_explicit = tail[True] + etale[True]
    assert mixed_explicit == Q(319, 180)
    exponent_offset = Q(4189, 180) - mixed_explicit
    assert exponent_offset == Q(43, 2)
    assert equal_main - exponent_offset == Q(16, 3)
    tail_defect_pairing = -Q(3, 48) * (2 * finite[False] - 2 * finite[True])
    assert tail_defect_pairing == 8
    assert Q(16, 3) + 2 * tail_defect_pairing == Q(64, 3)
    print("Cubic main calibration: I_equal=161/6; n_3,6=I_mixed-43/2.")
    print("Known elliptic-tail defect pairing: 8. Unknown main integral remains unevaluated.")


if __name__ == "__main__":
    main()
