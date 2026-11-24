"""Exact algebra for the single and simultaneous smoothing formulas (Notes 120–122).

Run with uv run --with sympy python checks/ramified_cubic_smoothing_degrees.py.
The analytic estimates and the determinant orientations on the smooth covers are
inputs to the formulas, and are established separately in the research notes.
"""

from itertools import combinations, product

from sympy import I, Matrix, Rational, conjugate, diag, simplify, zeros


def check_sheet_derivative(root, velocity):
    """Differentiate the conjugation interpolant at three distinct central roots."""
    w = root**2
    roots = (0, root, -root)
    increments = [velocity / (w - 3 * z**2) for z in roots]
    assert all(simplify(a - b) == 0 for a, b in zip(
        increments, [velocity / w, -velocity / (2 * w), -velocity / (2 * w)]))
    old_slope = conjugate(root) / root
    rhs = Matrix([conjugate(dz) - old_slope * dz for dz in increments])
    interpolation = Matrix([[1, z, z**2] for z in roots])
    constant, linear, quadratic = interpolation.inv() * rhs
    residue = conjugate(velocity) / conjugate(w) - old_slope * velocity / w
    assert simplify(linear) == 0
    assert simplify(quadratic + 3 * residue / (2 * w)) == 0
    assert simplify(constant + w * quadratic + residue / 2) == 0


def antilinear_matrix(a, b):
    """Real matrix of q ↦ (a+ib) conjugate(q)."""
    return Matrix([[a, b], [b, -a]])


def check_subset(blocks, evaluations, rows):
    matching = diag(*blocks)
    diagonal_height = Matrix.vstack(*evaluations)
    adjoint = Matrix.hstack(*rows)
    weighted_sum = (adjoint * matching.inv() * diagonal_height)[0]
    individual_sum = sum((h * m.inv() * e)[0] for m, e, h in zip(blocks, evaluations, rows))
    assert individual_sum == weighted_sum
    # Restore the common height as one real column before forming the matching quotient.
    full = matching.row_join(diagonal_height).col_join(adjoint.row_join(zeros(1, 1)))
    assert full.det() == -matching.det() * weighted_sum
    return full.det()


def main():
    for root, velocity in product([1 + I, 2 - I, 1 + 3 * I], [1, I, 2 - I]):
        check_sheet_derivative(root, velocity)
    print("Nine exact interpolations verified both the singular coefficient and its half-trace.")

    blocks = [antilinear_matrix(2, 1), antilinear_matrix(-1, 2),
              antilinear_matrix(1, -3), antilinear_matrix(3, 1)]
    evaluations = [Matrix([1, 2]), Matrix([1, -1]), Matrix([3, 1]), Matrix([0, 2])]
    rows = [Matrix([[1, 3]]), Matrix([[3, -1]]), Matrix([[2, 1]]), Matrix([[1, 2]])]
    for size in range(1, 5):
        for subset in combinations(range(4), size):
            check_subset([blocks[i] for i in subset], [evaluations[i] for i in subset],
                         [rows[i] for i in subset])
    print("Verified the signed determinant identity for all 15 nonempty subsets of four flags.")

    # In fixed orientations positive coefficient rescaling at the second flag can pass
    # through the two-node degeneracy while both one-node determinants retain sign.
    lists = []
    for positive_scale in (2, Rational(1, 2)):
        blocks = [antilinear_matrix(1, 0), positive_scale * antilinear_matrix(1, 0)]
        ev = [Matrix([1, 0]), Matrix([1, 0])]
        tests = [Matrix([[1, 0]]), Matrix([[-1, 0]])]
        singles = [check_subset([blocks[i]], [ev[i]], [tests[i]]) for i in range(2)]
        lists.append([*singles, check_subset(blocks, ev, tests)])
    assert lists[0][0] * lists[1][0] > 0 and lists[0][1] * lists[1][1] > 0
    assert lists[0][2] * lists[1][2] < 0
    assert {a + b for a, b in product([-1, 1], repeat=2)} == {-2, 0, 2}
    print("Both one-node signs agree between examples and the two-node sign reverses.")


if __name__ == "__main__":
    main()
