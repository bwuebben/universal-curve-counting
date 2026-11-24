"""Exact finite-dimensional checks for Note 072. Requires SymPy."""

from itertools import permutations
from math import prod

from sympy import I, Matrix, Rational, diag, eye, zeros


def residue_coordinates(contact, clutch):
    """Extract the leading poles and project them to the positive real form."""
    width = contact - 1
    # Residue of the product pairs orders j and contact-j on each branch.
    reversal = eye(width)[::-1, :]
    gram = diag(reversal, reversal)
    poles = zeros(2 * width, width)
    orders = []
    for j in range(1, contact):
        # Z^{-j} in the left frame is Z^{contact-j}/clutch in the right.
        exponents = (j, contact - j)
        order = min(exponents)
        orders.append(order)
        if exponents[0] == order:
            poles[j - 1, j - 1] = 1
        if exponents[1] == order:
            poles[width + contact - j - 1, j - 1] = 1 / clutch
    assert poles.T * gram * poles == zeros(width, width)

    # Fixed real vectors satisfy v_{k-j} = conjugate(v_j) on each branch.
    # Read real and imaginary coordinates for a pair and the real middle entry.
    columns = []
    for j in range(width):
        for scalar in (1, I):
            v = scalar * poles[:, j]
            real_projection = (v + gram * v.conjugate()) / 2
            coordinates = []
            for offset in (0, width):
                for p in range(1, contact):
                    entry = real_projection[offset + p - 1]
                    if 2 * p < contact:
                        coordinates.extend(entry.as_real_imag())
                    elif 2 * p == contact:
                        coordinates.append(entry.as_real_imag()[0])
            columns.append(Matrix(coordinates))
    determinant = Matrix.hstack(*columns).det() if columns else 1
    assert determinant != 0
    return poles, determinant, sum(orders)


def cycles(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            length += 1
            point = permutation[point]
        lengths.append(length)
    return lengths


def main():
    signs = {}
    for contact in range(1, 17):
        plus, det_plus, total_order = residue_coordinates(contact, I)
        minus, det_minus, _ = residue_coordinates(contact, -I)
        signs[contact] = int(det_minus / det_plus)
        codimension = plus.row_join(minus).rank() - (contact - 1)
        assert codimension == int(contact % 2 == 0)
        assert signs[contact] == (-1) ** codimension
        assert total_order == contact**2 // 4

    for degree in range(1, 8):
        for permutation in permutations(range(degree)):
            inversions = sum(left > right for i, left in enumerate(permutation)
                             for right in permutation[i + 1:])
            assert prod(signs[k] for k in cycles(permutation)) == (-1) ** inversions

    for degree in range(2, 18, 2):
        for genus in range(2, 11):
            denominator = 24 * degree
            main_degree = 2 * degree * (genus - 2)
            bridge_degree = 2 * degree
            assert Rational(bridge_degree, denominator) == Rational(1, 12)
            assert Rational(main_degree + bridge_degree, denominator) == Rational(genus - 1, 12)
    print("Checked residue planes and mode orders for contacts 1–16,")
    print("permutation signs through degree 7 and cyclic elliptic-tail degrees through 16.")


if __name__ == "__main__":
    main()
