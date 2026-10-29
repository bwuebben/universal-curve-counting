"""Exact divisor-coordinate check for Paper 1; not a flop-invariance proof."""

import json
from x9_plane_geometry import det

NEF = [[0, 0, 1, 1, 2, 2], [1, 1, 2, 4, 5, 5],
       [0, 1, 2, 4, 4, 5], [0, -1, -1, -3, -3, -3],
       [0, -1, -1, -2, -2, -2], [0, 0, 1, 2, 2, 2]]
CHARGES = [[0, 0, -1, 0, -1, 1, 0, 1, 0, 0],
           [0, 0, 0, 1, 1, 0, 0, 0, -1, -1],
           [0, 1, 0, 0, 0, 0, 0, 2, -2, 1],
           [1, 0, -1, 0, 0, 0, 0, 0, 0, 1],
           [0, 0, 1, 1, 0, 0, -1, -1, 0, 0],
           [0, 0, 0, -1, 0, 0, 1, 0, 1, -1]]
RAYS = [(1, 0, 0, 0), (0, 1, 0, 0), (-1, -1, 0, 0),
        (0, 0, 1, 0), (0, 0, 0, 1), (-4, -2, -1, 0),
        (-4, -2, 0, -1), (3, 1, 1, 1), (2, 1, 1, 1), (-2, -1, 0, 0)]


def main():
    assert det(NEF) == -1
    assert all(sum(row[i] * RAYS[i][j] for i in range(10)) == 0
               for row in CHARGES for j in range(4))
    basis_indices = (2, 5, 6, 7, 8, 9)
    recovered = [[sum(NEF[i][k] * CHARGES[k][j] for k in range(6))
                  for j in basis_indices] for i in range(6)]
    assert recovered == [[int(i == j) for j in range(6)] for i in range(6)]
    original = [(0, 0, 0, 0, 0, 1), (0, 0, 0, 0, -1, 0),
                (1, 0, 0, 0, 1, 1), (1, 0, 0, 0, 0, 0)]
    signed_rays = [({6, 8}, {3, 9}), ({6, 7}, {2, 3}),
                   ({5, 8}, {4, 9}), ({5, 7}, {2, 4})]
    for curve, (positive, negative) in zip(original, signed_rays):
        degrees = [sum(curve[k] * CHARGES[k][j] for k in range(6)) for j in range(10)]
        assert degrees == [int(j in positive) - int(j in negative) for j in range(10)]
    assert all(sum(s * c[j] for s, c in zip((1, -1, -1, 1), original)) == 0
               for j in range(6))
    kernel_minor = [[original[col][row] for col in (0, 1, 3)] for row in (0, 4, 5)]
    assert abs(det(kernel_minor)) == 1
    print(json.dumps({"nef_determinant": det(NEF), "original_exceptional_classes": original,
                      "flopped_node": "q2", "kernel_minor_determinant": det(kernel_minor),
                      "specialization_matrix": [[0, 1, 0, 0, 0, 0],
                                                [0, 0, 1, 0, 0, 0],
                                                [0, 0, 0, 1, 0, 0]]}, indent=2))


if __name__ == "__main__":
    main()
