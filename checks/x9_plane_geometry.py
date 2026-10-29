"""Exact fan calculations used in Paper 1; Python standard library only.

This certifies the displayed lattice data and nef inequalities, not the
geometric or enumerative conclusions drawn from them. Run from the project:
    python3 checks/x9_plane_geometry.py
Raw input: the seven-ray smoothing fan printed in cy_smoothing Paper 6 and
global_higgs_branch/checks/check_x9_neutral_quantum_periods.sage, pinned in
evidence/x9-inputs-2026-10-28.json.
"""

from fractions import Fraction as F
from itertools import combinations, product
import json

RAYS = [(-1, -3, -1, 0), (-1, -3, 0, -1), (-1, -1, 0, 0),
        (-1, 0, 0, 0), (0, 1, 0, 0), (1, -1, 0, 0), (1, 3, 1, 1)]
CHARGES = [(1, 1, -2, 0, 0, -1, 1), (0, 0, -1, 2, 0, 1, 0),
           (0, 0, 1, -1, 1, 0, 0)]
FIXED = [(1, 4, 5, 6), (1, 3, 4, 6), (0, 4, 5, 6),
         (0, 3, 4, 6), (0, 1, 4, 5), (0, 1, 3, 4),
         (1, 2, 3, 6), (0, 2, 3, 6), (0, 1, 2, 3)]
SINGULAR = (0, 1, 2, 5, 6)
PLUS = FIXED + [(1, 2, 5, 6), (0, 2, 5, 6), (0, 1, 2, 5)]
MINUS = FIXED + [(0, 1, 5, 6), (0, 1, 2, 6)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def det(a):
    if len(a) == 1:
        return a[0][0]
    return sum((-1)**j * a[0][j] * det([row[:j] + row[j+1:] for row in a[1:]])
               for j in range(len(a)))


def solve(a, b):
    denominator = det(a)
    if denominator == 0:
        return None
    return tuple(F(det([row[:j] + [b[i]] + row[j+1:]
                        for i, row in enumerate(a)]), denominator)
                 for j in range(len(a)))


def nef_certificate(cones, divisor):
    result = []
    for cone in cones:
        m = solve([list(RAYS[i]) for i in cone], [-int(i == divisor) for i in cone])
        assert m is not None
        slack = [dot(m, r) + int(i == divisor) for i, r in enumerate(RAYS)]
        assert min(slack) >= 0, (cone, divisor, m, slack)
        result.append({"cone": cone, "linear_form": list(map(str, m)),
                       "slack": list(map(str, slack))})
    return result


def main():
    assert all(sum(row[i] * RAYS[i][j] for i in range(7)) == 0
               for row in CHARGES for j in range(4))
    assert abs(det([[row[i] for i in (0, 4, 5)] for row in CHARGES])) == 1
    assert tuple(sum(c * RAYS[i][j] for i, c in zip((0, 1, 2, 3, 6), (1, 1, -3, 2, 1)))
                 for j in range(4)) == (0, 0, 0, 0)
    local_cones = [c for c in FIXED if 4 not in c and 5 not in c]
    assert local_cones == [(1, 2, 3, 6), (0, 2, 3, 6), (0, 1, 2, 3)]
    local_determinants = [det([list(RAYS[i]) for i in c]) for c in local_cones]
    assert local_determinants == [-1, 1, -1]

    facets = set()
    for triple in combinations(SINGULAR, 3):
        a = [list(RAYS[i]) for i in triple]
        normal = [(-1)**j * det([row[:j] + row[j+1:] for row in a]) for j in range(4)]
        if not any(normal):
            continue
        values = [dot(normal, RAYS[i]) for i in SINGULAR]
        if min(values) >= 0 or max(values) <= 0:
            facets.add(tuple(i for i, value in zip(SINGULAR, values) if value == 0))
    assert facets == {(0, 1, 2), (0, 1, 5), (0, 2, 6),
                      (0, 5, 6), (1, 2, 6), (1, 5, 6)}

    vertices = set()
    for subset in combinations(range(7), 4):
        m = solve([list(RAYS[i]) for i in subset], [-1] * 4)
        if m is not None and all(dot(m, r) >= -1 for r in RAYS):
            vertices.add(m)
    assert vertices
    bounds = [(int(min(v[j] for v in vertices)), int(max(v[j] for v in vertices)))
              for j in range(4)]
    assert all(value.denominator == 1 for v in vertices for value in v)
    # The complete toric fan makes this anticanonical polytope bounded.
    face = [m for m in product(*(range(lo, hi + 1) for lo, hi in bounds))
            if dot(m, RAYS[2]) == -1 and all(dot(m, r) >= -1 for r in RAYS)]
    exponents = [tuple(1 + dot(m, RAYS[i]) for i in (0, 1, 2, 3, 6)) for m in face]
    expected = {(0, 0, 0, 1, 0), (2, 0, 0, 0, 0), (1, 1, 0, 0, 0),
                (1, 0, 0, 0, 1), (0, 2, 0, 0, 0), (0, 1, 0, 0, 1), (0, 0, 0, 0, 2)}
    assert set(exponents) == expected and len(face) == 7
    vertex = (0, 1, -2, -2)
    assert [i for i, r in enumerate(RAYS) if dot(vertex, r) == -1] == list(SINGULAR)
    assert exponents[face.index(vertex)] == (0, 0, 0, 1, 0)
    result = {"rays": RAYS, "charges": CHARGES, "local_cones": local_cones,
              "local_determinants": local_determinants, "singular_facets": sorted(facets),
              "anticanonical_face": [{"lattice_point": m, "exponents_01236": e}
                                      for m, e in zip(face, exponents)],
              "nef": {"plus_E0": nef_certificate(PLUS, 0),
                      "plus_E4": nef_certificate(PLUS, 4),
                      "minus_E5": nef_certificate(MINUS, 5),
                      "minus_E4": nef_certificate(MINUS, 4)}}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
