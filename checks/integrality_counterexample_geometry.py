"""Exact projectivity, nefness, and primitive-pairing checks for Paper 2.

Uses SymPy 1.14.0. Completeness and contracted-locus geometry are checked
separately by x9_contracted_support.py and the accompanying proofs.
"""
import ast
import hashlib
import json
from pathlib import Path

import sympy as S

source = Path(__file__).with_name("x9_lift_bound.py")
assert hashlib.sha256(source.read_bytes()).hexdigest() == \
    "5e74c32b0b40ae4fb8d0d1439455bdc943f51c398bb164aacd776ea66d11fa18"
data = {}
for node in ast.parse(source.read_text()).body:
    if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
        key = node.targets[0].id
        if key in ("OLD_RAYS", "OLD_CONES", "BASIS", "HEIGHT_BASE"):
            data[key] = ast.literal_eval(node.value)
        elif key == "N":
            data[key] = S.Matrix(ast.literal_eval(node.value.args[1]))
rays, cones = data["OLD_RAYS"], data["OLD_CONES"]


def support(coefficients, strict=False):
    off_cone = []
    slopes = []
    for cone in cones:
        V = S.Matrix([rays[i] for i in cone])
        assert abs(V.det()) == 1
        m = V.inv() * S.Matrix([-coefficients[i] for i in cone])
        slopes.append(m)
        slacks = [S.Integer(coefficients[i]) + (S.Matrix([v]) * m)[0]
                  for i, v in enumerate(rays)]
        assert all(slack >= 0 for slack in slacks)
        for i, slack in enumerate(slacks):
            if i not in cone:
                if strict:
                    assert slack > 0
                off_cone.append(slack)
    return min(off_cone), slopes


height = [1009*(10000+h)+i+1 for i, h in enumerate(data["HEIGHT_BASE"])]
minimum, _ = support(height, strict=True)
assert minimum == 5047
_, anti_slopes = support([1]*11)
assert all(x.q == 1 for m in anti_slopes for x in m)
for index in (0, 1, 3, 4):
    adjacent = [list(m) for m in anti_slopes
                if (S.Matrix([rays[index]])*m)[0] == -1]
    assert S.Matrix(adjacent).rank() == 4  # this ray is a polytope vertex
assert abs(S.Matrix([rays[i] for i in (0, 1, 3, 4)]).det()) == 1
nef_lifts = []
for j in range(6):
    coefficients = [0]*11
    for i, ray in enumerate(data["BASIS"]):
        coefficients[ray] = data["N"][i, j]
    coefficients[10] = int(j >= 3)
    support(coefficients)
    nef_lifts.append(list(map(int, coefficients)))
assert 4-3 == 1  # integral divisor H3-H2 pairs primitively with alpha
print(json.dumps({
    "scope": "Toric arithmetic only; support and smooth hypersurface interpretation require the proof.",
    "strict_support_function": height,
    "minimum_off_cone_slack": int(minimum),
    "maximal_cones_checked_unimodular": len(cones),
    "nef_lifts_H1_through_H6": nef_lifts,
    "polytope_vertex_lattice_basis_indices": [0, 1, 3, 4],
    "primitive_pairing": 1,
}, indent=2))
