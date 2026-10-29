"""Exact arithmetic certificate for Paper 2.

This checks the local trace and the finite F1 arithmetic. It does not replace
the proofs of compact support, PT scheme structure, or imported GW comparison.
Run with SymPy 1.14.0; optionally pass --output PATH for a JSON certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as S
from local_curve_degree_three import build, in_w, r, w, T


def decompositions(a, b, smallest=(0, 1)):
    """All unordered nonzero effective vector decompositions."""
    if a == b == 0:
        yield ()
        return
    for x in range(a + 1):
        for y in range(b + 1):
            if (x, y) == (0, 0) or (x, y) < smallest:
                continue
            for tail in decompositions(a - x, b - y, (x, y)):
                yield ((x, y),) + tail


def phi(a, b):
    return a * (a + 1) // 2 + b * (1 - a)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    _, _, _, _, _, _, B = build(3)
    charpoly = [in_w(c) for c in B.charpoly(T).all_coeffs()]
    f2 = S.expand(in_w(S.trace(B)))
    f3 = S.expand(in_w(S.trace(B**2)))
    assert f2 == 81*w**3 + 50*w**4 + w**5
    assert f3 == (2673*w**6 + 9612*w**7 + 2534*w**8
                  + S.Rational(292, 3)*w**9 + w**10)
    U, V = w**3*(4-w)*(3-w)**2, w**3*(3-w)**2
    delta2 = S.expand(f2 - (2*U + V) - w*(2*w**3))
    delta3 = S.expand(f3 - (2*U**2 + V**2)
                      - w**2*(24*w**5 - 2*w**6))
    assert delta2 == 120*w**4 - 20*w**5 + 2*w**6
    assert delta3.coeff(w, 9) == S.Rational(4720, 3)
    z = S.symbols("z")
    local_pt = S.expand(f3.subs(w, (1+z)**2/z))
    assert local_pt.coeff(z, -9) == S.Rational(352, 3)

    cases = list(decompositions(3, 4))
    split_cases = [parts for parts in cases if len(parts) >= 2]
    minimum = min(sum(phi(*part) for part in parts) for parts in split_cases)
    assert minimum == 0
    assert phi(3, 4) == -2
    h0 = sum(4 - i + 1 for i in range(4))
    assert h0 == 14
    signed_euler = (-1)**(h0 - 1) * h0
    assert signed_euler == -14
    shift = -delta3.coeff(w, 9) * signed_euler
    assert shift == S.Rational(66080, 3)
    assert shift - S.floor(shift) == S.Rational(2, 3)

    root = Path(__file__).resolve().parents[1]
    inputs = ["checks/local_curve_degree_three.py",
              "checks/degree_three_quotient_check.py",
              "evidence/x9-contracted-support-certificate.json"]
    result = {
        "scope": "Exact arithmetic only; geometric and virtual-class proofs are in Paper 2.",
        "variable": "w=2-r-r^-1; PT z=-r",
        "degree_three_scalar_characteristic_coefficients": list(map(str, charpoly)),
        "scalar_F2_degree_three": str(f2),
        "scalar_F3_degree_three": str(f3),
        "connected_difference_base_genus_two": str(delta2),
        "connected_difference_base_genus_three": str(delta3),
        "scalar_local_PT_z_minus_9": str(local_pt.coeff(z, -9)),
        "F1_class_s_f": [3, 4],
        "nontrivial_effective_decompositions_checked": len(split_cases),
        "minimum_sum_surface_chi_over_nontrivial_decompositions": minimum,
        "surface_curve_chi": phi(3, 4),
        "linear_system_dimension": h0 - 1,
        "signed_projective_space_euler": signed_euler,
        "primitive_divisor_pairing_H3_minus_H2": 4 - 3,
        "geometric_exponent_integer_shift": str(shift),
        "geometric_exponent_fractional_part": str(shift - S.floor(shift)),
        "input_sha256": {p: hashlib.sha256((root/p).read_bytes()).hexdigest()
                         for p in inputs},
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
