"""Reconstruct the complete finite X9 line table from the exact GW input.

Run: python3 checks/x9_universal_line.py --output certificate.json
Mathematical prerequisites: Paper 1, especially the all-genus PT
bound and primitive reconstruction. This script certifies only the finite
coefficient arithmetic; it does not prove those prerequisites.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pins = {
        "evidence/x9-source-genus-zero-box.json":
            "f9cc83fd2133ff6eebc621dfc48db95b44400492e5a570d760b34e6e114a8b2f",
        "evidence/x9-lift-bound-certificate.json":
            "e46c52b377512f9043606c49254892f081d1d5f4ea8a9b7650ce1305ac720121",
    }
    data = {}
    for path, digest in pins.items():
        raw = (root / path).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == digest, path
        data[path] = json.loads(raw)
    numerical = data["evidence/x9-source-genus-zero-box.json"]
    bound = data["evidence/x9-lift-bound-certificate.json"]
    sectors = {(1, 0): {}, (0, 1): {}, (1, 1): {}}
    for row in numerical["all_box_bps_and_gw"]:
        degree = tuple(row["degree"])
        bidegree = (degree[1], degree[2])
        if degree[3] == 0 and bidegree in sectors:
            value = Fraction(row["genus_zero_gw"])
            assert row["primitive"]
            if value:
                sectors[bidegree][degree] = value
    products = defaultdict(Fraction)
    for a, left in sectors[(1, 0)].items():
        for b, right in sectors[(0, 1)].items():
            products[tuple(x + y for x, y in zip(a, b))] += left * right
    products = {degree: value for degree, value in products.items() if value}
    candidates = [tuple(n) for n in bound["lift_candidates"]["1"]["candidates"]]
    assert len(candidates) == len(set(candidates)) == 37
    assert set(products).union(sectors[(1, 1)]) <= set(candidates)
    rows = []
    for degree in candidates:
        rows.append({
            "degree": degree,
            "x_coefficient": str(sectors[(1, 1)].get(degree, Fraction(0))),
            "x_squared_coefficient": str(products.get(degree, Fraction(0))),
        })
    sums = {
        "x": str(sum(sectors[(1, 1)].values(), Fraction(0))),
        "x_squared": str(sum(products.values(), Fraction(0))),
    }
    output = {
        "scope": "Finite arithmetic for Paper 1; geometric and all-genus inputs are separate proofs.",
        "input_sha256": pins,
        "generator": "x=x_{0,1}",
        "coefficient_formula": "A_11=x*P_11+x^2*P_10*P_01",
        "candidates": rows,
        "nonzero_count": sum(bool(Fraction(r["x_coefficient"]) or Fraction(r["x_squared_coefficient"])) for r in rows),
        "specialization_sum": sums,
        "lower_connected_sums": {str(k): str(sum(v.values(), Fraction(0))) for k, v in sectors.items()},
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"PASS: {len(rows)} candidates, {output['nonzero_count']} nonzero universal polynomials; sum {sums}.")


if __name__ == "__main__":
    main()
