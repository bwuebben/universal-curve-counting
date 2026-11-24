"""Replay the numerical entry points cited by the five-article curve-counting series."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

PROJECT = Path(__file__).resolve().parent
SHARED = PROJECT / "supplementary_calculations"
PYTHON = [sys.executable]
SYMPY = ["uv", "run", "--no-project", "--with", "sympy==1.14.0", "python"]
SAGE = ["sage", "-python"]


def commands(group, temporary):
    checks = PROJECT / "checks"
    if group == "I":
        return [
            (SAGE, checks / "x9_source_genus_zero_box.py", ["--source-root", "vendor/global-higgs-branch", "--output", str(temporary / "source.json")]),
            (SAGE, checks / "x9_lift_bound.py", ["--output", str(temporary / "lift.json")]),
            (SAGE, checks / "x9_contracted_support.py", ["--output", str(temporary / "support.json")]),
            (PYTHON, checks / "x9_universal_line.py", ["--output", str(temporary / "line.json")]),
            (PYTHON, checks / "x9_plane_geometry.py", []),
            (PYTHON, checks / "x9_source_chamber.py", []),
            (SYMPY, checks / "local_series_comparison.py", []),
            (SYMPY, checks / "cubic_genus_ten_modes.py", []),
            (SYMPY, checks / "integrality_counterexample.py", ["--output", str(temporary / "counter.json")]),
            (SYMPY, checks / "integrality_counterexample_geometry.py", []),
        ]
    if group == "II":
        names = (
            "local_curve_degree_three degree_three_quotient_check integrality_counterexample "
            "integrality_counterexample_geometry local_series_comparison partial_integrality "
            "cubic_genus_ten_modes cubic_modal_eigenbranches cubic_szego_construction genus_five_cubic_correction"
        ).split()
        jobs = [(SYMPY, checks / (name + ".py"), []) for name in names]
        jobs.append((SYMPY, checks / "cubic_denominator_rows.py", ["--max-power", "104"]))
        for name in ("prime_five_genus_two", "quartic_weight_one_residues", "sextic_weight_one_residues"):
            jobs.append((SAGE, checks / (name + ".py"), []))
        jobs.append((SAGE, checks / "x9_contracted_support.py", ["--output", str(temporary / "ii-support.json")]))
        return jobs
    if group == "III":
        names = (
            "elliptic_face_indices quaternionic_horizontal_orbits etale_spin_characters "
            "etale_elementary_exponents fixed_excess_denominators dihedral_wall_transfer "
            "polyhedral_wall_transfer intersecting_unramified_walls etale_transport_and_prime_powers etale_cubic_transfer"
        ).split()
        return [(PYTHON, checks / (name + ".py"), []) for name in names] + [
            (SYMPY, checks / "etale_double_wall.py", [])
        ]
    return [(SYMPY, SHARED / "code" / name, []) for name in (SHARED / "arithmetic-checks.txt").read_text().splitlines() if name]


def verify(group, work):
    jobs = commands(group, work)
    for interpreter, script, args in jobs:
        log = work / f"{group}-{script.stem}.log"
        with log.open("w") as output:
            result = subprocess.run([*interpreter, str(script), *args], cwd=PROJECT,
                                    env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                                    stdout=output, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{script.name} failed:\n" + log.read_text()[-3000:])
        print(f"{group} {script.name} passed", flush=True)
    if group == "I":
        evidence = PROJECT / "evidence"
        for name, archived in [
            ("source.json", "x9-source-genus-zero-box.json"),
            ("lift.json", "x9-lift-bound-certificate.json"),
            ("support.json", "x9-contracted-support-certificate.json"),
            ("line.json", "x9-universal-line-certificate.json"),
            ("counter.json", "integrality-counterexample-certificate.json"),
            ("I-x9_plane_geometry.log", "x9-plane-fan-certificate.json"),
            ("I-x9_source_chamber.log", "x9-source-chamber-certificate.json"),
            ("I-integrality_counterexample_geometry.log", "integrality-counterexample-geometry.json"),
        ]:
            if json.loads((work / name).read_text()) != json.loads((evidence / archived).read_text()):
                raise ValueError(f"Certificate mismatch: {name}")
    print(f"{group}: {len(jobs)} completed", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("groups", nargs="*", help="I, II, III, IV-V (default: all)")
    parser.add_argument("--logs", type=Path, help="Keep individual command logs and generated certificates here")
    options = parser.parse_args()
    groups = options.groups or ["I", "II", "III", "IV-V"]
    if any(group not in ("I", "II", "III", "IV-V") for group in groups):
        parser.error("Unknown computation group")
    with tempfile.TemporaryDirectory(prefix="ucc-exact-") as temporary:
        work = options.logs or Path(temporary)
        work.mkdir(parents=True, exist_ok=True)
        for group in groups:
            verify(group, work)
