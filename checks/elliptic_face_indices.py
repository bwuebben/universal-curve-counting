"""Check the quotient's section labels, Dold coefficients, and elliptic partition."""

from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
import re
import subprocess


def divisors(n):
    return [k for k in range(1, n + 1) if n % k == 0]


def check_indices(cutoff):
    labels = list(product((0, 1), repeat=4))
    assert len(set(labels)) == 16
    assert {tuple(2 * v % 2 for v in label) for label in labels} == {(0,) * 4}
    odd_sections = sum(sum(label) % 2 for label in labels)
    assert odd_sections == 8
    exponents = {}
    indices = {}
    for k in range(1, cutoff + 1):
        # Alternating trace of negation^k acting on exterior powers of H^1.
        indices[k] = sum((-1) ** (j * (1 + k)) * comb(4, j) for j in range(5))
        assert indices[k] == sum((-1) ** ((k + 1) * sum(v)) for v in labels)
        ds = divisors(k)
        residual = indices[k] - sum(d * exponents[d] for d in ds if d < k)
        exponent, remainder = divmod(residual, k)
        assert remainder == 0
        exponents[k] = exponent
        assert exponent == {1: 16, 2: -8}.get(k, 0)

        # log of the elliptic partition has coefficient sum_{a|k} 1/a.
        log_partition = sum(
            exponents[d] * sum(Fraction(1, a) for a in divisors(k // d))
            for d in ds
        )
        assert log_partition == sum(Fraction(indices[a], a) for a in ds)


def check_printed_formula():
    pdf = Path(__file__).resolve().parents[1] / "unramified_wall_crossing" / "main.pdf"
    body = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True)
    prose = " ".join(body.rsplit("References", 1)[0].split())
    patterns = [
        r"sixteen classes βv are distinct and generate Γ0 . They satisfy 2βv = 2β0",
        r"exponents at degrees one and two are 16 and −8",
        r"fixed point indices 16 on the odd iterates and 0 on the even iterates",
        r"F1\s*\(Q2β0\s*\)−8\s*F1\s*\(Qβv\s*\)",
    ]
    for pattern in patterns:
        if re.search(pattern, prose, re.IGNORECASE) is None:
            raise AssertionError(f"Missing quotient calculation in the PDF: {pattern}")


if __name__ == "__main__":
    check_indices(128)
    check_printed_formula()
    print("Verified 16 torsion labels, both printed exponents, and 128 iterate and partition coefficients.")
