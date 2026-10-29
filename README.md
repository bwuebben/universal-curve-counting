# Universal curve counting

Two preprints by [Bernd Johannes Wuebben](https://orcid.org/0009-0002-4588-3284),
28 October 2026, with LaTeX sources and reproducible exact calculations.

1. **Universal curve counting under conifold transitions** — [PDF](paper/main.pdf),
   [source](paper/main.tex), [abstract and citation](https://bwuebben.github.io/papers/universal-conifold-counting/).
   18 pages. A rational universal specialization theorem for projective conifold
   transitions with smooth algebraic total space and torsion-free second homology,
   together with a complete compact calculation: 37 candidate coefficients,
   eight nonzero, and sum `3 x_{0,1}`. Integrality transfer is conditional.
2. **A counterexample to integral product exponents in universal curve counting** —
   [PDF](paper2/main.pdf), [source](paper2/main.tex),
   [abstract and citation](https://bwuebben.github.io/papers/universal-integrality-counterexample/).
   13 pages. A compact smooth Calabi–Yau threefold has a scalar local-curve product
   exponent with fractional part `2/3`. The paper also proves integrality in
   classes of divisibility at most two. The counterexample concerns the specified
   scalar generators; the integer elementary expansion is an input to the proof.

Read Paper 1 before Paper 2 for the transition setting. The counterexample paper
can also be read independently. Both are preprints; arXiv submission is pending.
See the [Universal curve counting programme](https://bwuebben.github.io/research/calabi-yau/#universal-curve-counting).

## Reproduce the calculations

The computations use Python 3, SymPy 1.14.0 and SageMath (verified with 10.9).
All finite inputs required by these checks are included. The two upstream toric
input files are byte-pinned under `vendor/global-higgs-branch/`.

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python reproduce.py
```

Alternatively, with `uv` and `sage` on the executable path:

```sh
uv run --no-project --with sympy==1.14.0 python reproduce.py
```

The runner compares all eight regenerated certificates with `evidence/`, checks
the agreement of independent matrix and quotient-algebra calculations, and
reconstructs Paper 1's printed tables. It uses a temporary directory and does not
replace the supplied certificates or manuscripts. `--sage /path/to/sage` selects
another SageMath executable. Individual scripts can also be run directly; for
`x9_source_genus_zero_box.py`, use `--source-root vendor/global-higgs-branch`.
The manuscript's original relative input path is replaced by this bundled path.

These checks establish finite rational identities, toric data and table entries.
The geometric support arguments, stable-pair moduli identification and imported
comparison theorems are given in the papers.

## Build the manuscripts

Each manuscript is independent. From `paper/` or `paper2/`, run:

```sh
latexmk -pdf main.tex
```

`paper/generate_data.py` regenerates the printed data appendix from `evidence/`.
The supplied PDFs are the release originals; rebuilding may change PDF metadata.
Check all supplied files with `shasum -a 256 -c SHA256SUMS`.

## Rights

See [RIGHTS.md](RIGHTS.md). The vendored inputs retain their MIT license.
