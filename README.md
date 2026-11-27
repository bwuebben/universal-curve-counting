# Universal curve counting

Preprints by [Bernd Johannes Wuebben](https://orcid.org/0009-0002-4588-3284),
with LaTeX sources and reproducible exact calculations.

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
3. **Antidiagonal local curves and integral exponents in universal curve counting** —
   [PDF](paper3/main.pdf), [source](paper3/main.tex),
   [abstract and citation](https://bwuebben.github.io/papers/antidiagonal-universal-curve-counting/).
   Assuming Pardon's universal Gromov–Witten evaluation for analytic families,
   the balanced antidiagonal local factors give a unique product with integral
   elementary-cluster exponents indexed by effective classes in full integral
   homology on every smooth projective Calabi–Yau threefold. The cluster
   evaluation characterizes the compatible universal factor, and agreement
   in degree two forces the splitting degrees to be balanced.
4. **A geometric product formula for universal curve counting on Calabi–Yau threefolds** — [PDF](paper4/main.pdf), [source](paper4/main.tex), [abstract and citation](https://bwuebben.github.io/papers/geometric-universal-curve-product/).
   For smooth projective threefolds with trivial canonical bundle and vanishing first and second structure-sheaf cohomology, the balanced factors give the universal series in Pardon’s rational cycle algebra. Deformation and excision identify the exponents with signed elementary-cluster multiplicities.

Read Paper 1 before Paper 2 for the transition setting. The counterexample paper
can also be read independently. Paper 3 returns to the same compact series and compares the two choices of local factors. Paper 4 gives the geometric proof of the antidiagonal product.

Zenodo archives and current editions:

- [Conifold transitions](https://doi.org/10.5281/zenodo.23750046) · [edition 3.0](https://doi.org/10.5281/zenodo.24433173).
- [Integrality counterexample](https://doi.org/10.5281/zenodo.23750048) · [edition 3.0](https://doi.org/10.5281/zenodo.24433174).
- [Antidiagonal local curves](https://doi.org/10.5281/zenodo.24432700) · [edition 2.0, revised 22 November 2026](https://doi.org/10.5281/zenodo.24447953).
- [Geometric product formula](https://doi.org/10.5281/zenodo.24455518) · [edition 1.0](https://doi.org/10.5281/zenodo.24455519).

The 21 November revision of the conifold and counterexample papers supplies Zenodo DOIs for the author’s cited preprints.
BibTeX citations are in [CITATION.bib](CITATION.bib).
See the [Universal curve counting programme](https://bwuebben.github.io/research/calabi-yau/#universal-curve-counting).

## Integral generators and local curve counts

1. **Integral generators for curve counting on threefolds**. Opposite circle weights give explicit integral polynomial generators and a unique product on the zero-Chern face of a compact Kähler threefold with nonnegative first Chern class on curves. [PDF](integral_factorization/main.pdf) · [LaTeX](integral_factorization/main.tex) · [Abstract and citation](https://bwuebben.github.io/papers/integral-generators-threefold-curves/) · [Zenodo](https://doi.org/10.5281/zenodo.24510962).
2. **Triple curves and unbounded denominators in threefold curve counting**. Scalar exponents have unbounded denominators on one smooth projective Calabi–Yau threefold, even in classes of divisibility three. [PDF](scalar_denominators/main.pdf) · [LaTeX](scalar_denominators/main.tex) · [Abstract and citation](https://bwuebben.github.io/papers/triple-curves-unbounded-denominators/) · [Zenodo](https://doi.org/10.5281/zenodo.24510966).
3. **Unramified covers and wall crossing for local Gromov–Witten invariants**. At two independent double-cover walls the alternating change in all source genera is the local invariant of their common fourfold cover. [PDF](unramified_wall_crossing/main.pdf) · [LaTeX](unramified_wall_crossing/main.tex) · [Abstract and citation](https://bwuebben.github.io/papers/unramified-local-gw-wall-crossing/) · [Zenodo](https://doi.org/10.5281/zenodo.24510972).
4. **Branched triple covers and wall crossing for local curve counts**. A branched cubic crossing changes the genus-five multiplicity while all unramified determinant signs through degree seven agree. [PDF](ramification_and_compact_counts/main.pdf) · [LaTeX](ramification_and_compact_counts/main.tex) · [Abstract and citation](https://bwuebben.github.io/papers/branched-triples-local-wall-crossing/) · [Zenodo](https://doi.org/10.5281/zenodo.24510974).
5. **Spin structures and local Gromov–Witten invariants of curves**. For positive maximal flat symplectic normal holonomy, a spin gluing sign determines the partition sum in every degree and source genus. [PDF](spin_partition_functions/main.pdf) · [LaTeX](spin_partition_functions/main.tex) · [Abstract and citation](https://bwuebben.github.io/papers/spin-local-gw-of-curves/) · [Zenodo](https://doi.org/10.5281/zenodo.24510976).

The [shared supplement](supplementary_calculations/README.md) contains the calculations for the branched-cover and spin papers. For the numerical checks accompanying the five papers, run `python3 verify_series.py` with SageMath and uv available.

## Topological strings and curve counting

**The topological string on balanced local curves from two circle actions on a Calabi–Yau fivefold**. On a balanced local curve the variable grading the Euler characteristic of a stable pair becomes a plane rotation weight when the external circle of the fivefold is exchanged. The equality holds on the common space of cycles in every positive curve degree and throughout the Picard family. The November 26 revision states the genus-two first variation at the lowest original grading as Theorem 1.2. [PDF](circle_exchange/main.pdf) · [Source and calculations](circle_exchange/README.md) · [Abstract and citation](https://bwuebben.github.io/papers/circle-exchange-local-curves/) · [Zenodo](https://doi.org/10.5281/zenodo.24549724).


## Refined circle exchange

**The first variation and the one-mark coefficient of circle exchange for balanced local curves**. For a balanced local curve of positive genus, the first variation of the two stable-pair indices agrees in integral Picard-family K-theory in every degree and both gradings. At unrestricted volume, their entire coefficient of exchanged grading one agrees for each Young diagram and is a sum over removable boxes. [PDF](circle_exchange_variation/main.pdf) · [Source and exact checks](circle_exchange_variation/README.md) · [Abstract and citation](https://bwuebben.github.io/papers/refined-circle-exchange-local-curves/) · [Zenodo](https://doi.org/10.5281/zenodo.24570158).

**The quadratic and cubic coefficients of circle exchange for balanced local curves**. The two modified stable-pair indices agree through exchanged degree three, over the Picard variety, for every positive genus and curve degree at arbitrary volume character. The disk comparison is proved separately for each Young diagram. [PDF](circle_exchange_coefficients/main.pdf) · [Source and exact calculations](circle_exchange_coefficients/README.md) · [Paper page](https://bwuebben.github.io/papers/quadratic-cubic-circle-exchange/) · [Zenodo](https://doi.org/10.5281/zenodo.24579580).

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

These checks establish finite rational identities, toric data and table entries.
The geometric support arguments, stable-pair moduli identification and imported
comparison theorems are given in the papers.

## Build the manuscripts

Each manuscript is independent. From `paper/`, `paper2/`, `paper3/`, or `paper4/`, run:

```sh
latexmk -pdf main.tex
```

`paper/generate_data.py` regenerates the printed data appendix from `evidence/`.
The supplied PDFs are the release originals; rebuilding may change PDF metadata.
Check all supplied files with `shasum -a 256 -c SHA256SUMS`.

## Rights

See [RIGHTS.md](RIGHTS.md). The vendored inputs retain their MIT license.
