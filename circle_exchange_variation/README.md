# The first variation and the one-mark coefficient of circle exchange for balanced local curves

Bernd Johannes Wuebben · November 26, 2026 · [Manuscript](main.pdf) · [LaTeX source](main.tex) · [Abstract and citation](https://bwuebben.github.io/papers/refined-circle-exchange-local-curves/) · [Zenodo](https://doi.org/10.5281/zenodo.24570158)

The two main statements compare the Picard-family stable-pair indices in arbitrary positive curve degree and genus. One fixes the order in the volume parameter, the other fixes the exchanged grading and retains the volume parameter. The paper follows [the two-circle comparison at volume one](../circle_exchange/README.md).

## Build the manuscript

From this directory, with TeX Live and latexmk installed, run

~~~sh
SOURCE_DATE_EPOCH=1795694400 FORCE_SOURCE_DATE=1 TZ=UTC latexmk -pdf -halt-on-error -interaction=nonstopmode main.tex
qpdf --check main.pdf
~~~

## Reproduce the calculations

The checks use Python 3, SymPy 1.14.0 and SageMath 10.9. From this directory run

~~~sh
python3 checks/walls.py
uv run --no-project --with sympy==1.14.0 python checks/genus_two.py
PYTHONDONTWRITEBYTECODE=1 sage -python checks/one_mark.py
~~~

- The wall check enumerates every partition through degree 13 and compares its hook lengths and legs with the contributing reflection intervals. It checks the weight and determinant sums and the cancellation under transposition.
- The genus-two check computes the first 20 ribbon coefficients, the fixed-divisor weights and the Adams character through divisor length seven.
- The one-mark check tests the cleared telescoping identity and endpoint substitutions, the nonlinear genus-two coefficient, and the individual removable-box residues. It compares 48 shapes at two exact rational triples, five shapes over the rational function field in two variables at fixed volume, and four shapes from independent height arrays up to Euler excess six.

These are finite checks of formulas in the proofs. The properness, specialization, integral descent and unrestricted statements depend on the mathematical arguments in the manuscript.
