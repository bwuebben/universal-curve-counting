# The topological string on balanced local curves from two circle actions on a Calabi–Yau fivefold

Bernd Johannes Wuebben · November 25, 2026 · Revised November 26, 2026 · [Publication and citation](https://bwuebben.github.io/papers/circle-exchange-local-curves/) · [Zenodo](https://doi.org/10.5281/zenodo.24549724).

[Read the manuscript](main.pdf) · [LaTeX source](main.tex) · [Bibliography](references.bib)

The central theorem compares two external circles of a local Calabi–Yau fivefold over a fixed smooth positive-genus curve. Its new content is the integral comparison on the common space of cycles, in every degree and across the balanced Picard family. The known fixed-target hook formula enters once, to evaluate the coefficient in the original presentation. The rotated calculation derives it on the surviving plane Hilbert schemes. In genus two the applications identify the connected class on the collision diagonal and evaluate a first variation of its underlying degree-two family index at the lowest grading.

Theorem 1.2 gives the integral first variation at the lowest original grading for every positive exchanged grading. At the second exchanged grading its Chern character is twelve times the principal polarization of the Picard variety.

## Build

Run from this directory with TeX Live and latexmk.

~~~sh
latexmk -pdf -halt-on-error -interaction=nonstopmode main.tex
~~~

The JHEP style and bibliography files are included. The title, contents and introduction begin on separate pages.

From the project root, run the two exact-algebra checks.

~~~sh
uv run --no-project --with sympy==1.14.0 python checks/rotated_circle_and_collision.py
uv run --no-project --with sympy==1.14.0 python checks/first_volume_variation.py
~~~

The checks cover split levels through curve degree seven, connected Euler products through degree six, and the first-variation coefficients through degree twenty in the second counting variable. The all-degree statements rely on the proofs in the manuscript.
