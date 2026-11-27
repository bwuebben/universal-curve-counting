# The quadratic and cubic coefficients of circle exchange for balanced local curves

Bernd Johannes Wuebben · November 27, 2026

[Manuscript](main.pdf) · [LaTeX source](main.tex)

Two external circle actions give different fixed threefolds over a balanced local curve. The paper compares their quadratic and cubic exchanged coefficients with the volume character free. Its local statement is a rational identity for each Young diagram. Descent gives the determinant-normalized integral family comparison over the Picard variety in every positive genus and curve degree.

## Build

From this directory run

~~~sh
SOURCE_DATE_EPOCH=1795755600 FORCE_SOURCE_DATE=1 latexmk -pdf -halt-on-error main.tex
~~~

The source uses an 11-point AMS article and standard TeX Live packages. The PDF trailer omits the path-dependent identifier.

## Exact calculations

The two programs use SageMath. Each implements the finite identity from the indicated proof.

~~~sh
sage -python checks/quadratic_interpolation.py
sage -python checks/cubic_slices.py --all-characters --workers 4
~~~

The quadratic computation checks the four-label identity at the off-diagonal pairs of five interpolation nodes over a rational function field. The cubic computation checks all 540 pairs of Macdonald parameters, seven slice ratios at each pair, and the required choices of six labels. The degree bounds and the repeated-step specialization are given in the manuscript. Running the cubic program without options checks one pair and reports the ranks of the seven row spaces. These calculations prove the stated rational identities. The contour calculation and the descent to the Picard variety are given separately in the text.

Copyright 2026 Bernd Johannes Wuebben. All rights reserved.
