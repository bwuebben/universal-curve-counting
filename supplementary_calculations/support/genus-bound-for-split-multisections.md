# A1. A genus bound for multisections of split real line equations

## Statement

Let $C$ be a compact Riemann surface of genus $g$, and let
$N=\bigoplus_{j=1}^s L_j$ be a sum of smooth complex line bundles.
Put $d_j=\deg L_j$ and $d=\max_j d_j$. Give each line a smooth real
Cauchy--Riemann operator $D_j$. The graph almost-complex structure on
$N$ is the one for which a graph over a holomorphic map $f$ solves
$f^*D_j u_j=0$ in each coordinate; pullback includes the differential
of $f$ on the $(0,1)$-form factor.

Let $X_0$ be a compact connected nodal complex curve and
$\nu:X=\coprod_{i=1}^k X_i\to X_0$ its normalization. Suppose that

1. $f:X\to C$ is holomorphic and nonconstant on every component,
   with total degree $m$;
2. $u_j\in\Gamma(X,f^*L_j)$ solves $f^*D_j u_j=0$;
3. the graph $F=(f,u_1,\ldots,u_s)$ matches at the two preimages
   of every node of $X_0$;
4. $F$ is generically one-to-one onto its image, and the images of
   distinct components of $X$ are distinct.

The last condition means that the support is reduced and its
parametrization has multiplicity one. It excludes a multiple cover
of a lower-degree graph. No immersion or transversality assumption
at a node is required for the bound.

**Proposition.** Under these hypotheses, if $b$ is the number of nodes
of $X_0$ and $R(f)=\sum_{x\in X}(e_x(f)-1)$ is the total ramification,
then

$$
                  R(f)+2b\ \le\ m(m-1)d.                 \tag{1}
$$

Consequently

$$
          p_a(X_0)\ \le\ 1+m(g-1)+\frac{m(m-1)}2d.       \tag{2}
$$

For balanced split lines, $d_j=g-1$, this becomes

$$
          p_a(X_0)\ \le\ 1+\frac{m(m+1)}2(g-1).          \tag{3}
$$

In particular a genus-two base gives bounds $2,4,7,11,\ldots$ in
degrees $1,2,3,4,\ldots$. The degree-two case recovers the bound in
[the line-operator note](double-cover-line-operator-bound.md), while
allowing arbitrary branch profiles and arbitrary numbers of line
summands in higher degree.

## Proof by differences of sheets

Form the analytic fiber product $X\times_C X$, take its reduced
structure, and normalize it. Discard the components which are
diagonals of components $X_i$. Denote the remaining disjoint union
of smooth compact Riemann surfaces by $Z$. There are maps

$$
       p_1,p_2:Z\longrightarrow X,\qquad
       h=f p_1=f p_2:Z\longrightarrow C.
$$

All components of $Z$ map finitely onto $C$: a finite fiber product
of the finite maps in question has no component supported over a
point. Above a regular value of $f$, its points are exactly the
ordered pairs of distinct sheets. Therefore

$$
                 \sum_{Z_\alpha\subset Z}\deg h_\alpha
                         =m(m-1).                       \tag{4}
$$

On a connected component $Z_\alpha$ define the difference sections

$$
           v_{\alpha j}=p_1^*u_j-p_2^*u_j
                     \in\Gamma(Z_\alpha,h_\alpha^*L_j).
$$

The two bundles being subtracted are canonically the same pullback
from $C$. Functoriality of pullback, including its differential
factor, and real linearity give

$$
                       h_\alpha^*D_j v_{\alpha j}=0.     \tag{5}
$$

At least one difference is not identically zero on each
$Z_\alpha$. Otherwise the two distinct sheets represented by
$Z_\alpha$ would have identical graph values over an open set of
regular values. This contradicts generic injectivity, or the
distinctness of the images of two components of $X$.
Choose one such index $j(\alpha)$, once for each component.

A nonzero solution of a real line Cauchy--Riemann equation has
isolated zeros of positive integral order, whose total is the
degree of the line. Thus

$$
 \#\operatorname{zeros}(v_{\alpha,j(\alpha)})
       = d_{j(\alpha)}\deg h_\alpha
       \le d\deg h_\alpha.                              \tag{6}
$$

The zero count includes multiplicities. The coefficients in (5)
are smooth even at a branch point of $h_\alpha$; vanishing of the
pulled-back zeroth-order coefficient there causes no difficulty.
The input is the scalar similarity principle, for example
[Wendl, *Lectures on Symplectic Field Theory*, Theorem 2.5.3,
printed pp. 52--53](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf#page=64).
The total-order statement is its usual first-Chern-class consequence.

We next count forced zeros in (6). At a point $x\in X$ with
ramification index $e$, choose holomorphic coordinates in which
$f(z)=z^e$. Near $(x,x)$ the fiber product has equation

$$
                        z_1^e=z_2^e.
$$

It has $e$ distinct smooth branches $z_2=\zeta z_1$, with
$\zeta^e=1$. One is the discarded diagonal. Each of the other
$e-1$ branches gives a distinct point of the normalization $Z$.
At every such point all difference sections vanish. Whichever
nonzero coordinate was chosen on its global component therefore
has a zero of order at least one. The contribution from these
points is at least $e-1$. Distinct ramification points give
distinct points in this count. Summing gives at least $R(f)$.

Now let a node of $X_0$ have distinct normalization preimages
$x,y$. Matching gives $F(x)=F(y)$, hence $f(x)=f(y)$. The ordered
pairs $(x,y)$ and $(y,x)$ lie off the diagonal. There is at least
one point of the normalization over each pair: locally the fiber
product has an equation $z_1^{e_x}=z_2^{e_y}$ and has
$\gcd(e_x,e_y)\ge1$ branches. Each point lies on a component
retained in $Z$, and every difference vanishes there. The two
ordered pairs thus contribute at least two zeros. They cannot
coincide with a point counted over $(x,x)$, and different nodes
have different pairs of normalization preimages. This contributes
at least $2b$ further zeros.

Adding (6) over all components and using (4) proves (1).
Riemann--Hurwitz on the components and connectedness of the dual
graph give

$$
 \begin{aligned}
 R(f)&=2\sum_i g(X_i)-2k-m(2g-2),\\
 p_a(X_0)&=\sum_i g(X_i)+b-k+1.
 \end{aligned}
$$

It follows that
$R(f)+2b=2p_a(X_0)-2-m(2g-2)$, which proves (2).
Substitution of $d=g-1$ proves (3).

## Sharpness and relation to adjunction

The coefficient of $m(m-1)d$ cannot be lowered in this class.
Fix $m\ge2$ and take a holomorphic line $L$ of positive degree $d$ for which
$L^m$ has a section $s$ with a nonempty reduced zero divisor.
The spectral curve

$$
                         w^m=s
$$

in the total space of $L$ is smooth and connected. Smoothness
follows at a zero of $s$ from its nonzero derivative; elsewhere
the derivative with respect to $w$ is nonzero. Local monodromy
around a simple zero is an $m$-cycle, so the cover is connected.
Its $md$ branch values each have ramification index $m$; hence
$R(f)=m(m-1)d$. Its tautological graph is embedded and attains
(2). Such sections exist, for example, whenever $md\ge2g+1$,
by very ampleness and a general choice of section. Extra line
coordinates may be set equal to zero.

For a single line, (2) is also the usual adjunction bound in the
four-dimensional total space: the class $m[C]$ has square
$m^2d$ and first Chern number $m(2-2g+d)$. See
[Wendl, *Holomorphic Curves in Low Dimensions*, Theorems 2.49 and 2.51](https://www.mathematik.hu-berlin.de/~wendl/pub/rationalRuled.pdf)
for intersection positivity and adjunction, respectively.
The proof above does not project the full graph to a single
simple curve. A coordinate projection can be multiply covered;
the nonzero coordinate may vary with the component $Z_\alpha$.
This is why the argument applies to a split bundle of any rank
without a simplicity hypothesis on any one coordinate projection.
No claim of literature priority is made here.

## Fiberwise nonlinear equations with independent coordinates

The same bound holds for the following larger class. Fix complex
linear Dolbeault operators $\bar\partial_{L_j}$ and smooth
fiberwise maps

$$
 A_j:L_j\longrightarrow
              \Lambda^{0,1}T^*C\otimes L_j
$$

covering the identity of $C$. The graph equations are

$$
                  \bar\partial_{f^*L_j}u_j+f^*A_j(u_j)=0. \tag{7}
$$

Here $f^*A_j$ again pulls back the form factor. The equation for
coordinate $j$ depends only on that coordinate. The maps $A_j$
need not be linear, and need not vanish at the zero section.

On $Z_\alpha$, put $a=p_1^*u_j$ and $b=p_2^*u_j$.
The difference of the two equations (7) is a smooth real line
Cauchy--Riemann equation for $a-b$, because

$$
 A_j(a)-A_j(b)=B_j(a,b)(a-b),\qquad
 B_j(a,b)=\int_0^1(d_{\rm fib}A_j)_{b+t(a-b)}\,dt.        \tag{8}
$$

The integral is a real-linear map on the common fiber $L_j$;
the construction is intrinsic under linear changes of bundle
frame. Pulling back the form factor gives a smooth zeroth-order
coefficient on $Z_\alpha$. Consequently the same positivity and
degree count applies, and the proof of (1)--(2) is unchanged.
This extension does not include coupling between different
coordinates: in that case the difference in one coordinate
need not solve a homogeneous scalar equation.

## Compactness for the controlled family

There is a useful compactness consequence with explicit spatial
hypotheses. Suppose $g\ge2$, $d\ge0$, and a compact smooth parameter
family of the equations above is realized in a fixed tube in a closed
symplectic manifold. Suppose the corresponding ambient almost-complex
structures form a compact smooth family tamed by its fixed symplectic
form. Fix a compact subset $K$ of the tube and a degree bound $M$.
Consider connected reduced nodal graph maps satisfying the proposition,
with images in $K$, no contracted components, and degree at most $M$.

**Corollary.** Every sequence of these maps has a stable-map convergent
subsequence, after reparametrization. Every limiting connected stable
map has source arithmetic genus at most

$$
                   G=1+M(g-1)+\frac{M(M-1)}2d.
$$

**Proof.** The degree and arithmetic genus take only finitely many
values by (2); pass to a subsequence on which both are fixed, and on
which the parameters converge smoothly. The class in the tube is
$m[C]$, so the symplectic area is $m\int_C\omega$ and is bounded.
The nodal graph maps are stable: an automorphism of a nonconstant
component preserving a generically injective map is the identity.
Gromov compactness for stable maps of fixed genus and bounded area in
the closed ambient target now supplies the subsequence. Its image
remains in $K$, and its arithmetic genus is the fixed source genus,
at most $G$. Equivalently, to reduce this to smooth-domain compactness,
the number of normalization components is at most $M$, and (1) bounds
the number of nodes. Pass to a fixed labeled dual graph and apply
marked Gromov compactness on each normalization component, treating
the node preimages as marked points. Continuity of evaluation permits
regluing each paired pair of marks in the limit; this preserves the
total arithmetic genus. The limit can have contracted components or multiple
covers, but their genera are included in this same bound. $\square$

The compactness theorem used in this corollary is the usual fixed-genus
theorem, not a new analytic compactness result. For the marked smooth-domain
theorem with varying tame structures and stable nodal limits see
[Wendl, *Holomorphic Curves in Low Dimensions*, Theorem 2.35,
printed p. 38](https://www.mathematik.hu-berlin.de/~wendl/pub/rationalRuled.pdf#page=48).
The extra input here is that the genus hypothesis has been deduced from
the split equations before applying that theorem. The spatial trapping
and common tame realization are assumptions, not consequences of (2).

## Consequences and limits for the comparison

For fixed $m,g$ and line degrees, the bound is uniform over all
smooth operators in the split class, and over all smooth independent
fiber equations (7). It bounds ramification and the number of nodal
identifications of reduced, multiplicity-one supports, not just
the genus of smooth double-cover graphs. Thus it supplies a concrete
genus estimate in arbitrary covering multiplicity for this class
of local geometries.

It does not by itself supply compactness in a noncompact tube:
sections of a linear equation can be scaled without bound. Nor
does it bound the arithmetic genus of stable maps with contracted
components, or of maps which cover a lower-degree support with
higher multiplicity. For example, a contracted curve of arbitrary
genus can be attached at a point of any fixed graph.

Finally, arbitrary perturbations used for transversality can couple
the line coordinates and destroy (5) or (8). A uniform bound in
the split class does not prove a uniform bound for these
perturbations. A theorem comparing this controlled class with
elementary-cluster factors, including contracted components and
the required transverse walls, is still needed before asserting
any vanishing of universal exponents. The separate
[wall comparison investigation](double-cover-wall-comparison.md)
examines that obstruction.
