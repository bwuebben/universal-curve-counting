# A8. First derivatives of scalar projections in triangular families

**Conclusion.** The first derivative of the scalar projection is transverse
to zero in the universal moduli space of marked, simple, noncore triangular
graphs with smooth domain. This holds in arbitrary covering degree and base
genus, under the parameter freedom specified below. The proof uses a local
deformation of the finite covering map, local interpolation for a linear
Cauchy--Riemann equation, and perturbations supported on an injective annulus.
It does not require the scalar projection to be immersed beforehand.

For balanced line degrees, the marked moduli over a path has Fredholm index
three. Its projected critical-point constraint has real codimension four.
Consequently generic triangular paths have immersed scalar projections on
all their simple noncore graphs. The result addresses the final analytic
qualification in [Calculation A2](all-genus-double-cover-comparison.md).
It does not prove an all-degree cluster formula, a statement about stable-map
boundaries, or transversality for multiple covers of the core.

## 1. Setting and the precise perturbation freedom

Let $C$ be a closed Riemann surface of genus $g$, and let $L_1,L_2$ be
complex line bundles. Fix a real Cauchy--Riemann operator $D_2$ satisfying

$$
                 \ker h^*D_2=0                                      \tag{1}
$$

for every nonconstant holomorphic map $h$ from a smooth closed connected
curve to $C$. Consider the equations

$$
 \bar\partial f=0,\qquad
 \bar\partial_{f^*L_1}u_1+f^*A(u_1)=0,\qquad
 f^*D_2u_2+f^*F(u_1)=0.                                               \tag{2}
$$

Here $A,F$ are smooth maps over $C$, with values in
$\Lambda^{0,1}T^*C\otimes L_1$ and
$\Lambda^{0,1}T^*C\otimes L_2$, respectively. Pullback of each form factor
uses $df$. Write $U=(f,u_1,u_2)$ and $v=(f,u_1)$. The first equation for
$u_1$ defines an almost complex structure $J_A$ on the surface $L_1$;
$v$ is $J_A$-holomorphic.

Assume $A(z,0)=F(z,0)=0$. Arbitrarily many prescribed fiber jets of $A,F$
along the zero section may be held fixed. The required parameter space must
contain every smooth perturbation $(\dot A,\dot F)$ with compact support
away from that section and inside the region where the triangular equations
are valid. It is enough to contain the local perturbations constructed below.
For paths, they may additionally be multiplied by a parameter bump supported
in the interior of the permitted parameter interval.

Prescribing jets does **not** mean forbidding perturbations on a fixed open
neighborhood of the entire core. The latter stronger restriction need not
provide the annuli required for curves lying arbitrarily close to the core.
Each perturbation used here vanishes on some neighborhood of the core, but
that neighborhood may depend on the curve and the perturbation.

The domain $Y$ is closed, connected and smooth. Its complex structure varies.
The map $f:Y\to C$ has positive degree $m$. We restrict to simple full graphs
$U$ whose image is not the core. The covering map $f$ itself need not be
simple as a holomorphic curve into $C$ and need not have simple ramification.

The scalar-simplicity argument in Proposition 1 of the
[triangular investigation](triangular-transversality.md)
applies in every degree: a simple full graph has simple scalar projection.
Indeed, on a component of the normalized off-diagonal fiber product, equality
of the scalar coordinates makes the second-coordinate difference a section
in the kernel of a pullback of $D_2$. Equation (1) makes that difference zero,
contradicting simplicity of the full graph. Thus no immersion assumption is
being inserted at this point.

## 2. An injective annulus exists even at a cusp

**Lemma 1.** Let $q\in Y$. For sufficiently small $\rho>0$, every point of
the punctured coordinate disk $0<|w|<\rho$ centered at $q$ is a globally
injective point of $v$, is unramified for $f$, and has $u_1(w)\ne0$.

**Proof.** A simple closed pseudoholomorphic curve is embedded away from
finitely many critical points and self-intersections. Apply this to $v$,
and shrink the disk to remove all exceptional domain points other than
possibly $q$. The ramification points of the nonconstant holomorphic map
$f$ are finite. Finally $u_1$ is nonzero and has isolated zeros: using
$A(z,0)=0$, the integral of the fiber derivative of $A$ along the segment
from $0$ to $u_1$ rewrites its equation as a homogeneous scalar real
Cauchy--Riemann equation. If $u_1$ were identically zero, (1) and the
second equation would force $u_2=0$, contrary to the noncore assumption.
Shrink the disk once more. $\square$

The finite exceptional-set statement is
[Wendl, *Lectures on Holomorphic Curves*, Corollary 2.15.4](https://www2.mathematik.hu-berlin.de/~wendl/pub/jhol_bookv33_2.pdf).
Its local factorization input is Theorem 2.14.7. In particular a nontrivial
local covering of a branch germ would give infinitely many noninjective
points and cannot occur for the simple closed map used here. A primitive
cusp is allowed: its derivative may vanish at its center, while its punctured
disk supplies the annulus.

For later extension of perturbations, take a compact annulus $K$ inside a
slightly larger annulus $V\Subset\{0<|w|<\rho\}$. Global injectivity gives

$$
 v(K)\cap v(Y\setminus V)=\varnothing.
$$

These are disjoint compact sets. Hence there is a target neighborhood of
$v(K)$ meeting the curve only in $v(V)$ and disjoint from the zero section.
Smooth values prescribed along $v(V)$ with support in $v(K)$ extend by a
tubular neighborhood and cutoff to permitted $\dot A,\dot F$.

## 3. Local interpolation for the vertical equation

**Lemma 2.** Let

$$
 D\xi=\partial_{\bar w}\xi+a(w)\xi+b(w)\overline\xi
$$

be a smooth real Cauchy--Riemann operator on a trivial rank-$r$ complex
bundle near $0$. After shrinking the disk, for every smooth $h$ and every
$c,d\in\mathbb C^r$ there is a smooth solution of

$$
                 D\xi=h,\qquad \xi(0)=c,
                 \partial_w\xi(0)=d.                                \tag{3}
$$

Only the holomorphic derivative is freely prescribed. The equation
determines the antiholomorphic derivative at the center.

**Proof.** Fix $p>2$. On the unit disk the map

$$
 W^{2,p}\longrightarrow W^{1,p}\oplus\mathbb C^r\oplus\mathbb C^r,
 \qquad \zeta\longmapsto
 (\partial_{\bar z}\zeta,\zeta(0),\partial_z\zeta(0))                   \tag{4}
$$

has a bounded right inverse. Start with a bounded right inverse
$T:W^{1,p}\to W^{2,p}$ for $\bar\partial$ on a disk. For data $(k,c,d)$,
replace $Tk$ by

$$
                  Tk-Tk(0)-z\partial_zTk(0)+c+zd.
$$

Point and derivative evaluation are bounded on $W^{2,p}$. Rescaling
$w=\epsilon z$ replaces $D$ by

$$
 \partial_{\bar z}+\epsilon a(\epsilon z)
            +\epsilon b(\epsilon z)\,\mathrm{conjugation}.
$$

The zeroth-order addition has operator norm $O(\epsilon)$ from
$W^{2,p}$ to $W^{1,p}$. Thus the augmented map remains split-surjective
for small $\epsilon$, by the Neumann-series argument for a bounded right
inverse. The derivative datum in rescaled coordinates is $\epsilon d$,
which is still arbitrary. Interior elliptic regularity makes the solution
smooth. $\square$

This is the first-derivative version of the local right-inverse argument in
[Wendl, *Lectures on Symplectic Field Theory*, Lemma 2.5.1](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf).
The additional interpolation condition is proved explicitly above; value
interpolation alone would not establish the claim needed here.

## 4. Moving the covering derivative

Put

$$
 R=2h-2-m(2g-2),\qquad h=g(Y).
$$

The stack of finite maps from smooth genus-$h$ curves to the fixed $C$,
with arbitrary ramification allowed, is smooth of complex dimension $R$.
Its tangent space at $f$ is $H^0(Y,N_f)$, where

$$
 0\longrightarrow T_Y\xrightarrow{df}f^*T_C
      \longrightarrow N_f\longrightarrow0.                           \tag{5}
$$

The sheaf $N_f$ is torsion of length $R$, so $H^1(Y,N_f)=0$.
This is the nondegenerate-map deformation theorem of
[Horikawa, *On deformations of holomorphic maps*, Theorems 1--2](https://www.jstage.jst.go.jp/article/pjab1945/48/2/48_2_52/_pdf/-char/ja),
with the full treatment in
[Horikawa, *On deformations of holomorphic maps I*](https://doi.org/10.2969/jmsj/02530372).
Finite deck groups give orbifold charts; the full simple graph has trivial
automorphism group. No fixed ramification-profile stratum is used.

At a ramification point choose coordinates with $f(w)=w^e$, $e\ge2$.
Then

$$
                    (N_f)_q\cong\mathbb C\{w\}/(w^{e-1}).
$$

For $e\ge3$, the class of $\alpha w$ supplies an arbitrary first-order
change $\dot f_w(0)=\alpha$. For $e=2$, this class is zero in the
unmarked normal sheaf: the missing freedom is motion of the marked point,
since $f''(0)\dot q=2\dot q$. Thus the derivative map on the *marked*
covering space has its expected transverse zero even at a simple branch.

The following concrete realization is useful because it keeps the ensuing
forcing local. Choose $\chi$ supported in a small coordinate disk and equal
to one near $0$. For real $s$ set, in the target coordinate,

$$
                    f_s(w)=f(w)+s\alpha\chi(w)w,                       \tag{6}
$$

and set $f_s=f$ outside that disk. On the inner disk $f_s$ is holomorphic
for the original $j$. On the transition annulus, $df$ is invertible, and
for small $s$, $f_s$ is an orientation-preserving local diffeomorphism.
Define

$$
                   j_s=(df_s)^{-1}j_Cdf_s
$$

on that annulus, and $j_s=j$ elsewhere. These definitions agree smoothly
near its boundary, where $f_s$ is holomorphic. Every almost complex
structure on a real surface is integrable. Therefore $(j_s,f_s)$ is an
actual family of degree-$m$ holomorphic maps, unchanged off the chosen
disk, with $j_s=j$ near $q$ and $\dot f_w(q)=\alpha$.

This construction does not keep the local branching order fixed. It must
not be used inside a stratum that forbids splitting higher ramification.

## 5. Surjectivity of the projected first derivative

**Proposition 3.** At any solution of (2) with simple noncore graph and a
marked point $q$ satisfying $dv(q)=0$, the derivative of the scalar
holomorphic first derivative is surjective on the tangent space to the
universal marked solution space. Its target has complex dimension two.
All parameter variations used for this surjectivity preserve the prescribed
core jets.

**Proof.** Since $dv(q)=0$, necessarily $df(q)=0$, so use the coordinates
of Section 4. Trivialize both line bundles over the corresponding target
disk. In these coordinates the vertical equations have the form

$$
 \partial_{\bar w}u_1+a_1(f,u_1)\overline{f_w}=0,
 \qquad
 \partial_{\bar w}u_2+a_2(f,u_1,u_2)\overline{f_w}=0,                   \tag{7}
$$

where $a_2$ is affine real-linear in $u_2$. The reference line-bundle
connection terms are included in $a_1,a_2$.

Specify arbitrary $\alpha,\beta\in\mathbb C$. Construct the compactly
supported $(\dot f,\dot j)$ of (6) with $\dot f_w(q)=\alpha$.
Choose its support inside an inner disk. Holding $A,F$ fixed for the
moment, the vertical linearization is

$$
                  D^{\mathrm{vert}}\xi+K(\dot f,\dot j)=0,             \tag{8}
$$

where $D^{\mathrm{vert}}$ is a smooth rank-two real Cauchy--Riemann
operator, triangular in these trivializations. The smooth forcing $K$
is supported in the inner disk. Lemma 2 solves (8) on a slightly larger
disk with

$$
                      \xi(q)=0,\qquad \partial_w\xi_1(q)=\beta.
$$

There is no restriction on $\alpha$ coming from this local equation.
For example (7) determines $\partial_{\bar w}\xi_1(q)$ from $\alpha$
and the values at $q$; it leaves $\partial_w\xi_1(q)$ free.

Multiply $\xi$ by a second cutoff $\rho$, equal to one on a neighborhood
of the support of $(\dot f,\dot j)$ and near $q$, and supported in the
larger disk. Extend $\rho\xi$ by zero. The residual in (8) is precisely

$$
                    r=[D^{\mathrm{vert}},\rho]\xi,                    \tag{9}
$$

supported in an outer annulus where every point is globally injective
for $v$, where $df\ne0$, and where $u_1\ne0$, by Lemma 1. Along this
annulus prescribe

$$
      \dot A(f(w),u_1(w))=-\frac{r_1(w)}{\overline{f_w(w)}},
 \qquad
      \dot F(f(w),u_1(w))=-\frac{r_2(w)}{\overline{f_w(w)}}.            \tag{10}
$$

The notation suppresses the chosen form and line frames. Since $f_w$
does not vanish on this compact support, these are smooth data. Lemma 1
extends them to allowed parameters supported away from the core and
meeting the scalar curve only on the annulus. Their pullbacks cancel
(9) exactly. The resulting global variation is tangent to the universal
solution space and has

$$
                    (\dot f_w(q),\dot u_{1,w}(q))=(\alpha,\beta).
                                                                         \tag{11}
$$

At $dv(q)=0$, variations of frames or of the target almost complex
structure contribute no additional terms to the linearization of the
first derivative. Near $q$ the chosen $j_s$ is fixed. Thus (11) is
exactly the required vertical derivative in any local trivialization of
the complex rank-two jet bundle.

The construction was written using parametrized maps. To pass to a
Hurwitz or Teichmüller slice, subtract the infinitesimal reparametrization
generated by a vector field $X$ and transform the marked point by the
same change of parametrization. The derivative contribution
$-d^2v(q)X(q)$ to the map variation is canceled by
$d^2v(q)X(q)$ from the marked-point variation. Frame changes multiply
the derivative by invertible maps and have no inhomogeneous term at its
zero. Hence (11) descends to the marked quotient. This explicitly retains
the marked-point freedom needed when $e=2$. $\square$

If $u_1(q)=0$, no parameter is varied at that point. Its zero is isolated,
and the outer annulus has positive distance from the zero section.
Consequently (10) preserves even the entire formal jet of the parameters
along the core. No division by $u_1(q)$ occurs. Likewise, higher cusp
orders cause no loss: the only denominator is $f_w$ on a compact annulus,
not at the critical point.

## 6. Fredholm and regularity details

Fix a compact scalar disk bundle, fix the coefficients on a collar of its
outer boundary, and restrict attention to curves whose scalar images lie
in the interior of the complement of that collar.
For paths, take smooth coefficients on its product with $[0,1]$, with the
prescribed core jets and fixed endpoint collars. These constraints define
a closed affine subspace $\mathcal P^\infty$ of the natural $C^\infty$
Fréchet space, equipped with the seminorms of uniform convergence of all
derivatives on that fixed compact domain. Thus $\mathcal P^\infty$ is a
Baire space. The same statement holds if all core jets are prescribed,
since this imposes a countable collection of closed linear conditions on
the difference from a reference family. Any required open tameness or
smallness condition is imposed within this space.

The ambient parameter space is not the union of spaces of perturbations
vanishing on unspecified neighborhoods of the core. Such perturbations
are only the local directions used in the surjectivity proofs; their
vanishing neighborhoods may vary. General elements of
$\mathcal P^\infty$ need satisfy only the stated jet and collar constraints.

The preceding tangent computation is usable because the universal
solution space is smooth even before scalar immersion has been proved.
Work over a smooth local chart of the finite-cover stack. At fixed $f$,
the vertical linearization is the elliptic operator

$$
 D^{\mathrm{vert}}:
 W^{1,p}(f^*(L_1\oplus L_2))
 \longrightarrow L^p(\Lambda^{0,1}T^*Y\otimes f^*(L_1\oplus L_2)).
$$

An adjoint cokernel section annihilating all parameter variations vanishes
on any globally injective, unbranched scalar disk off the core: arbitrary
values of $\dot A,\dot F$ give both output coordinates there. Adjoint unique
continuation then makes the section zero everywhere. The universal
vertical equation is therefore split-surjective. A finite collection of
smooth parameter variations spans the finite-dimensional cokernel of
$D^{\mathrm{vert}}$.

For derivative evaluation use $W^{k,p}$ charts with $k\ge3$, $p>2$,
and auxiliary local $C^\ell$ Banach parameter charts, with $\ell$
sufficiently high for the required Sard--Smale differentiability.
One may complete the smooth allowed directions in the $C^\ell$ norm;
the perturbations proving surjectivity remain smooth and preserve every
prescribed jet. The first-derivative evaluation with a
moving mark is then at least continuously differentiable. Increasing
$k$ and the coefficient regularity gives any prescribed finite
differentiability. One must not treat derivative evaluation as a bounded
functional on $W^{1,p}$.

Universal surjectivity at high regularity follows from the preceding
$W^{1,p}\to L^p$ argument and elliptic regularity: use the finite smooth
parameter directions spanning the cokernel, solve at low regularity,
and upgrade the solution for a $W^{k-1,p}$ right-hand side. The four real
directions in Proposition 3 are themselves smooth and lie in the kernel
of the universal equation. Adding them proves surjectivity of the
augmented equation-plus-jet map in these high-regularity charts. This
avoids a distribution argument involving derivatives of point masses.

For $\deg L_1=\deg L_2=g-1$, Riemann--Roch gives

$$
 \operatorname{ind}_{\mathbb R}D^{\mathrm{vert}}
       =4(1-h)+4m(g-1)=-2R.                                           \tag{12}
$$

The covering chart has real dimension $2R$. The projection from the
universal unmarked graph moduli to the static parameter space therefore
has Fredholm index zero. Adding a path parameter and one free marked
point gives index three. The complex rank-two projected derivative
constraint lowers the index to

$$
                              0+1+2-4=-1.                            \tag{13}
$$

Relative parametric transversality thus excludes scalar critical points
for generic paths on the parameter interval where perturbations are
allowed. A fixed endpoint collar must already have the desired property;
the positive-spin collars with no noncore graphs meet this condition.
Countably many covering degrees, domain genera and local Fredholm charts
can be handled by countable smooth local perturbation arguments: in each
chart the finite collections of smooth permitted directions above supply
transversality, and the resulting open dense local conditions are
intersected in $\mathcal P^\infty$. This yields a residual subset in its
natural $C^\infty$ topology, rather than merely a residual set in one
finite-regularity completion. No compactness assertion follows from this
countable intersection.

For comparison, the general ambient first-jet transversality theorem is
[Oh--Zhu, *Embedding property of $J$-holomorphic curves in Calabi--Yau manifolds for generic $J$](https://arxiv.org/abs/0805.3581).
The construction above verifies the restricted triangular perturbation
freedom directly; genericity in the full space of almost complex structures
would not by itself imply this conclusion.
