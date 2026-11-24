# A3. Prime-degree genus cutoffs for split-spin local contributions

Let $C$ be a smooth curve of genus $g\ge2$, and let $a,b$ be
non-effective theta characteristics on $C$. Use positive spin operators
on both lines, so the embedded core has sign $+1$. Write $n_{m,h}$ for
the balanced elementary exponents of this local problem, as defined in
Article III, Proposition 3.1 and Corollary 3.2.

**Theorem.** For every prime $p$,

$$
 n_{p,h}=0\qquad\text{if}\qquad
 h>1+\binom{p+1}{2}(g-1).                                      \tag{1}
$$

The lower bound is $n_{p,h}=0$ for $h<1+p(g-1)$, as already proved
in Article III, Proposition 3.1. In particular the only possible degree-three exponents
over a genus-two core have genera four through seven. For distinct
even spin lines on such a core, the complete connected degree-three
series is

$$
 [q^3]\mathsf L_{a,b}
 =\frac{28}{3}z^3-6z^4+(21+A)z^5+(B-2)z^6,
 \quad z=4\sin^2(t/2),\quad A=n_{3,6},\ B=n_{3,7}.                \tag{2}
$$

The integers $A,B$ are not determined here. The previously established
entries are $n_{3,4}=-64$ and $n_{3,5}=64$. Every even theta
characteristic in genus two is non-effective, so no additional
restriction is needed in this example. In higher genus the stated
non-effectivity hypothesis is retained.

The new analytic input is first-derivative transversality for the
scalar projection of a triangular graph. It extends the controlled
actual-curve argument of [Calculation A2](all-genus-double-cover-comparison.md)
to arbitrary fixed degree cutoffs. Prime-degree filtering then removes
proper-multiple correction clusters. The theorem is a uniform bound
on elementary support genera, not a vanishing theorem for source
genera of stable maps. All contracted components are retained.

## 1. A triangular path with fixed regular core

Set $c=g-1$. On smooth lines $L_1,L_2$ of degree $c$, fix a positive
spin operator $D_b$ on $L_2$. The spin/Stokes identity implies

$$
                 \ker f^*D_b=0                                    \tag{3}
$$

for every nonconstant finite holomorphic map from a smooth connected
curve to $C$. Sections over connected stable sources also vanish:
they vanish on the nonconstant normalized components, and matching
propagates that vanishing over the contracted components.

The equations for a graph $U=(f,u_1,u_2)$ are

$$
 \bar\partial f=0,\qquad
 \bar\partial_{f^*L_1}u_1+f^*A_t(u_1)=0,\qquad
 f^*D_bu_2+f^*F_t(u_1)=0.                                         \tag{4}
$$

The form factors in the pullbacks include $df$. The smooth maps $A_t,F_t$
are defined on $L_1$, with values in the $(0,1)$-form bundles on $C$
with coefficients in $L_1,L_2$, respectively. These equations define
an almost complex structure whose projections to $C$ and to the
almost complex surface $L_1$ are holomorphic. The fibers have their
standard complex structure. Hence a compact nonconstant curve in
the tube cannot be contained in a fiber.

Join the first-line operators $D_b$ and $D_a$ through invertible
operators. Lower the positive mass on $b$ to zero; connect its
Dolbeault structure to that of $a$ in
$\operatorname{Pic}^{g-1}(C)\setminus\Theta$; restore the mass on $a$.
The theta divisor has real codimension two, and its complement in
the connected complex torus is path connected. Non-effectivity
ensures invertibility at the two zero-mass endpoints. Positive mass
ensures injectivity on the other pieces, and the base index is zero.
Choose smooth representatives $D_{1,t}$ and constant endpoint collars.

Prescribe

$$
 A_t(0)=F_t(0)=0,\qquad
 d_{\rm fib}A_t|_0=D_{1,t}-\bar\partial_{L_1},\qquad
 d_{\rm fib}F_t|_0=0.                                            \tag{5}
$$

Keep the entire positive-spin equations fixed in the endpoint collars.
They have no noncore support. The core operator throughout the path
is $D_{1,t}\oplus D_b$, and is invertible with positive sign.

Work in the closed affine space of smooth coefficient paths on a
fixed compact disk bundle in $L_1$, subject to (5) and the fixed
endpoint collars. It carries its $C^\infty$ Fréchet topology and is
Baire. Use high-$C^\ell$ local Banach charts for transversality and
then the countable smooth local perturbation argument in this space.
The tangent space contains every smooth perturbation supported off
the zero section and in the interior parameter interval. Such local
directions do not define the ambient topology: one does not take
the union of spaces vanishing on unspecified open core neighborhoods.

The compact realization in the proof of Article III, Proposition 3.1, places this path on a common tube
$V$ in a closed symplectic six-torus, with one taming form, primitive
class $A_0=[C]$, and $c_1=0$. The construction uses compatibility on
the zero section and openness of tameness on a sufficiently small
tube, uniformly over the compact parameter interval. Perturbations
may be chosen small enough to retain it. Local curve classes are
$mA_0$, and degree is their positive integral projection degree.
Fix any finite cutoff $M$ and take
$M\int_C\omega\le\Lambda<(M+1)\int_C\omega$.

## 2. The genus estimate and simple scalar projections

For a reduced nodal graph with normalization $Y$ and total degree
$m$, normalize the off-diagonal fiber product $Y\times_CY$. On any
retained component the first-coordinate difference satisfies a
homogeneous scalar real Cauchy--Riemann equation, by integrating the
fiber derivative of $A_t$ between the two values.

If this first difference were identically zero, the difference of
the second equations would be a pullback of $D_b$ applied to the
second difference. Equation (3) would make that difference zero as
well. This contradicts the generic injectivity of the reduced full
graph on the chosen off-diagonal component. Thus the first difference
is nonzero on every retained component. The zero count and node count
of [Calculation A1](genus-bound-for-split-multisections.md) give

$$
 R(f)+2b\le m(m-1)c,\qquad
 p_a\le G_m:=1+\frac{m(m+1)}2c.                                  \tag{6}
$$

Here $b$ counts nodes, and $R$ is the total ramification of the
normalized covering. For a smooth simple graph of genus $h$,
Riemann--Hurwitz also gives $h\ge H_m:=1+mc$.

The same difference argument has two further consequences in every
degree. A simple full graph has simple scalar projection
$v=(f,u_1)$ to $L_1$. Two distinct simple full graphs have distinct
scalar images. A putative coincidence would give the same scalar
values on an off-diagonal fiber-product component, and then (3)
would identify the full graphs. These statements do not require
scalar immersion.

## 3. Removing critical points of the scalar projection

We give the essential argument; the full local proof and regularity
details are in
[the projected-jet investigation](triangular-projected-jet.md).
It is important to work before immersion in covering-map and vertical
section charts, rather than assert a locally free normal quotient.

The deformation space of a finite map $f:Y\to C$ with smooth domain
and arbitrary ramification is a smooth complex orbifold of dimension
$R=2h-2-m(2g-2)$. Its deformation sheaf is the torsion cokernel of
$T_Y\to f^*T_C$, of length $R$ and with zero first cohomology.
This is the curve case of
[Horikawa's nondegenerate-map deformation theorem](https://www.jstage.jst.go.jp/article/pjab1945/48/2/48_2_52/_pdf/-char/ja).
Finite deck transformations cause no stabilizer of a full simple graph.
Together with the vertical sections these charts have real Fredholm
index $2R-2R=0$, since

$$
 \operatorname{ind}_{\mathbb R}D^{\rm vert}
 =4(1-h)+4mc=-2R.                                                 \tag{7}
$$

Universal regularity already holds on the simple noncore locus.
A simple scalar curve has an injective unbranched disk off the core.
Arbitrary values of $(\dot A,\dot F)$ on that disk kill any adjoint
cokernel section there, and unique continuation kills it globally.

Let $q$ be a marked point with $dv(q)=0$. A simple closed curve in
an almost complex surface is embedded away from finitely many
critical and self-intersection points, by
[Wendl, Corollary 2.15.4](https://www2.mathematik.hu-berlin.de/~wendl/pub/jhol_bookv33_2.pdf).
The zeros of $u_1$ and the ramification points of $f$ are also finite.
Hence a sufficiently small punctured disk about $q$ consists of
globally injective scalar points, is unbranched over $C$, and misses
the core. This supplies an annulus for parameter extensions even
when $v$ has a cusp or $u_1(q)=0$.

Write $f(w)=w^e$, $e\ge2$, near $q$. For arbitrary $\alpha\in\mathbb C$
the localized variation $f_s=f+s\alpha\chi(w)w$ is holomorphic near
$q$. On the cutoff annulus choose
$j_s=(df_s)^{-1}j_Cdf_s$; outside it retain $j$. For small $s$ this
gives an actual family of finite holomorphic maps with
$\dot f_w(q)=\alpha$. Its induced vertical forcing is supported on
an inner disk.

A smooth real Cauchy--Riemann equation on a small disk admits a local
inhomogeneous solution with arbitrary value and holomorphic first
derivative at the center. Indeed the augmented operator

$$
 W^{2,p}\longrightarrow W^{1,p}\oplus\mathbb C^2\oplus\mathbb C^2,
 \qquad \xi\longmapsto(\bar\partial\xi,\xi(0),\partial_w\xi(0))
                                                                    \tag{8}
$$

has a bounded right inverse for $p>2$, obtained from a disk
$\bar\partial$ right inverse and a holomorphic affine correction.
Rescaling makes the smooth zeroth-order terms a small perturbation.
Solve the vertical linearized equation on a larger disk with
$\xi(q)=0$ and $\partial_w\xi_1(q)=\beta$, for arbitrary
$\beta\in\mathbb C$. Cut it off on a further outer annulus, where
the covering variation and its forcing are already zero.

The remaining residual $(r_1,r_2)$ is supported on that outer
injective annulus. Prescribe there
$\dot A=-r_1/\overline{f_w}$ and $\dot F=-r_2/\overline{f_w}$,
and extend from the scalar curve by a tubular neighborhood and cutoff.
Global scalar injectivity on this annulus prevents an effect on
other points of the same curve. Its image has positive distance
from the core, so these perturbations preserve (5).

The resulting universal tangent has arbitrary
$(\dot f_w(q),\dot u_{1,w}(q))=(\alpha,\beta)$, which supplies all
four real directions of the projected derivative. If the target
almost complex splitting is not the product splitting, these
coordinates differ from complex tangent coordinates by an invertible
real shear. The conclusion is unchanged. Passing to a slice for
reparametrizations transforms the mark at the same time; its Hessian
term cancels the map's gauge term at the zero derivative. This
retains the required direction at a simple branch point $e=2$.

Use high-$W^{k,p}$ charts, $k\ge3$, for moving-mark derivative
evaluation. First prove surjectivity of the ordinary universal
equation at $W^{1,p}\to L^p$, then use finite smooth parameter
directions and elliptic regularity to upgrade it. The four smooth
kernel directions just constructed give the augmented jet
surjectivity at high regularity. Derivative evaluation is never
treated as a bounded functional on $W^{1,p}$.

The path parameter and mark increase index zero to three. The
projected derivative constraint has real codimension four. Relative
parametric transversality therefore excludes its zero set along
generic paths. The fixed endpoint collars already have no noncore
graphs. Countably many degrees, genera and local charts are handled
in the stated smooth Baire space. For the conclusions below it is
enough to impose this through the fixed cutoff $M$.

## 4. The remaining actual-curve isolation conditions

**Lemma.** Such paths can be chosen so that every simple curve in
$V$ of degree at most $M$ is embedded; distinct such curves are
disjoint; each fixed-parameter simple moduli space is discrete; and
a dense set of parameters makes all these curves unobstructed.

**Proof.** The core is regular by (5). For a noncore graph, scalar
immersion from Section 3 gives the exact normal-operator sequence

$$
 0\longrightarrow(f^*L_2,f^*D_b)
 \longrightarrow(N_U,D_U^N)
 \longrightarrow(N_v,Q_v)\longrightarrow0.                         \tag{9}
$$

In particular $\ker D_U^N\hookrightarrow\ker Q_v$. This is precisely
the structural input in the restricted Petri and fold argument of
Calculation A2, Section 3. Its other inputs hold in arbitrary degree:
the full normal operator has index zero; scalar normal sections form
a real Cauchy--Riemann line equation; and there are unbranched scalar
injective disks off the core and off any finitely many other scalar
images.

For completeness, perturbation values on such disks give universal
normal regularity. Keeping those values zero on the curve, arbitrary
first normal derivatives give all real maps from the scalar normal
line to the full rank-two normal output. At nullity one this makes
the map to $\operatorname{Hom}(\ker D_U^N,\operatorname{coker}D_U^N)$
onto. At nullity $k\ge2$, two independent kernel sections project
to sections $q_1,q_2$ that are pointwise real independent on some
disk: otherwise their real ratio is holomorphic and constant,
contradicting unique continuation. The two-column annihilator
$q_1\otimes\eta_1+q_2\otimes\eta_2=0$ then vanishes. Thus the
derivative onto these two columns has rank $2k$, and the nullity-$k$
locus is locally contained in a submanifold of codimension at least
$2k\ge4$. A generic path avoids it.

At nullity one, second normal derivatives vary the quadratic
Kuranishi coefficient independently of values and first derivatives.
A nonzero scalar kernel section and an adjoint generator overlap
on a disk, and a quadratic real map on that scalar vector has
arbitrary output. Therefore the one-dimensional path moduli have
only Morse critical points of the parameter projection. Their
slices are discrete, and their critical values form a countable set.
Its complement is dense.

Two-point evaluation is transverse with the same restricted
parameters. At fixed marks use $W^{1,p}\to L^p$ augmented by point
values; parameter perturbations off the marks and adjoint unique
continuation kill the $L^q$ annihilator. Normal values and mark
motions then kill the evaluation covectors. For distinct curves
use separate scalar disks; for the fixed core use a disk of the
other curve off that core. Upgrade to high Sobolev regularity before
varying the marks, as in Calculation A2. The collision locus has expected
path dimension $1+4-6=-1$. Hence there are no self-intersections,
intersections of distinct curves, or intersections with the core.
Intersecting these residual conditions with Section 3 proves the
lemma. $\square$

All uses of the normal sequence occur after scalar immersion has
been established. The argument concerns actual simple curves;
it imposes no genericity on branched covers of the core or on
compactified covering spaces. In particular it does not assume
a simple-wall theorem inside the split operator locus.

## 5. Prime-degree cluster transport

Fix a prime $p$ and apply the preceding construction with $M=p$.
We use the cycle compactness and cluster results of
[Doan--Ionel--Walpuski](https://arxiv.org/html/2103.08221v2), abbreviated
DIW. Only their localization to curves whose supports lie in $V$
is needed. The argument is the one of Calculation A2, with the following
degree filtering made explicit.

At a fixed parameter, all simple local supports through the cutoff
are embedded and pairwise disjoint, and their spaces are discrete
and countable. Hence all cycles through this cutoff are countable:
they are finite sums with positive integral multiplicities. Choose
small closed Hausdorff balls about $C$ whose compact sets lie in
$V'\Subset V$. The cycle space in such a ball is compact by DIW
Theorem 2.11. The set of its distances from $C$ is a countable
compact set; it includes every multiple of $C$ at radius zero.
Choose arbitrarily small radii outside it. This selects compact
open-and-closed cycle insertions without ignoring disconnected
cycles or multiple supports.

Invertibility makes the core the unique degree-one support in a
sufficiently small ball. Restrict each boundary-free insertion to
DIW's connected-cycle space. It is relatively open and closed there,
so Proposition 3.10 transports its connected Gromov--Witten contribution
over a short interval. The stronger all-cycle boundary exclusion
also follows from the same properness argument. The parametrized
implicit function theorem preserves uniqueness in degree one after
shortening the interval.
Cover the path by finitely many intervals. At each overlap select
a slice from the dense set where all simple curves are regular,
and refine the two neighborhoods to a common smaller one, as in
DIW Proposition 3.19. The refinement differences are compact and
open and closed in the selected cycle spaces. They contain no
degree-one support and no multiple of such a support.

Now retain connected degree-$p$ cycles in a refinement difference.
Connectedness is closed under Hausdorff convergence of nonempty
compact sets, so this is a compact subset; degree is fixed by
the integral class. Disjointness of the embedded supports means
any such cycle is $rD$ for one support $D$. Its projection degree
satisfies $r\deg(D)=p$. Since degree one is excluded, primality
gives $r=1$ and $\deg(D)=p$. Thus this compact subset consists
of simple degree-$p$ supports. DIW Propositions 2.21 and 2.24
identify their cycle topology with the embedded-curve topology.
It is discrete and therefore finite. Every selected core is an
actual triangular curve, so its genus lies in $[H_p,G_p]$ by (6).

Each such core is primitive relative to the degree-$p$ coefficient.
Transport a sufficiently small insertion around the regular core
to a nearby globally generic ambient almost complex structure,
using DIW Propositions 3.20--3.21. The implicit function theorem
and openness of the embedded multiplicity-one locus preserve its
core genus and sign. The primitive case of DIW Proposition 3.54
then gives its complete degree-$p$ contribution

$$
                         \sigma_D(2\sin(t/2))^{2h(D)-2}.            \tag{10}
$$

This final deformation may leave the triangular class: the genus
of the selected core is already fixed and any proper multiple of
it has degree at least $2p$. Small insertions also exclude unrelated
cores of degree $p$. The same argument is available for the
degree-two and degree-three corrections simultaneously at cutoff
three, since their proper multiples have degree at least four.

Connected stable maps always have connected image support. They
therefore contribute only on the connected subset just considered;
no assumption that connectedness is an open condition is needed.
The compact clopen insertions are the objects transported, before
restricting to connected maps and the coefficient of $q^p$.
All stable maps over each selected support, including their
contracted components, are retained in (10). The same finite set
of supports works coefficientwise for every source genus.

Adding the refinement differences gives integers $c_{p,h}$ with

$$
 q^p
   =\sum_{h=H_p}^{G_p}c_{p,h}(2\sin(t/2))^{2h-2}.                   \tag{11}
$$

This is the localized simple-isotopy argument of DIW Proposition
3.30, with its assumptions verified above. No single isolating
tube for the whole path is asserted; the insertions vary along
the finite interval cover inside the common geometric region $V$.

The equal-spin endpoint has elementary series $G_g(q,t)$ and only
its degree-one elementary factor. In the logarithmic elementary
expansion at prime degree $p$, only factors of degrees one and $p$
contribute. The endpoint degree-one factors agree. The triangular
identification in Article III, Proposition 3.1 and Corollary 3.2 therefore gives
$c_{p,h}=n_{p,h}$, proving (1). This uses the fixed-target
identification of universal exponents; it introduces no numerical
evaluation on arbitrary analytic families. Composite degrees below
$p$ do not enter this coefficient.
