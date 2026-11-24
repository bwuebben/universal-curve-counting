# C4. Ramification derivatives on a cubic Hurwitz divisor

## 1. Statement and dependencies

Let

$$
 C:\quad y^2=\prod_{i=1}^6(x-i),\qquad
 \theta_T=\mathcal O_C\left(\sum_{i\in T}W_i\right)K_C^{-1},
 \quad |T|=3,
                                                               \tag{1}
$$

where $W_i$ is the Weierstrass point over $i$. Complementary triples
give the same line, and these are the ten noneffective even theta
characteristics. For one such line write $\theta=\theta_T$.
Let $\mathcal F$ be the proper finite genus-five cubic Hurwitz
surface and $D_\theta$ its determinantal kernel divisor from
[Calculation C3](single-spin-cubic-kernels-and-spectral-incidence.md).
Let $R\to\mathcal F$ be the degree-two ramification scheme.

The following is the additional conclusion of this calculation.

**Proposition.** For each of the ten lines in (1), only finitely many
smooth connected genus-five cubic kernel covers have a ramification
point $r$ at which the unique spin section $s$ satisfies $ds(r)=0$.
Consequently $ds(r)\ne0$ at the generic point of every reduced
component of $R\times_{\mathcal F}D_\theta$.

At a ramification point the derivative is independent of the local
frame of $\theta$, since a change of frame changes it by a multiple
of $df(r)=0$. Multiplying $s$ by a nonzero scalar also preserves its
nonvanishing.

The component calculation in §7 below proves that, when $|\theta^3|$ has no total ramification, the
divisor $D_\theta$ has no nodal-boundary component, every compact
elliptic-bridge restriction of $\Delta_\theta$ vanishes, and the
only remaining kernel surfaces are the finite-source surface and
ramification-tail surfaces. It also proves that the generic derivative
condition in the proposition identifies the tail cycle with its
determinantal test family, including multiplicities. Those are inputs
to the corollary

$$
 S_\theta=J_{\mathrm{fin},\theta}+1,
 \qquad
 J_{\mathrm{fin},\theta}
 =\left\langle\Delta_\theta,[Z_{\mathrm{fin},\theta}]\right\rangle.
                                                               \tag{2}
$$

Here $Z_{\mathrm{fin},\theta}$ is the closure of the smooth-source
genus-six kernel surface, and $S_\theta=\int_M\Delta_\theta^2$
on the main genus-six four-cycle. Article V, Theorem 7.3 computes the complete square $S_\theta=6$ independently, so (2) gives $J_{\mathrm{fin},\theta}=5$ on this curve. Calculation C3 obtains the same number on an auxiliary spectral incidence. A direct comparison of those compactifications at singular spectra remains open.

## 2. A meromorphic frame for the spectral equation

Put

$$
 P=\prod_{i\in T}(x-i),\quad Q=\prod_{i\notin T}(x-i),
 \quad \beta=\frac yP,\quad \beta^2=\frac QP.
                                                               \tag{3}
$$

If $D_P=\sum_{i\in T}W_i$ and $K_\infty$ is the sum of the
two points at infinity, then $\theta=\mathcal O(D_P-K_\infty)$.
A meromorphic frame $t$ with this divisor can be normalized so that
$t^2=P\,dx/y$. The two sections of $\theta^3$ are
$t^3/P$ and $\beta t^3/P$. In the fiber coordinate $z=s/t$,
the trace-free spectral equation is therefore

$$
 G(x,z)=Pz^3+(a+bx)z+c+d\beta=0.                              \tag{4}
$$

It is written on the complement of the six Weierstrass points and
infinity; regular frames will be used at those points. Its parameter
scaling is $(a,b,c,d)\mapsto(\lambda^2a,\lambda^2b,
\lambda^3c,\lambda^3d)$ with $z\mapsto\lambda z$.
The proper singular-point incidence is the surface $I$ of Calculation C3.

The pencil $|\theta^3|$ is basepoint-free because the zero divisors
of its two displayed sections are $D_P,D_Q$. It has no total
ramification in these examples. Indeed the polynomial

$$
                       Q'P-QP'                               \tag{5}
$$

has degree four, is squarefree and is coprime to $PQ$ for every
one of the ten partitions. The map $\beta$ is unramified at the
Weierstrass points, where it has a simple zero or pole. It is
unramified at infinity since the sums of the roots of $P$ and $Q$
are distinct. Thus (5) gives eight simple ramification points on
$C$. For $T=123$ it equals

$$
 9x^4-126x^3+621x^2-1260x+876,
 \qquad \operatorname{disc}(Q'P-QP')=2317870080.
                                                               \tag{6}
$$

## 3. Irreducibility of the nonnodal singular incidence

At a singular point $(p,u)$ with $u\ne0$, the three equations
$G=G_z=G_x=0$ give

$$
 a+bx=-3Pu^2,\qquad
 b=-P'u^2-\frac{d\beta'}u,\qquad
 c=2Pu^3-d\beta.
                                                               \tag{7}
$$

On this open set use the scaling to set $u=1$. Define

$$
 r=\frac{\beta'}\beta
   =\frac12\left(\frac{Q'}Q-\frac{P'}P\right),\qquad
 A=6PP''-4(P')^2,\qquad
 B=6P(r'+r^2)+4P'r,
 \qquad e=d\beta(p).
                                                               \tag{8}
$$

The Hessian determinant of (4) vanishes precisely when

$$
                         A+Be-r^2e^2=0.                       \tag{9}
$$

For example its three entries before setting $u=1$ are
$G_{zz}=6Pu$, $G_{xz}=2P'u^2-d\beta'/u$, and
$G_{xx}=P''u^3+d\beta''$. These give (9) directly.
The discriminant of the quadratic (9) is $B^2+4Ar^2$.

For $T=123$,

$$
 A=-4(3x^2-12x+13),\qquad
 B^2+4Ar^2=\frac{27N(x)}{4Q^4P^2},                            \tag{10}
$$

where

$$
\begin{aligned}
N={}&99x^{16}-5544x^{15}+142164x^{14}-2212056x^{13}
 +23324514x^{12}-176206968x^{11}\\
 &+982686852x^{10}-4104078888x^9+12867095475x^8
 -29998678752x^7+50704505488x^6\\
 &-58818910080x^5+40881679264x^4-8703496704x^3
 -9613590272x^2+7238762496x-1211897088.
\end{aligned}                                                   \tag{11}
$$

The polynomial $N$ is squarefree and coprime to $PQ$. More generally,
for all ten partitions the reduced numerator of $B^2+4Ar^2$ has
degree sixteen, is squarefree and is coprime to $PQ$, and its
denominator is a square up to a nonzero constant. These are exact
polynomial checks using (8), with the prescribed integer roots.
Each numerator root has two unramified preimages on $C$, at both
of which the discriminant has odd valuation. It is therefore not
a square in $\mathbb C(C)$. Hence (9) is irreducible over that
function field and gives one irreducible curve dominating $C$.

There is no additional curve over an excluded point of $C$. At a
simple ramification point of $\beta$, one has $r=0,r'\ne0$,
so (9) is affine in $e$ with nonzero coefficient $6Pr'$. At an
ordinary point with $r\ne0$ its quadratic coefficient is nonzero.
At infinity the same argument applies in the coordinate $1/x$,
with the corresponding regular frames; $\beta$ is unramified
there. At a Weierstrass point use the local square-root coordinate
and a regular theta frame. The canonical sections have even local
expansions, while one of the two cubic-pencil sections has a simple
zero. Solving the singularity equations leaves one canonical
coefficient free. In the Hessian its coefficient is a nonzero
multiple of $u^2$, so again the Hessian is not identically zero
on the incidence fiber.

Finally consider $u=0$. The incidence equations require
$A_{\rm sp}(p)=B_{\rm sp}(p)=B_{\rm sp}'(p)=0$, where the
subscript distinguishes the spectral coefficients from (8).
If $B_{\rm sp}\ne0$, the point is one of the eight simple
ramification points of the cubic pencil. It is not Weierstrass,
and the Hessian is $-(A_{\rm sp}'(p))^2$. It vanishes only for
$A_{\rm sp}=0$, giving finitely many points. If $B_{\rm sp}=0$,
the Hessian can vanish only when $A_{\rm sp}$ has a double zero,
at one of the six Weierstrass points. These too are isolated
incidence points. Thus the nonnodal singular incidence consists of
the single dominating curve, together with at most finitely many
additional points.

## 4. A general member has exactly one ordinary cusp

We exhibit a point of the curve (9) whose spectrum has one ordinary
cusp and no other singularity. Fix either point $p$ over $x=0$,
normalize $u=1$, and let $e$ satisfy (9) at $x=0$. Put

$$
 a=-3P(0),\quad b=-P'(0)-e r(0),\quad
 c=2P(0)-e,\quad d^2=e^2\frac{P(0)}{Q(0)}.                    \tag{12}
$$

Choose $d=e/\beta(p)$. The primitive integral quadratic for $e$
and the exact checks for the ten partitions are as follows.

| $T$ | Quadratic $H_T(e)$ |
| --- | --- |
| 123 | $5329e^2+125436e+748800$ |
| 124 | $189e^2+7568e+76800$ |
| 125 | $1083e^2+66700e+940800$ |
| 126 | $2809e^2+245112e+4377600$ |
| 134 | $1849e^2+184632e+4204800$ |
| 135 | $1369e^2+240150e+7142400$ |
| 136 | $121e^2+33492e+1209600$ |
| 145 | $243e^2+109400e+4627200$ |
| 146 | $529e^2+415344e+20966400$ |
| 156 | $289e^2+552300e+34617600$ |

Each displayed quadratic is irreducible over $\mathbb Q$. The
cubic discriminant, regarded as a section of $K_C^3$, has coefficient
$U+54cd\,y$ in the frame $(dx/y)^3$, where

$$
 U=4(a+bx)^3+27(Pc^2+Qd^2),\qquad
 \mathcal R=U^2-2916PQc^2d^2.                                 \tag{13}
$$

In the field $\mathbb Q[e]/(H_T)$, substitution of (12) gives
$\mathcal R=x^3J_T(x)$ with

$$
 \deg J_T=3,\quad J_T(0)\ne0,\quad
 \gcd(J_T,J_T')=\gcd(J_T,PQ)=1,\quad cd\ne0.                 \tag{14}
$$

These statements follow by polynomial reduction modulo the displayed
quadratics. They supply exact data for a finite verification, without
an approximate root calculation. The nonzero degree-six leading
coefficient of $\mathcal R$ excludes zeros of the discriminant at
infinity; coprimality with $PQ$ excludes them at Weierstrass points.
Since $cd\ne0$, over each root of $\mathcal R$ exactly one of
the two points of $C$ is a discriminant zero.

At $p$ the cubic has a double root $u=1$ and a distinct third
root. Its local double-root factor has discriminant of order three
in the local coordinate $x$. Completing the square gives an
ordinary cusp. The three remaining discriminant zeros are simple,
so they give smooth simple ramification points of the spectrum,
not additional singularities. Also $d\ne0$, so the spectrum is
integral by Calculation C3's reducibility criterion. Its normalization
has genus six.

Having precisely one ordinary cusp and otherwise being smooth is
open in the nonnodal incidence near this member. The single curve
of §3 is therefore generically of this type. Its locus of spectra
whose total delta invariant is at least two is finite.

## 5. Ramification derivatives and the tail cycle

Every component of $D_\theta$ is generically a smooth-source
cover, by the no-total-ramification statement of §2 and the
component reduction. Its section gives a finite birational map to
an integral spectral curve of arithmetic genus seven. For a smooth
source of genus five this is its normalization, and the total delta
invariant of the spectrum is two.

If $df(r)=ds(r)=0$, the normalization map to the smooth surface
$\operatorname{Tot}\theta$ is not immersive at $r$. Its image
singularity is not an ordinary node, so the marked spectrum belongs
to the nonnodal singular incidence. Sections 3–4 show that only
finitely many such spectra have delta invariant two. Each has a
unique normalization, and the section is unique up to scale. Thus
only finitely many smooth kernel covers can fail the derivative
condition. The finite ramification marking preserves this finiteness.
It cannot contain a component of $R\times_{\mathcal F}D_\theta$.
This proves the proposition.

For completeness, the imported local tail comparison uses a node
smoothing parameter $\tau$ and the pulled-back determinantal
equation $g_R$. Its two kernel equations reduce generically to
$(g_R,\tau h)$, where the restriction of $h$ is the nonzero
derivative $ds(r)$ times the elliptic principal-part isomorphism.
Thus the ideal is $(g_R,\tau)$ and retains every power in $g_R$.
The supported oriented multiplicity is positive, and gluing has
generic stack degree one. The actual tail cycle is consequently
$j_*[(R\times_{\mathcal F}D_\theta)\times
\overline{\mathcal M}_{1,1}]$. Its pairing with
$\Delta_\theta$ is $1$ by Calculation C3. Together with the compact
bridge vanishing, this proves (2) using the component calculation in §7.

## 6. Consequence for stable replacement and its remaining gap

On these explicit spin curves the residual spectral family of
Calculation C3 is nodal away from finitely many points of $I$. A spectrum
with a unique ordinary node or cusp gives a smooth partial
normalization after blowing up its marked singular point. Any
nonnodal residual fiber therefore requires either a more singular
marked point or an additional nonnodal unmarked point. Such a
spectrum has total delta at least two and lies in the finite set
of §§3–4. Reducible spectra cause no exception: two simple zeros
of their canonical coefficient give ordinary nodes, and a double
zero at a Weierstrass point gives the retained nodal rational tripod.

An integral residual fiber has arithmetic genus six and normalization
genus at least four. Its total delta is at most two. If that delta
is two, the normalization is an étale cubic, so all its graph
branches are immersed. Thus the only possible nonnodal singularity
in that case is an ordinary tacnode. If the delta is one, the only
nonnodal possibility is an ordinary cusp. This reduces the remaining
stable replacement problem to cusp and tacnode fibers over finitely
many incidence points.

This finiteness does not establish equality of compact characteristic
degrees. Stable reduction changes cohomology lattices on exceptional
divisors, and tacnode replacement can itself acquire a cuspidal
elliptic component. A comparison must control the actual index
sheaf and any successive modifications, including their intersections.
Article V, §7 instead obtains the smooth-source degree 5 by subtracting the ramification-tail degree from the individual square.

## 7. Components of the kernel cycle and the elliptic-tail multiplicity

Write $M$ for the reduced main four-cycle and $W=-R\pi_*f^*\theta$. The constant-kernel formula of [Calculation C3, §6](single-spin-cubic-kernels-and-spectral-incidence.md) gives
$$
j^*\Delta_\theta=c_2(j^*W)-\xi c_1(j^*W)+\xi^2,\qquad \xi=-c_1(\ker\bar\partial).
$$
We use the pencil condition $\theta^3\ne\mathcal O_C(3p)$ for every $p$, verified for the ten lines in §2.

**Unramified normalizations.** On a cyclic cubic the two trace-free summands of the space of spin sections are the sections of $\theta\eta^{\pm1}$ for a nontrivial line $\eta$ of order three. They have equal dimensions by Serre duality. Any nonzero section gives $\theta\eta=\mathcal O_C(p)$, contrary to the pencil condition. On a non-Galois unramified cubic let $\pi:D\to C$ be its discriminant double. The trace-free summand of the direct image is $\pi_*L$ for an order-three line $L$ on $D$, by induction of a nontrivial character of the cyclic subgroup in $S_3$. The degree-two bundle $\pi^*\theta\otimes L$ on $D$ has at most one section. Indeed, two sections would give a hyperelliptic pencil, whose square on this genus-three curve is $K_D$. Squaring the displayed line would give $L^2=1$, a contradiction. In partition $2+1$ the double contributes at most one section of $\theta\eta$ for a two-torsion line $\eta$ and the identity contributes none. In partition $1+1+1$ there are no sections.

**Matching on a finite nodal source.** On a connected smooth cubic a nonzero section is primitive over the function field, since the covering degree is prime and $\theta$ has no sections on $C$. Its three values are generically distinct. On a $2+1$ normalization they are $s,-s,0$. Thus a fixed kernel section matches two distinct sheets over only finitely many points. In genus five a finite nodal source has one node and unramified normalization. There are finitely many such normalizations, each with at most one projective section and finitely many matching pairs. No component of $D_\theta$ lies in the nodal boundary. In genus six a finite nodal source has one node and half-ramification one, or two nodes and unramified normalization. In the former case, the connected smooth genus-five cores form at most a divisor in the Hurwitz surface by C3. On a disconnected $2+1$ normalization a branched double with trace-free summand $A^{-1}$, $\deg A=1$, can have a section only for $A=\theta$. Its branch section varies in $\mathbb P H^0(K_C)$ of dimension one. Finite matching therefore bounds the finite nodal kernel locus by dimension one.

**Contracted blocks.** If $r$ is half the total ramification of the nonconstant normalization, $t$ the number of direct nodes, and $q_i,a_i$ the genus and number of external attachments of each connected contracted block, the genus budget is
$$
r+t+\sum_i(q_i+a_i-1)=2 .
$$
A rational tripod has three attachments and an unramified normalization. The possible kernel maps are finite by matching and $\overline{\mathcal M}_{0,3}$ is a point. A genus-two contracted block with one attachment, or two one-pointed elliptic tails, leaves a connected unramified cubic. Its full unramified summand has already been removed from the main cycle, as in [Calculation C2](cubic-common-kernel-excess-intersection.md). A single elliptic tail with one attachment lies over the ramification scheme of a genus-five core. Its surface components have generically smooth core by the preceding paragraph.

For an elliptic bridge, $r=t=0$ and the two sheets have equal section values at the attachments. There are finitely many choices of normalization and attachment pair. Order the attachments and rigidify that normalization. The compact bridge family is then $\overline{\mathcal M}_{1,2}$, with fixed image $q\in C$ and constant kernel line. Normalization at the two external nodes gives
$$
j^*W=W_{\rm norm}+\theta_q\otimes(\mathcal O+\mathbb E^\vee).
$$
Here $W_{\rm norm}$ is a constant virtual vector space of rank zero and $\mathbb E$ is the elliptic Hodge line. Thus $\xi=0$ and $c_2(j^*W)=0$, even when the two markings lie on a rational bubble. Finite push-pull descends the vanishing of $j^*\Delta_\theta$ to the unordered family. Its contribution to the square vanishes regardless of the multiplicity of this bridge surface in the supported cycle.

**Tail multiplicity.** At a generic ramification-marked genus-five core, let $g_R$ be the pulled-back local determinantal equation of $D_\theta$, and $\tau$ the smoothing of the elliptic attachment. The one-tail chart of C2 gives kernel equations restricting at $\tau=0$ to $(g_R,0)$. The second equation is $\tau h$ with $h|_{\tau=g_R=0}=ds(r)\,\delta$, where $\delta:T_eE\to H^1(E,\mathcal O_E)$ is the principal-part isomorphism at the elliptic attaching point. The derivative proposition above makes $h$ a unit generically on every reduced component. The kernel ideal is therefore $(\tau,g_R)$, including any multiplicity in $g_R$. The complex determinantal orientation has positive sign by C2 and C3. The smooth tail map recovers the core, marked ramification point and elliptic tail, so the generic gluing degree as stacks is one. The supported tail cycle thus has the multiplicities of the compact test family of C3, §7, whose pairing is $1$. Only the smooth-source closure remains. This proves (2).
