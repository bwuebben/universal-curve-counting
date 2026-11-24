# B3. The second genus excess in every covering degree

Let $C$ be a smooth curve of genus $g\ge2$, let $a,b$ be theta
characteristics, and use their positive real spin operators. Set
$c=g-1$ and $H_m=1+mc$. Write $L_{m,h}$ for the complete connected
local invariant, $c_m$ for the signed stack count of connected étale
degree-$m$ covers, and $F_m$ for the finite-map Euler integral in source
genus $H_m+1$, as in Article V, §§4–6.
Let $M_m$ be the reduced closure of the smooth connected degree-$m$
covers with four distinct simple branch points in
$\overline{\mathcal M}_{H_m+2}(C,m)$, and put

$$
                   I_m(a,b)=\int_{[M_m]}e(V_a)e(V_b).                 \tag{1}
$$

The real spin cokernels $V_a,V_b$ have rank four at this genus excess.
The notation $[M_m]$ means the sum of the fundamental cycles of its
four-dimensional irreducible components; irreducibility is not asserted.

**Theorem.** For every $m\ge2$,

$$
 \boxed{L_{m,H_m+2}
   =I_m(a,b)-\frac{mc}{12}F_m
          +\frac{mc(5mc-1)}{1440}c_m.}                              \tag{2}
$$

This extends the degree-two reduction of
[Calculation B1](second-excess-double-cover-components.md). The new
geometric issue is that higher degrees admit rational components joining
three sheets, ramified branches at direct nodes, and several nodes over
one target point. All these strata belong to $M_m$. After removing only
the summands with degree one over an étale target cover, the virtual
four-cycle is $[M_m]$ plus the elliptic-tail image of the first-excess
finite-map surface. Intermediate étale factorizations remain included.

For distinct even spins over a genus-two core, (2) gives the particularly
concrete remaining cubic target

$$
 L_{3,6}=I_3(a,b)+\frac{319}{180},\qquad
 n_{3,6}=I_3(a,b)-\frac{43}{2}.                                    \tag{3}
$$

## 1. Separating the full étale-lift summands

Let $\mathcal X_m=\overline{\mathcal M}_{H_m+2}(C,m)$. For every
connected étale degree-$m$ cover $p:Y\to C$, its full degree-one
lift summand is

$$
 \left[\overline{\mathcal M}_{H_m+2}(Y,1)
                         /\operatorname{Aut}_C(p)\right].           \tag{4}
$$

The image of the source fundamental group is locally constant in
stable-map families. Thus their union $\mathcal E_m$ and its complement
$\mathcal X_m^\circ$ are proper open-and-closed substacks. This is
the decomposition by étale lifts of
[Bryan--Pandharipande, Definition 2.2 and Proposition 2.3](https://arxiv.org/html/math/0009025v3#S2.SS2).
For composite $m$, $\mathcal X_m^\circ$ still contains maps factoring
through étale covers of smaller degree. It must not be called the
primitive summand. At prime $m$ it is that summand.

A degree-one map to the smooth curve $Y$ has a unique nonconstant
component, mapping isomorphically, and every connected contracted block
attaches to it just once. Distinct points cannot be identified over
the same target point of $Y$. This characterizes (4) among maps whose
nonconstant normalization is an étale degree-$m$ cover of $C$.

## 2. All second-excess strata

Normalize the nonconstant components of a connected stable map and
write their disjoint union as $Y=\coprod_{j=1}^sY_j$. Let $2r$ be
their total ramification and let $t$ count direct nodes between their
branches. For each maximal connected contracted subcurve let $q_i$
be its arithmetic genus and $a_i$ its number of attachments to $Y$.
Riemann--Hurwitz and the connected dual-graph formula give

$$
 \sum_j g(Y_j)=mc+s+r,\qquad
 \boxed{r+t+\sum_i(q_i+a_i-1)=2.}                                  \tag{5}
$$

Stability gives $2q_i-2+a_i>0$. Every term in (5) is nonnegative,
and a stable contracted block has positive cost. Therefore the only
possible block types are

$$
                   (q,a)=(1,1),(2,1),(1,2),(0,3).                  \tag{6}
$$

In particular three is the attachment bound at this excess in every
covering degree. The $(0,3)$ block is a three-pointed rational curve.
Positive-genus blocks include their full stable nodal boundaries.

The list of possibilities and largest stratum dimensions is as follows.
The normalization has as many connected components as its degree
partition and the indicated connected graph permit; no restriction to
irreducible normalization is made.

| $r,t$ | Contracted blocks | Largest dimension |
|---|---|---|
| $2,0$ | none | $4$ |
| $1,1$ | none | $3$ |
| $0,2$ | none | $2$ |
| $1,0$ | one elliptic tail $(1,1)$ | $4$ |
| $0,1$ | one elliptic tail $(1,1)$ | $3$ |
| $0,0$ | two elliptic tails | $4$ |
| $0,0$ | one genus-two tail $(2,1)$ | $5$ |
| $0,0$ | one elliptic bridge $(1,2)$ | $3$ |
| $0,0$ | one rational tripod $(0,3)$ | $1$ |

These dimensions count branch positions, node or attachment images,
and $\dim\overline{\mathcal M}_{q,a}=3q-3+a$. Coincident images and
higher ramification lower them. A finite cover with prescribed branch
data has only finitely many monodromy possibilities in fixed degree.

The two-tail and genus-two-tail rows necessarily have connected étale
nonconstant normalization, since singly attached blocks cannot connect
different components. They belong to $\mathcal E_m$. Their five-dimensional
excess stratum is consequently removed before using a top-dimensional
fundamental-cycle argument. All other rows lie in $\mathcal X_m^\circ$:
a direct node or a multiply attached contracted block is incompatible
with a degree-one lift, while ramification is incompatible with an
étale normalization.

## 3. The finite, bridge and tripod loci lie in the main closure

We use admissible covers only to produce smoothings; no isomorphism
between their compactification and the stable-map compactification is
assumed. Expand the target at a point $x$ to $C\cup_x\mathbb P^1$.
Attach covering components over the bubble, with pole orders equal to
the local degrees of the branches on $Y$. Stabilizing the composite
map to $C$ then recovers the prescribed contracted component or node.

The balanced equations at a target node are
$uv=\tau$, $xy=\tau^e$, $x=u^e$, $y=v^e$. The deformation and
admissible-cover comparison theorems of
[Abramovich--Corti--Vistoli, Theorem 3.0.1, Corollary 3.0.4 and Proposition 4.2.2](https://arxiv.org/html/math/0106211v1)
allow smoothing after finite base change. Choose the branch markings
on the fixed $C$ to collide in the indicated bubbles. The stabilized
unmarked target remains exactly $C$ throughout. Hence these are
fixed-target smoothings. Connectedness of the constructed admissible
source ensures connected nearby smooth sources.

For a rational tripod choose the three simple poles $0,1,\infty$ of

$$
                    \phi(z)=z+z^{-1}+(z-1)^{-1}.                    \tag{7}
$$

It has four distinct simple finite branch values. Its fiber polynomial
has discriminant $w^4-2w^3-5w^2+6w-23$, whose discriminant is
$-5619712$. Glue its poles to the three attachments; unused sheets
receive degree-one bubble components. Stabilization retains exactly
the given tripod. This works equally when its attachments lie on
one, two or three nonconstant components.

For a smooth elliptic bridge $(E,p,q)$, the pencil $|p+q|$ gives a
degree-two map to the bubble, with simple poles at the attachments
and four simple finite branch points. Add the unused degree-one
sheets. Stabilization recovers the bridge. Limits of smooth pointed
elliptic curves include every nodal bridge, so the complete bridge
stratum lies in $M_m$.

For finite nodal maps, a node in the $r=1,t=1$ row can join local
degrees only of types

$$
                         (1,1),(2,1),(3,1),(2,2).                  \tag{8}
$$

Their rational bubble maps can be taken as, respectively,

$$
 z+z^{-1},\quad z^2+z^{-1},\quad z^3+z^{-1},\quad z^2+z+z^{-2}.
                                                                    \tag{9}
$$

Each has two poles of the indicated orders and distinct simple
finite branch values. The twice-attached rational component contracts
to the required node after stabilization. A pole profile $(e_1,e_2)$
uses $e_1+e_2$ finite bubble ramification units and replaces
$e_1+e_2-2$ units on the old branches, thereby adding exactly two.

For $m\ge4$, two direct nodes may share a target image. Give each
pair its own rational covering component over the same target bubble,
and choose their branch sets disjoint by affine target changes.
Any unused sheet which is itself ramified at that image needs a
polynomial bubble cover of the corresponding degree, rather than a
degree-one component. Its single pole attaches with the required local
degree; its finite simple branch points replace the old ramification
without changing the total. Those one-attachment components disappear
on stabilization. All bubble branch sets can be chosen disjoint.

This treats the $r=1,t=1$ and $r=0,t=2$ rows. Higher ramification
on smooth finite covers also deforms to distinct simple branching,
since the finite-map deformation sheaf is torsion with no first
cohomology. Thus every finite map in (5), every elliptic bridge and
every rational tripod belongs to $M_m$.

## 4. The virtual four-cycle

Let $\mathcal F_m$ be the proper first-excess finite-map surface from
Article V, §§4–6 and let $\mathcal U_m\to\mathcal F_m$ be its universal
one-pointed stable curve. A marking at a finite node is carried by
an inserted rational component. Gluing a pointed elliptic curve gives
a proper map

$$
 \gamma:\mathcal U_m\times\overline{\mathcal M}_{1,1}
                                      \longrightarrow\mathcal X_m^\circ.
                                                                    \tag{10}
$$

Its image $G_m$ contains both elliptic-tail rows of the table and
their limits, including marked-node configurations. Conversely those
rows are obtained in this way. The generic stack degree onto each
image component is one: the finite core and its attached elliptic
tail are intrinsic. No uniform factor $1/m$ is inserted.

Sections 2–3 show that the support of $\mathcal X_m^\circ$ is
$M_m\cup G_m$. Both are pure four-dimensional closures. At a smooth
four-simple-branch cover, finite-map deformation theory is unobstructed.
At a generic elliptic tail, its attaching point on the two-branch
finite core is unramified. The node-smoothing line maps isomorphically
to $T_{f(x)}C\otimes H^1(E,\mathcal O_E)$ by the differential of
the core map and the elliptic principal-part map. Thus this locus
is also unobstructed of dimension four. Both generic virtual
multiplicities are one.

A dimension-four virtual cycle on this pure four-dimensional stack
is determined by its generic multiplicities. Hence

$$
 [\mathcal X_m^\circ]^{\rm vir}
  =[M_m]+\gamma_*[\mathcal U_m\times\overline{\mathcal M}_{1,1}].      \tag{11}
$$

The identity does not require global smoothness or disjointness of
the displayed components. It does not apply the same argument to
the removed five-dimensional étale-lift stratum.

## 5. Normal orientations and the two explicit terms

The Lee--Parker spin bundles and their determinant orientations are
the same as in Article V, §§4–6 and Calculation B1. On the source of (10), the
normalization sequence cancels the constant tail sections against
the matching fiber. The remaining tail factor is
$\mathbb E^\vee\otimes\operatorname{ev}^*(a\oplus b)$ with complex
orientation. This continues over a marked finite node: the inserted
rational component has no first cohomology and its matching sequence
recovers the original finite core.

With $\lambda=c_1(\mathbb E)$ and $\kappa=c_1(K_C)$,

$$
 c_2(\mathbb E^\vee\otimes\operatorname{ev}^*(a\oplus b))
                    =-\lambda\operatorname{ev}^*\kappa.
$$

The squared Hodge class vanishes on $\overline{\mathcal M}_{1,1}$;
the product of two classes from $C$ vanishes as well. Since
$\int\lambda=1/24$ and $\pi_*\operatorname{ev}^*\kappa=2mc$,
the elliptic-tail contribution is $-mcF_m/12$. The stack integrals
already retain all stabilizers, including the pointed elliptic one.

On (4), the differential of the étale cover identifies tangential
obstruction theories. Each contracted block attaches once, so the
holomorphic kernels with coefficients in $p^*a,p^*b$ are constant
over the proper degree-one moduli. The constant-rank quotient and
orientation formula in
[Lee--Parker, §12](https://users.math.msu.edu/users/parker/LeeParkerObstructionBundle.pdf)
and the degree-one calculation of
[Pandharipande, Theorem 2](https://arxiv.org/abs/math/9811140)
give the complete signed contribution $\sigma(p)s(t)^{2mc}$, where
$s(t)=2\sin(t/2)$. Summing with the deck weights gives $c_m s(t)^{2mc}$.
This includes the excess-dimensional genus-two-tail locus with its
virtual contribution. Finally

$$
 [t^{2mc+4}]s(t)^{2mc}=\frac{mc(5mc-1)}{1440},                       \tag{12}
$$

which proves (2).

## 6. Cubic calibration

In base genus two and degree three, the established first-excess and
étale counts are

$$
 \begin{array}{c|cc}
 & a\ne b\text{ even}& a=b\text{ even}\\\hline
 c_3&28/3&220/3\\
 F_3&-6&-70.
 \end{array}                                                       \tag{13}
$$

Thus the mixed tail and étale terms are $3/2$ and $49/180$,
respectively, proving the first identity in (3). The known elementary
extraction in [Calculation A3](../support/prime-degree-genus-cutoffs.md) gives
$L_{3,6}=4189/180+n_{3,6}$, proving the second. For equal spins,
the complete elementary coefficient is $1673/36$; subtracting its
tail contribution $35/2$ and étale contribution $77/36$ gives

$$
                            I_3(\theta,\theta)=\frac{161}{6}.       \tag{14}
$$

The common complex-index class is $c_2(-R\pi_*f^*\theta)$,
independent of the theta characteristic. Define only as cohomology
classes

$$
 \Delta_\theta=\tfrac12\bigl(c_2(-R\pi_*f^*\theta)-e(V_\theta)\bigr).
$$

They are supported on the holomorphic-kernel loci. No identification
with reduced kernel surfaces is implicit in this definition. The
ordinary identity comparing the square of a difference with the two
squares, together with (14), gives

$$
 \boxed{n_{3,6}=\frac{16}{3}
       -2\int_{[M_3]}(\Delta_a-\Delta_b)^2.}                        \tag{15}
$$

The corresponding defect pairing on the separate elliptic-tail
component is already determined: it is $8$. Indeed, on the source
of (10), the difference of the two Euler classes is the difference
of the core Euler classes times the common tail first Chern class
$\operatorname{ev}^*c_1(\theta)-\lambda$. Squaring and pushing forward
gives

$$
 \int_{[G_3]}(\Delta_a-\Delta_b)^2
   =-\frac{3}{48}\bigl((-70)+(-70)-2(-6)\bigr)=8.                  \tag{16}
$$

This explains the equivalent primitive-summand expression
$n_{3,6}=64/3-2\langle(\Delta_a-\Delta_b)^2,
[\mathcal X_3^\circ]^{\rm vir}\rangle$ on the primitive summand.
Article V, Corollary 1.2 gives the exponent independently. Taking
$n_{3,6}=-20$ would require that pairing to be $38/3$; this is a
consistency check.

The distinction between real Euler classes and complex index classes
is essential. At virtual complex rank two, higher Chern classes need
not vanish; even the stable Pontryagin identity involves
$c_2^2-2c_1c_3+2c_4$. Smooth-cover disjointness of two kernel loci
therefore supplies neither their compact self-intersections nor
their intersections on contracted boundaries. Calculations C1–C4 and Article V, §7 evaluate these intersections.

The exact check
`code/second_excess_all_degrees.py`
uses Sylvester determinants to verify every rational bubble's branch
polynomial and its squarefreeness, and reproduces the cubic calibration
from the established hook-count code. It also checks the known tail
defect pairing. The mixed integral is evaluated by Theorem 3.1 together with this reduction.
