# C2. The common spin-kernel intersection on the cubic main locus

## 1. Statement and notation

Let $C$ be a smooth genus-two curve and let $a,b$ be distinct even
theta characteristics. They are noneffective. Write $\mathcal F$ for
the proper surface of connected finite cubic maps of arithmetic genus
five, with its reduced nodal boundary, and let $M=M_3$ be the reduced
closure of smooth connected cubic covers with four distinct simple
branch points in $\overline{\mathcal M}_6(C,3)$. Thus $M$ is the
proper four-cycle of [Calculation B3](../excess/second-excess-reduction-in-all-degrees.md),
with all of its components retained. It lies in the primitive summand;
the full degree-one-upstairs étale summands are separate.

For the universal map $f$ and curve $\pi$ over $M$, put

$$
 W_\theta=-R\pi_*f^*\theta,\qquad
 c=c_2(W_\theta),\qquad
 \Delta_\theta=\frac{c-e(V_\theta)}2\in H^4(M;\mathbb Q),
 \qquad\theta=a,b.                                               \tag{1}
$$

Here $V_\theta$ is the oriented real rank-four positive-spin-mass
cokernel bundle. Rational Grothendieck--Riemann--Roch makes $c$
independent of the theta characteristic. The virtual complex bundle
$W_\theta$ has rank two; this does not assert that its higher Chern
classes vanish. The classes in (1) are defined before making any
identification with holomorphic kernel cycles.

Let $S_\theta$ be the residue-one Szegő kernel and set

$$
 Q=\left.(S_a^2-S_b^2)\right|_\Delta\in H^0(C,K_C^2).
                                                                  \tag{2}
$$

**Theorem.** The determinantal kernel loci $K_a,K_b$ on $\mathcal F$
have scheme-theoretic intersection isomorphic to the zero scheme of
$Q$. Its length is four. Every common finite map is smooth and is
one of the covers $f_p:Y_p\to C$ of
[Calculation C1](cubic-covers-with-two-spin-sections.md).
On $M$, the common holomorphic-kernel support consists exactly of the
disjoint closed curves

$$
 E_p=\{f_p:Y_p\to C\text{ with a contracted pointed genus-one tail
 at its unique triple ramification point }r_p\}
 \simeq\overline{\mathcal M}_{1,1},
 \qquad p\in Z(Q)_{\mathrm{red}}.                               \tag{3}
$$

If $Q$ has four simple zeros, then

$$
                  \boxed{\displaystyle
                  \int_{[M]}\Delta_a\Delta_b=-\frac13.}         \tag{4}
$$

Each of the four curves in (3) contributes $-1/12$. This includes
the nodal elliptic endpoint and its stabilizer. The numerical assertion
is stated under the simple-zero hypothesis; no specialization argument
at colliding zeros is needed here.

Define $S_a^{\mathrm{self}}=\int_{[M]}\Delta_a^2$ and
$S_b^{\mathrm{self}}=\int_{[M]}\Delta_b^2$. The equal-spin calibration
and component reduction of Calculation B3 give the following consequence:

$$
 \begin{aligned}
 I_3(a,b)&=\frac{51}{2}-2(S_a^{\mathrm{self}}+S_b^{\mathrm{self}}),\\
 n_{3,6}&=4-2(S_a^{\mathrm{self}}+S_b^{\mathrm{self}}).
 \end{aligned}                                                   \tag{5}
$$

## 2. The finite-cover intersection as a scheme

This section records the new local argument; the construction of the
individual covers remains that of Calculation C1. Near any common cover,
the surface $\mathcal F$ is smooth. For $L=f^*\theta$, choose a large
relative reduced divisor $D$ away from the ramification and resolve
$R\pi_*L$ by

$$
 \pi_*L(D)\longrightarrow\pi_*(L(D)|_D).
$$

Cancel an invertible block. Since $h^0(L)=1$ and $\chi(L)=-1$,
the residual complex is a column $r:\mathcal O\to\mathcal O^2$.
Its domain generator lifts to a family of meromorphic sections $s$.
At the divisor sections write
$s=\alpha_i/x_i+\beta_i+\cdots$ and
$\sigma=df=q_i(x_i)\,dx_i$ in local spin frames. The residue theorem
gives an identity on the whole neighborhood, before imposing $r=0$:

$$
 0=\sum_i\operatorname{Res}(\sigma s^2)
  =\sum_i\bigl(2q_i(0)\beta_i+q_i'(0)\alpha_i\bigr)\alpha_i.
                                                                  \tag{6}
$$

The principal-parts vector is a holomorphic linear image of $r$.
Thus (6) gives a holomorphic row annihilating $r$. At the central
kernel section $u$ it is twice the Serre-duality functional of
$\sigma u\ne0$. The row is nonzero, so the residual complex becomes
$\mathcal O\xrightarrow{(g_\theta,0)}\mathcal O^2$.
The determinantal ideal is exactly $(g_\theta)$. On its zero scheme
the kernel is a free line, including after Artin base change.

For an Artin deformation with both kernels, the two lifted sections
give a trace-free map $a^{-1}\oplus b^{-1}\to T$, where
$f_*\mathcal O_Y=\mathcal O_C\oplus T$. Its determinant has a
simple zero on the central curve. The implicit-function theorem over
the Artin base therefore identifies its divisor with a unique relative
point $p$. Dualizing gives

$$
 0\longrightarrow E=T^*\longrightarrow a\oplus b
                  \longrightarrow\mathcal O_p\longrightarrow0. \tag{7}
$$

Both quotient coordinates are units. The polynomial cubic-algebra
construction of Miranda, applied over the same Artin ring, gives the
unique possible nonzero datum up to the two diagonal changes of frame:

$$
 \lambda\bigl(BX^3-3rAX^2Y+3r^2BXY^2-r^3AY^3\bigr),
 \quad A=S_a(\cdot,p),\quad B=S_b(\cdot,p),\quad r,\lambda\in R^*.
                                                                  \tag{8}
$$

In the modification frame $k=X-rY$, $l=zY$, its four coefficients
are $\lambda(zB,3r(B-A),6r^2(B-A)/z,4r^3(B-A)/z^2)$.
The residue-one kernels have zero constant term. Consequently the
only extension condition is $Q(p)=0$ in the Artin ring, with its
actual nilpotents. The diagonal action on $(r,\lambda)$ is free and
transitive, so no extra parameter or infinitesimal stabilizer remains.
The inverse construction recovers $p$ as the determinant divisor.
Hence

$$
 \widehat{\mathcal O}_{\mathcal F,f_p}/(g_a,g_b)
       \simeq\widehat{\mathcal O}_{C,p}/(Q).                    \tag{9}
$$

Finite colength in a regular surface makes each $g_\theta$ a
nonzero Cartier equation near the common points. At a simple zero
of $Q$, their differentials are independent. The local multiplicity
is $\operatorname{ord}_pQ$, and the total is four. Calculation C1 gives
trivial deck groups, so no additional weight modifies this length.

We also need transversality to the triple-ramification divisor
$\mathcal T\subset\mathcal F$. On this divisor the discriminant
line is $\eta=ab^{-1}$. At its branch point $q$, the trace pairing
on $T$, in the algebra $w=z^3$, has matrix
$\left(\begin{smallmatrix}0&3w\\3w&0\end{smallmatrix}\right)$.
After division by $q$ it is nondegenerate. Since
$\det T=\eta\mathcal O(-q)$, this identifies $T$ with $T\otimes\eta$
in families on $\mathcal T$. It identifies the two cohomology-jump
ideals restricted to $\mathcal T$. At a simple common point their
common restriction must have order one: otherwise both independent
differentials in (9) would annihilate $T\mathcal T$. Thus each
$K_\theta$ meets $\mathcal T$ transversely.

The relative deformation calculation immediately below establishes the local scheme structure. The global single-spin dimension bound is proved in C3.
No global reduced-divisor formula for the individual $K_\theta$
is being inferred from their behavior at the common points.

## 3. Exhaustion of the compact common support

The smooth classification extends to finite reduced sources as follows.
For a connected reduced finite cubic $u:Z\to C$,
$u_*\mathcal O_Z=\mathcal O_C\oplus T$ and
$\deg T=4-p_a(Z)$. Suppose both spin sections are nonzero.
If $a^{-1}\oplus b^{-1}\to T$ has rank one, its joint graph lies
in a line $N\subset a\oplus b$ of degree at most zero. The graph
is a reduced monic spectral cubic, and $Z$ maps finitely and
birationally to it. This remains true for normalization degrees
$2+1$: the generic section values are $s,-s,0$, whereas three
degree-one components have no sections. Therefore

$$
                    p_a(Z)\le4+3\deg N\le4.                    \tag{10}
$$

This is an arithmetic-genus inequality, not merely a bound on the
normalization. In the rank-two case, genus six forces
$T=a^{-1}\oplus b^{-1}$, whose Miranda cubic coefficient space is
$H^0(b\oplus a\oplus b\oplus a)=0$. Genus five gives (7).
An axial quotient permits only a pure cube, incompatible with a
generically reduced algebra. A nonaxial quotient gives (8), and the
actual multiplication algebra in Calculation C1 is smooth. These arguments
also cover reducible and non-Gorenstein finite sources: the trace-zero
rank-three multiplication identities are polynomial and do not require
irreducibility. Thus every finite common-kernel genus-five source is
smooth, and no finite genus-six source has both kernels.

For a stable map $X\to C$ of genus six, contract each maximal
connected ghost block to obtain the finite reduced curve
$Z=\operatorname{Spec}_C f_*\mathcal O_X$. If the blocks have genera
$q_i$ and attachment numbers $a_i$, then

$$
 r+t+\sum_i(q_i+a_i-1)=2,\qquad
 p_a(Z)=6-\sum_iq_i,
 \qquad H^0(X,f^*\theta)=H^0(Z,u^*\theta),                     \tag{11}
$$

where $2r$ is total ramification and $t$ counts direct nodes on the
nonconstant normalization. Stability allows only $(1,1),(2,1),(1,2)$
and $(0,3)$. The finite genus-six exclusion removes rational tripods
as well as finite maps. An elliptic bridge becomes a genus-five node,
which the smoothness conclusion removes. The same argument removes
a nodal core with an elliptic tail, including the rational insertion
when its marking reaches a node. Two elliptic tails or one genus-two
tail force a connected étale core and belong to the removed full
étale summands. The only primitive common kernels are therefore
$Y_p$ with a genus-one tail attached at an arbitrary point.

For membership in the main locus the attachment must be ramified.
Indeed, the necessary smoothing condition of
[Rezaee--Swaminathan, Theorem 1.7](https://arxiv.org/html/2407.01845v2#S1.Thm7)
applies to a nodal ghost as well. For a one-pointed stable genus-one
curve $(E,e)$ the connecting map
$\delta:T_{E,e}\to H^1(E,\mathcal O_E)$ is an isomorphism, by
Serre duality and the nonvanishing dualizing differential at $e$.
The theorem then forces $df$ at the attachment to vanish. This rules
out higher-order smoothings at an unramified point, not merely their
first-order tangent vectors. Each $Y_p$ has exactly one ramification
point, of index three. Calculation C1 proves sufficiency by a fixed-target
admissible-cover smoothing and includes the nodal tail by properness.
This proves (3).

## 4. The main chart and the two kernel surfaces

Assume henceforth that $Q$ has simple zeros and fix $f_p$. The
finite-cover deformation has the local $A_2$ form
$w^3+tw+v$ at its triple point. Its two parameters are $t,v$;
$\mathcal T$ has equation $t=0$. Mark the attachment at $w=x$,
and let $\tau$ smooth the node joining it to a pointed genus-one
curve. On an orbifold chart of $\overline{\mathcal M}_{1,1}$ these
give five smooth ambient parameters, including the elliptic parameter.
The normalization of the tangential deformation complex leaves one
obstruction, $H^1(E,\mathcal O_E)\otimes T_qC$.

The entire node-preserving family exists, so the one-equation versal
chart has equation $F=\tau B$. Its first derivative in the smoothing
direction on $\tau=0$ is the principal-part obstruction
$\delta\,df_x$. Thus, up to a unit,

$$
            B|_{\tau=0}=t+3x^2,\qquad
            \frac{\partial B}{\partial t}\ne0.                \tag{12}
$$

This calculation is valid also at the nodal pointed genus-one curve:
its one-parameter moduli chart is smooth and the same $\delta$ is an
isomorphism. The unique branch containing smooth sources is $B=0$,
a smooth fourfold chart of $M$. Its tail divisor $D$ is Cartier,
with equation $\tau=0$, and

$$
 D=R\times\overline{\mathcal M}_{1,1},\qquad
 R:\ t=-3x^2,\qquad R\longrightarrow\mathcal F\text{ has degree }2.
                                                                  \tag{13}
$$

These statements are local near the complete curve $E_p$, with
equivariant changes of charts. They make no assertion that the
whole of $M$ is smooth. The smoothing coordinate is intrinsic as a
section of the tensor product of the two tangent lines; consequently

$$
             N_{D/M}|_{E_p}=T_{r_p}Y_p\otimes\mathbb L_e^\vee,
             \qquad c_1(N_{D/M}|_{E_p})=-\psi_e=-\lambda.       \tag{14}
$$

Both pulled-back spin sections vanish simply at $r_p$. This follows
from the actual cubic algebra in Calculation C1: near $q=\iota p$, its
generator $U$ satisfies $U^3-3dU+2Ad=0$, and the other trace-free
generator has a nonzero linear term in $U$, since $A(q)/B(q)=\pm1$.
The two sections extend over the tail with value zero.

Let $Z_\theta$ denote the local holomorphic-kernel locus on $M$.
The divisor resolution and residue identity (6) work on the nodal
family as well, using $\sigma:f^*K_C\to\omega_\pi$. They reduce
the rank-one-to-rank-three cohomology column to two equations.
On the tail divisor the first equation is $g_\theta|_R$.
The first-order obstruction to extending the kernel section over
the smoothing is
$\tau u_\theta'(r_p)\delta$ in the elliptic $H^1$ summand, and
its coefficient is nonzero. The core variation and this elliptic
variation are independent. Since $K_\theta$ is transverse to
$\mathcal T$, its pullback to $R$ is smooth. It follows, including
scheme structure, that

$$
             Z_\theta:\quad\tau=0,\quad g_\theta|_R=0         \tag{15}
$$

is a smooth surface near $E_p$. One may also see exact equality as
follows: the two independent equations define a smooth surface, and
the displayed tail surface is a contained smooth surface of the same
dimension with the same tangent space. The implicit-function theorem
identifies their germs. This rules out an additional local kernel
branch leaving the tail divisor.

## 5. The supported Euler correction

The comparison needed here is a comparison of families, not just of
individual cokernel fibers. Use the nodal Fredholm construction in
[Calculation B2](../excess/genus-five-mixed-spin-double-covers.md).
It applies unchanged to degree three: Siebert's cohomology realization,
norm-continuity and uniform sup estimates identify the Dolbeault family
with its finite cohomology complex; the plumbing estimate for the
zeroth-order spin term makes $\bar\partial+sA(df)$ continuous also
at $s=0$. Stabilization by forms supported away from nodes identifies
the positive-mass quotient charts and their determinant orientations
with those of the actual real obstruction bundle. The construction
retains finite stabilizers. The relevant sources are
[Siebert, §4](https://arxiv.org/abs/dg-ga/9608005v2),
[Zinger, §2.4](https://arxiv.org/html/0807.0805v3#S2.SS4), and
[Lee--Parker, §§8, 11--12](https://users.math.msu.edu/users/parker/LeeParkerObstructionBundle.pdf).

We specify the supported class to avoid a possible ambiguity. The
mass homotopy gives a stable oriented real identification between
$V_\theta$ and the realification of $W_\theta$. Classify this datum
by the homotopy pullback

$$
                 P=BU\times^h_{BSO}BSO(4),                    \tag{16}
$$

where $BU$ is the virtual-rank-two component. Over the kernel-free
complement, the same homotopy gives a reduction to $BU(2)$: the
Dolbeault cokernel is an actual complex rank-two bundle. On this
subspace $c_2-e=0$. Since $H^3(BU(2);\mathbb Q)=0$, the relative
cohomology exact sequence gives a unique universal lift of $c_2-e$
to $H^4(P,BU(2);\mathbb Q)$. Pulling it back along the family and
dividing by two defines a supported lift

$$
       \widetilde\Delta_\theta\in
       H^4(M,M\setminus\mathcal K_\theta;\mathbb Q),            \tag{17}
$$

where $\mathcal K_\theta$ is the full holomorphic-kernel support.
The specified homotopy is part of this construction. Equality of
absolute characteristic classes away from the kernel, by itself,
would not determine a supported lift.

At every point of (15), the transverse complex variations have rank
two. The spin mass pairs nontrivially with $\sigma u$, while (6)
makes the two complex variations orthogonal to that functional.
After cancellation the real reduction is homotopic to

$$
             u\longmapsto(z_1u,z_2u,\overline u):
                         \mathbb C\longrightarrow\mathbb C^3.\tag{18}
$$

Its quotient has the opposite orientation to the first two complex
coordinates at the central point. The projected last-coordinate
section has local Euler degree $-1$, compared with $+1$ for the
complex reduction replacing $\overline u$ by $u$. The relative
Euler difference is $-2$ times the complex Thom class. Therefore
$\widetilde\Delta_\theta$ is the Thom class of $Z_\theta$ near
the whole of $E_p$, with coefficient $+1$. The nodal elliptic point
cannot carry an additional degree-four correction: it has complex
codimension at least three in this smooth neighborhood, so such
supported degree-four cohomology vanishes. Equivalently (18) itself
holds on its uniformizing chart.

The product of the two supported lifts is supported on the common
locus (3). Excision therefore computes its pairing with $[M]$
entirely in the smooth neighborhoods just described. No assertion
that every individual kernel component elsewhere is a smooth surface,
or that it has coefficient one, is used.

## 6. The compact excess intersection

On $\mathcal F$, $g_a,g_b$ are independent local parameters at
$f_p$. Pullback by the finite flat ramification-marking map in (13)
makes their intersection have length two. In suitable linear
coordinates their leading form is

$$
 g_a=v+\alpha t+\cdots,\qquad
 g_b=v+\beta t+\cdots,\qquad \alpha\ne\beta,
 \qquad t=-3x^2.                                                \tag{19}
$$

Thus the local intersection on $R$ has ideal equivalent to
$(v,x^2)$. More intrinsically, finite flat pullback doubles the
length-one intersection in (9). The two surfaces $Z_a,Z_b$ are
tangent inside $D$; their intersection is a doubled $E_p$, not a
transverse intersection along the reduced curve.

Let $i:D\hookrightarrow M$ on this neighborhood. The Cartier-divisor
self-intersection formula, applied to the two divisor classes on
$D$, gives

$$
 [Z_a]\,[Z_b]
 =i_*\left(c_1(N_{D/M})\,
          [\rho^*K_a]\,[\rho^*K_b]\right),                    \tag{20}
$$

where $\rho:D\to\mathcal F$ forgets the ramification marking and
the tail. This formula does not require a clean intersection of
$Z_a$ and $Z_b$ in $M$. Their proper divisor intersection inside
$D$ is $2[E_p]$ with its scheme multiplicity. Equations (14), (17)
and (20) consequently give

$$
 \operatorname{Contr}_{E_p}(\Delta_a\Delta_b)
       =-2\int_{\overline{\mathcal M}_{1,1}}\lambda
       =-\frac1{12}.                                          \tag{21}
$$

The left side denotes the contribution localized to $E_p$ by the
supported cup product in (17). The Hodge degree is the stack degree
$1/24$: the presentation $\overline{\mathcal M}_{1,1}=\mathbb P(4,6)$
identifies the Hodge line with $\mathcal O(1)$, of degree $1/(4\cdot6)$.
There is no core deck
factor; the pointed-tail stabilizers are exactly those already
included in this integral. Summing the four disjoint contributions
proves (4).

Finally the identity of classes
$e(V_a)e(V_b)=\tfrac12(e(V_a)^2+e(V_b)^2)-2(\Delta_a-\Delta_b)^2$
and the equal-spin main integral $161/6$ from Calculation B3 yield
$I_3=161/6-2(S_a^{\mathrm{self}}+S_b^{\mathrm{self}})+4(-1/3)$.
Subtracting $43/2$ proves (5).

## Deformation of the joint kernel scheme

### The joint deformation functor is the zero scheme of Q

Let $R$ be a local Artin $\mathbb C$-algebra. A deformation of a
common cover lying in $K_a\cap K_b$ has both kernel lines free by
(5). Choose generators, and split the finite flat algebra by its
trace:

$$
       f_*\mathcal O_{\mathcal Y_R}=\mathcal O_{C_R}\oplus T_R.
$$

The generators give $a_R^{-1}\oplus b_R^{-1}\to T_R$. Its
central determinant has one simple zero. Therefore its relative
determinant vanishes on the graph of a unique point $p_R\in C(R)$,
and its cokernel is an invertible sheaf on that graph. To see this
without a flatness assumption on a putative cokernel, use an
invertible matrix entry at the central point to diagonalize the
matrix locally to $(1,h)$; the derivative of $h$ in the curve
coordinate is a unit. The implicit-function theorem, including
over Artin bases, writes $h$ as a unit times $z-p_R$.

After dualizing there is the family elementary modification

$$
 0\longrightarrow E_R\longrightarrow a_R\oplus b_R
                    \longrightarrow\mathcal O_{p_R}\longrightarrow0,
 \tag{6}
$$

with a harmless invertible factor on $\operatorname{Spec}R$
suppressed. Both entries of the quotient are units because the
central quotient is nonaxial. Write their ratio as $r\in R^*$.

Miranda's cubic algebra formulas are polynomial identities with
denominators only $2$ and $3$. They therefore apply to (6) over
$R$, not only to its closed fiber. In particular the cubic datum
is a section of $\operatorname{Sym}^3E_R\otimes\det E_R^{-1}$.

For each noneffective theta characteristic,
$H^0(\theta(p_R))$ is a free rank-one $R$-module, generated by its
normalized Szegő section: this follows either by cohomology and
base change or by the residue exact sequence, since
$H^0(\theta)=H^1(\theta)=0$. Thus the residue argument of [Calculation C1](cubic-covers-with-two-spin-sections.md)
works over $R$ without an unaccounted infinitesimal coefficient.
The cubic is exactly

$$
 \lambda\bigl(BX^3-3rAX^2Y+3r^2BXY^2-r^3AY^3\bigr),
 \qquad \lambda\in R^*.
 \tag{7}
$$

In the modified frame $k=X-rY$, $l=zY$, its coefficients are

$$
 \lambda\left(zB,\;3r(B-A),\;
             6r^2\frac{B-A}{z},\;
             4r^3\frac{B-A}{z^2}\right).
 \tag{8}
$$

The normalized kernels have expansions
$A=z^{-1}+c_a(p_R)z+\cdots$ and
$B=z^{-1}+c_b(p_R)z+\cdots$. The only nonregular term of (8) is
a unit times $(c_b-c_a)/z$ in its last entry. Therefore extension
holds over $R$ if and only if $Q(p_R)=0$ in $R$. There is no
replacement of that equation by its radical.

Changing the two section generators is the diagonal torus action

$$
 r\longmapsto tr/s,\qquad
 \lambda\longmapsto\lambda s^2/t.
$$

It acts freely and transitively on the two units $(r,\lambda)$:
there is a unique choice making both equal to one. Hence no scalar
parameter or residual infinitesimal stabilizer survives. Conversely,
(6)--(8) construct a finite flat cover from every $p_R$ satisfying
$Q(p_R)=0$. Its source is a deformation of the smooth central
source, so it is a family in the smooth-source neighborhood under
consideration. Recovering the determinant point and reconstructing
the cubic are inverse operations, and commute with Artin base change.

This proves an equivalence of the unframed deformation functors,
not merely a bijection of closed points. Consequently, for the
cover $f_p$ corresponding to a zero $p$ of $Q$,

$$
 \widehat{\mathcal O}_{\mathcal F,f_p}/(g_a,g_b)
       \simeq \widehat{\mathcal O}_{C,p}/(Q).
 \tag{9}
$$

The use of Miranda here is justified by the displayed polynomial
construction. The irreducibility hypotheses in his original global
statement are not silently being applied to a nonreduced Artin base.

### A global one-dimensional-kernel bound

For every finite connected nodal cubic map of arithmetic genus
five and every noneffective theta characteristic $\theta$ on $C$,

$$
                         h^0(Y,f^*\theta)\le1.
 \tag{13}
$$

Indeed write $f_*\mathcal O_Y=\mathcal O_C\oplus T$ as before;
$\deg T=-1$. Suppose two independent sections give
$\theta^{-1}\oplus\theta^{-1}\to T$. If its generic rank is one,
let $N\subset T$ be its saturated image. Every line subbundle of
$f_*\mathcal O_Y$ has degree at most zero: a positive-degree
subbundle would, by adjunction, give a nonzero section of
$f^*N^{-1}$, which has negative degree on every positive-degree
component of the finite source. Hence $\deg(N\theta)\le1$,
and on a genus-two curve it has at most one independent section,
a contradiction.

If the generic rank is two, its determinant has degree one, so
dualizing gives an elementary modification of
$\theta\oplus\theta$ at one point $p$. The constant
$\mathrm{GL}_2$-action makes its quotient axial, and

$$
                         E\simeq\theta\oplus\theta(-p).
$$

The four cubic coefficient spaces are, in some order,
$H^0(\theta(p))$, $H^0(\theta)$,
$H^0(\theta(-p))$, and $H^0(\theta(-2p))$.
Only the first is nonzero. Every cubic datum is therefore a pure
cube, with identically zero discriminant, whereas the source is
reduced and the generic algebra is a product of separable field
extensions. This contradiction proves (13).

The proof covers reducible and irreducible finite nodal sources.
It is a bound on the fiber dimensions. It does not alone prove
that every individual global kernel component is a generically
reduced divisor, or that no entire connected component of
$\mathcal F$ has persistent kernel.
