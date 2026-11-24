# C3. Single-spin cubic kernels and singular spectral curves

## 1. The spectral calculation and its intersection classes

Fix a smooth genus-two curve $C$ and a noneffective even theta
characteristic $\theta$. Let $\mathcal F=\mathcal F_3(C)$ be the
proper smooth surface of connected finite cubic maps of arithmetic
genus five, and let $M=M_3$ be the proper reduced main four-cycle of
genus-six cubic maps, as in
[Calculation B3](../excess/second-excess-reduction-in-all-degrees.md) and
[the proof](cubic-common-kernel-excess-intersection.md).
The index $W_\theta=-R\pi_*f^*\theta$ has virtual complex rank one
on $\mathcal F$ and two on $M$. The positive-mass real spin cokernel
$V_\theta$ has the corresponding real ranks two and four.

On $M$ denote the individual square by

$$
 S_\theta=\int_{[M]}\Delta_\theta^2,
 \qquad \Delta_\theta=\frac{c_2(W_\theta)-e(V_\theta)}2.
                                                                  \tag{1}
$$

Calculation C2 gives $n_{3,6}=4-2(S_a+S_b)$ for distinct even spins when
their Szegő quadratic has four simple zeros. This calculation proves three
separate compact calculations relevant to (1).

**Theorem 1.** The cohomology-jump scheme on $\mathcal F$ is an
effective Cartier divisor $D_\theta$, with its determinantal
multiplicities. Globally on that surface,

$$
 e(V_\theta)=c_1(W_\theta)-2[D_\theta],\qquad
 D_aD_b=4\ (a\ne b),\qquad D_\theta^2=-12.                     \tag{2}
$$

Reducedness of $D_\theta$ is not assumed. The first-excess divisor
square in (2) is different from the genus-six number $S_\theta$.

**Theorem 2.** On the entire primitive genus-six stable-map summand,
$h^0(X,f^*\theta)\le1$. The proper weighted incidence surface $I$
of a singular point of a trace-free spectral cubic, defined in §3,
is reduced and irreducible and is a local complete intersection.
It carries a proper flat family $\rho:X_I\to I$ of connected
reduced planar curves of arithmetic genus six. Set
$\mathcal K=\mathcal O_I(-1)$, $\xi=c_1(\mathcal K^*)$ and
$h=c_1(\theta_p)$, where $p:I\to C$ marks the singular point.
Then

$$
 \begin{aligned}
 \int_I\xi^2&=\frac43,&\qquad \int_I\xi h&=\frac12,\\
 W_I=-R\rho_*f^*\theta
   &=\mathcal K+2\mathcal K^2-\theta_p^{-3}\mathcal K^2,\\
 c_1(W_I)&=-3\xi+3h,& c_2(W_I)&=2\xi^2-3\xi h,\\
 \int_I\bigl(c_2(W_I)-\xi c_1(W_I)+\xi^2\bigr)&=5.
 \end{aligned}                                                   \tag{3}
$$

The number $5$ is an intersection on this specified family, not a
value for $S_\theta$. The family is not proved nodal everywhere.
Moreover the full individual kernel locus has contracted-tail and
bridge components outside the closure of its smooth-source locus.
Article V, Corollary 1.2 computes the cubic exponents from the partition formula.

**Theorem 3.** Let $R\to\mathcal F$ be the degree-two ramification
scheme, including its length-two fibers at triple points and nodes,
and put $\Gamma_\theta=(R\times_{\mathcal F}D_\theta)
\times\overline{\mathcal M}_{1,1}$. Attaching an elliptic tail
at the marked ramification point gives a proper morphism
$j_\theta:\Gamma_\theta\to M$, with a rational insertion when
the marking reaches a node. Then

$$
                 \int_{[\Gamma_\theta]}j_\theta^*\Delta_\theta
                       =-\frac{D_\theta^2}{12}=1.              \tag{3a}
$$

This is a test-family integral. Identification of
$j_{\theta*}[\Gamma_\theta]$ as the relevant summand of the full
supported kernel cycle is an additional statement, not part of (3a).

## 2. Kernel bounds and the spectral graph

For a connected reduced finite cubic $f:Z\to C$, including reducible
or non-Gorenstein sources, write
$f_*\mathcal O_Z=\mathcal O_C\oplus T$. Then
$\deg T=4-p_a(Z)$. Two independent spin sections would define
$\theta^{-1}\oplus\theta^{-1}\to T$. In rank one, saturating the
image gives a line $N\subset f_*\mathcal O_Z$ with $\deg N\le0$:
a positive-degree line would pull back to a line of negative degree
on every component with a nonzero section, which is impossible.
Thus $h^0(N\theta)\le1$, a contradiction.

In rank two and genus six the map is an isomorphism. Its dual is
$E=\theta\oplus\theta$, whose cubic multiplication datum lies in
$H^0(\operatorname{Sym}^3E\otimes\det E^{-1})=H^0(\theta)^{\oplus4}=0$.
The zero datum gives a generically nonreduced algebra. In genus five
the dual is an elementary modification of $\theta\oplus\theta$
at one point. A constant change of basis makes it
$\theta\oplus\theta(-p)$. Only the coefficient space
$H^0(\theta(p))$ is nonzero, so the cubic is a pure cube, again
not generically reduced. These arguments use the trace-free
multiplication identities of Miranda, which are polynomial and apply
without an irreducible-source hypothesis.

For a primitive genus-six stable map, contract its maximal connected
ghost blocks to $Z=\operatorname{Spec}_C f_*\mathcal O_X$. The genus
budget and projection formula of Calculation C2 give

$$
 r+t+\sum_i(q_i+a_i-1)=2,\qquad
 p_a(Z)=6-\sum_iq_i,\qquad
 H^0(X,f^*\theta)=H^0(Z,f_Z^*\theta).                            \tag{4}
$$

If $\sum q_i\le1$, the preceding finite bounds apply. If it is
two, the budget forces a connected étale cubic core and only
one-attachment ghosts. This is a full degree-one-upstairs étale
summand, excluded from the primitive part. This proves the bound
in Theorem 2, including rational tripods and disconnected nonconstant
normalizations.

A nonzero spin section on a finite cubic has trace zero and has a
reduced monic graph

$$
 S:\quad y^3+A y+B=0,\qquad
 A\in V_2:=H^0(C,K_C),\quad B\in V_3:=H^0(C,\theta^3),
 \qquad \dim V_2=\dim V_3=2.                                  \tag{5}
$$

The map to the graph is finite and birational componentwise. On an
irreducible cubic this follows from prime degree and noneffectivity;
on a reducible cubic the nonzero generic values are $s,-s,0$ on the
degree partition $2+1$. The graph algebra has underlying module
$\mathcal O_C\oplus\theta^{-1}\oplus\theta^{-2}$, hence
$p_a(S)=7$. Every nonzero coefficient pair gives a reduced spectrum:
a repeated root would be an integral meromorphic section of $\theta$,
hence a holomorphic one, and therefore zero; the repeated-root
equations would then force $A=B=0$. The same argument shows that
reducibility is equivalent to $B=0$.

Near any primitive genus-six kernel point the residue identity of
Calculation C2 reduces a cohomology complex to
$\mathcal O\xrightarrow{(r_1,r_2,0)}\mathcal O^3$. Thus the
kernel scheme has two equations and a free kernel line on its zero
scheme. This does not, by itself, assert a regular sequence or a
reduced surface everywhere on the possibly singular main locus.

## 3. The proper singular-point incidence

Use the weighted projective stack of lines

$$
 P=\mathbb P_C\bigl(\theta[1]\oplus V_2[2]\oplus V_3[3]\bigr),
 \qquad\mathcal K=\mathcal O_P(-1),\quad\xi=c_1(\mathcal K^*).
                                                                  \tag{6}
$$

It parametrizes $(p,u,A,B)$ under
$(u,A,B)\mapsto(tu,t^2A,t^3B)$; in families
$u\in\theta_p\mathcal K^{-1}$. Let $I$ be the singular-point
incidence

$$
 A(p)+3u^2=0,\qquad B(p)-2u^3=0,
 \qquad A'(p)u+B'(p)=0.                                        \tag{7}
$$

The last expression is intrinsic on the zero scheme of the first
two. Indeed for $F=y^3+Ay+B$, those equations say $F=F_y=0$ at
$(p,u)$, so the horizontal differential is independent of the
frame and the horizontal lift. Its line is
$K_C|_p\otimes\theta_p^3\mathcal K^{-3}
=\theta_p^5\mathcal K^{-3}$.

Both $K_C$ and $\theta^3$ are basepoint-free pencils. A base point
of the latter would give $\theta\simeq\mathcal O_C(p)$, contrary
to noneffectivity. Their maps of degrees two and three jointly
embed $C$ in $\mathbb P^1\times\mathbb P^1$: a hyperelliptic
pair in one cubic-pencil fiber, or common ramification at a
Weierstrass point, would give the same contradiction.

After solving the first two equations of (7), the third has local
form

$$
                \beta(p)b+\alpha(p)u a+\gamma(p)u^3=0,         \tag{8}
$$

where $a,b$ are the remaining coefficient coordinates. The zeros
of $\alpha$ are the Weierstrass points, and those of $\beta$
are the ramification points of $|\theta^3|$. They do not coincide,
since $H^0(\theta^3(-2p))=H^0(\theta)=0$ at a Weierstrass point.
Thus the third equation is regular. For $u\ne0$ the coefficient
matrix in (7) has rank three; before quotient the incidence is an
affine-line bundle over $\operatorname{Tot}\theta\setminus0$.
The locus $u=0$ has dimension at most one: generically it forces
$B=0$ and one projective choice of $A$, and at the finitely many
zeros of $\beta$ it leaves a weighted projective line. Consequently
$I$ is an irreducible, generically reduced local complete intersection
surface, hence reduced. Its fundamental cycle is

$$
 [I]=(2\xi+2h)(3\xi+3h)(3\xi+5h)
       =18\xi^3+66h\xi^2.                                    \tag{9}
$$

The forgetful map to $\mathbb P(2,2,3,3)$ is proper and has finite
geometric fibers, since a reduced spectrum has only finitely many
singular points. Its image is an irreducible discriminant surface.
Its general spectrum has one ordinary node and is integral. To
check the first assertion, at general $p$ and $u\ne0$ the free
coefficient varies $A'(p)$ nontrivially; the Hessian determinant has
a nonzero quadratic term in that coefficient. To check uniqueness,
the normal to the discriminant at a node is

$$
              (\dot A,\dot B)\longmapsto u\dot A(p)+\dot B(p).
                                                                  \tag{10}
$$

Two nodes on a general discriminant member would give proportional
normals, contradicting the embedding by the two pencils. A cubic
cannot have two distinct singular roots over the same base point.
The reducible locus $B=0$ is proper as well. Thus a dense open
subset of $I$ parametrizes integral one-node spectra and their
smooth genus-six normalizations. The node's distinct branches make
normalization commute with Artin deformations preserving the node.
Conversely the determinant of the spectral-algebra inclusion for
a deformed section recovers its unique node. This identifies the
ordinary-node incidence locally with the determinantal kernel
scheme, not just with its closed points. It is a regular surface
there, with local coefficient $+1$ in the supported correction
$\Delta_\theta$.

## 4. The global divisor on the finite genus-five surface

A genus-five finite kernel source maps birationally to (5) with
algebra colength two. Its spectrum therefore lies outside the
one-node open subset just proved. Properness of the incidence makes
that complement a closed subset of the discriminant of dimension
at most one.

There are finitely many nodal intermediate sources over a fixed
reduced spectrum. Normalize the spectrum once. Every such source
identifies disjoint pairs of distinct points of this normalization
lying over the finite singular set. Its subalgebra consists exactly
of functions whose residue values agree at each selected pair.
There are no gluing scalars or jet parameters in a unital nodal
algebra. Thus only finitely many choices occur. General nonnodal
intermediate algebras would not obey this assertion. Since the
kernel is one-dimensional, its projectivization adds no parameter.
The entire genus-five kernel locus has dimension at most one;
no component of $\mathcal F$ has persistent kernel.

At any of its points, cancellation in a divisor resolution of
$R\pi_*f^*\theta$ gives a column $r:\mathcal O\to\mathcal O^2$.
The exact residue identity from Calculation C2 is valid also at finite
nodes. Locally $xy=s$ and
$df=(xF_x-yF_y)dx/x$ is dualizing; its branch residues cancel.
The identity supplies a holomorphic row annihilating $r$, nonzero
at the kernel because its Serre-duality functional is $2df\,u$.
The column becomes $(g,0)$. The dimension bound makes $g$ a
nonzero germ in the regular surface ring. Therefore the global
determinantal ideal defines an effective Cartier divisor
$D_\theta$.

The positive-mass Fredholm family, including finite nodal sources,
is the one compared in Calculation B2. Its stable real identification
with the rank-one complex index, together with the reduction to
an actual complex line off $D_\theta$, defines a supported lift
of $c_1-e$. In the rank-one version of Calculation C2's construction,
the classifying pair is
$BU_{[1]}\times^h_{BSO}BSO(2)$ with subspace $BU(1)$; uniqueness
of its universal relative lift follows from $H^1(BU(1);\mathbb Q)=0$.
Near a generic divisor point the real reduction is
$u\mapsto(g u,\overline u)$. If $g$ vanishes to order $m$, its
transverse Euler difference is $-2m$. Thus the divisor is taken
with exactly its determinantal multiplicities. Lower-dimensional
intersections cannot supply another degree-two supported class.
This proves the first identity of (2) globally.

The finite scheme intersection in Calculation C2 gives $D_aD_b=4$ for
distinct even spins, including colliding zeros. The rational index
class $c_1(W_\theta)$ is independent of $\theta$. Article V, §§4–6 gives
the finite integrals $\int e(V_\theta)^2=-70$ and
$\int e(V_a)e(V_b)=-6$. Consequently

$$
 (D_a-D_b)^2=\frac{-70-70-2(-6)}4=-32,
 \qquad D_a^2+D_b^2=-24.                                      \tag{11}
$$

Applying the last identity to any three distinct even theta
characteristics gives $D_\theta^2=-12$ for each. The identity applies to the determinantal divisors with their multiplicities.

## 5. A proper family of arithmetic-genus-six curves

Over $I$ consider the relative surface
$S_I=\operatorname{Tot}_{C\times I}(\theta\boxtimes\mathcal K^{-1})$
and the spectral divisor $S=(y^3+Ay+B=0)$. Blow up the section
$(p,u)$, writing $b:\widetilde S_I\to S_I$ and $E$ for its
exceptional divisor. Since the equation vanishes to order at least
two along that section, the residual divisor

$$
                         X_I=b^*S-2E                           \tag{12}
$$

is effective. Its equation is a nonzerodivisor in every smooth
surface fiber of the blowup, so $X_I\to I$ is flat. It is proper,
since the original monic spectrum is finite over the proper curve
$C$ and the blowup is projective.

At a double point this is the strict transform. At a triple point
it retains one copy of $E$. Every possible triple point is ordinary:
it has $u=0$, $\operatorname{ord}_pA\ge2$ and
$\operatorname{ord}_pB\ge3$. If $A\ne0$, then $p$ is Weierstrass,
$\operatorname{div}A=2p$, and $B=0$; the tangent cone is
$y(y^2+a x^2)$ with $a\ne0$. If $A=0$, then $B\ne0$ has divisor
$3p$, possible when $\theta^3\simeq\mathcal O_C(3p)$, and the
tangent cone is $y^3+b x^3$, $b\ne0$. In both cases the exceptional
component is a reduced rational tripod with three distinct
attachments. Dropping it would give the wrong stable-map boundary.

These observations prove reducedness of all fibers of (12).
They are connected: integral spectra remain connected, while the
reducible spectrum $y(y^2+A)=0$ has two intersections; resolving
one leaves the other, or its double intersection gives the retained
tripod. The genus and the family index follow from the exact
sequence below. The curves are planar, being Cartier divisors in
the smooth surface fibers. Nodality everywhere is a different
assertion and is not claimed.

Let $N$ be the normal bundle of the center section in $S_I/I$.
Its horizontal and vertical factors give

$$
 \det N=\theta_p^{-1}\mathcal K^{-1},\qquad
 \mathcal O(-S)|_{(p,u)}=\theta_p^{-3}\mathcal K^3.              \tag{13}
$$

For $E=\mathbb P_{\mathrm{lines}}(N)$,
$\mathcal O_E(E)=\mathcal O_E(-1)$ and
$R^1q_*\mathcal O_E(-2)=\det N$, by relative Serre duality.
The sequences for $E$ and $2E$, followed by the Cartier sequence
of (12), give

$$
 \begin{gathered}
 0\longrightarrow\mathcal O_S\longrightarrow b_*\mathcal O_{X_I}
       \longrightarrow\mathcal O_{(p,u)}\otimes Q\longrightarrow0,
 \qquad R^1b_*\mathcal O_{X_I}=0,\\
 Q=\det N\otimes\mathcal O(-S)|_{(p,u)}
                         =\theta_p^{-4}\mathcal K^2.            \tag{14}
 \end{gathered}
$$

This is an actual sheaf sequence. At triple points it is not full
normalization: the exceptional component retains agreement of the
three branch values. In every case it increases the Euler
characteristic by one, so the source arithmetic genus is six.

The spectral algebra before blowup is
$\mathcal O_C\oplus\theta^{-1}\mathcal K
\oplus\theta^{-2}\mathcal K^2$. Tensoring (14) by $\theta$ and
using $R\Gamma(\theta)=0$, $\chi(\mathcal O_C)=-1$ and
$\chi(\theta^{-1})=-2$ gives the index in (3). The quotient line
in (14) cannot be replaced by an untwisted skyscraper when computing
family Chern classes.

The tautological line $\mathcal K$ is the whole space of sections
on every fiber of (12). Indeed the dual of the connecting map
from $\theta_pQ$ in (14) restricts to evaluation at $p$ on
$H^0(C,\theta^3)\mathcal K^{-2}$, which is surjective. Thus the
connecting map is nonzero and no second section is acquired.
Cohomology and base change give a kernel line and a rank-three
first-cohomology bundle on $I$.

## 6. The characteristic number of this family

The rational Chow relation for the lines convention (6) is

$$
  (\xi+h)(2\xi)^2(3\xi)^2=0,
  \qquad p_*\xi^4=\frac1{36},\qquad \xi^5=-h\xi^4.             \tag{15}
$$

The denominator $36$ is the weighted stack degree, with all
stabilizers retained. Together with (9), it gives
$\int_I\xi^2=(66-18)/36=4/3$ and
$\int_I\xi h=18/36=1/2$. The index of (3) has Chern classes

$$
 c_1(W_I)=-3\xi+3h,\qquad c_2(W_I)=2\xi^2-3\xi h.             \tag{16}
$$

On a nodal family in the actual kernel locus, write $H=R^1\pi_*f^*\theta$
and $K=R^0\pi_*f^*\theta$. The constant-rank Lee--Parker comparison
gives $e(V_\theta)=-c_2(H-\overline K)$, while $W=H-K$.
For $c_1(K)=-\xi$ it follows that

$$
             \Delta_\theta=c_2(W)-\xi c_1(W)+\xi^2.             \tag{17}
$$

Hence the natural class on the whole algebraic family (12) that
agrees with this restriction where the source is nodal is
$\gamma_I=c_2(W_I)-\xi c_1(W_I)+\xi^2=6\xi^2-6\xi h$.
Its degree is

$$
                  \int_I\gamma_I=6\left(\frac43-\frac12\right)=5.
                                                                  \tag{18}
$$

This proves the last assertion of Theorem 2. It uses no conjectural
partition series and assigns no Euler class to a nonnodal source.

## 7. The ramification-tail test family

The differential of the universal finite map is a section of
$\omega_\pi\otimes f^*T_C$. Let $R$ be its zero scheme on the
universal source. At simple ramification it has a simple zero; at
triple ramification its local equation is $3z^2+t=0$. At a finite
node the universal chart is

$$
 uv=s,\qquad w=u+v+q,\qquad df=(u-v)\,du/u.
                                                                  \tag{18a}
$$

Thus $R$ has $u=v$ and $s=u^2$ there. These charts prove that
$R$ is smooth and proper and $R\to\mathcal F$ is finite flat of
degree two, with all stack stabilizers retained. The nodal special
fiber is $\mathbb C[u]/(u^2)$, not a reduced point. Stabilizing the
additional marking inserts a rational component at a node and
places the marking on that component.

Gluing the pointed elliptic curve first defines a morphism on the
whole reduced smooth stack $R\times\overline{\mathcal M}_{1,1}$.
Its image is in $M$: simple ramification uses an elliptic function
with sole pole $2e$ and an unused degree-one bubble sheet; triple
ramification uses a sole pole $3e$ as in Calculation C1. The fixed-target
admissible smoothing includes nodal cores and elliptic endpoints
by closure. The morphism from the reduced source therefore factors
through the reduced closed main locus. Restricting it to
$\Gamma_\theta$ gives the stated factorization even when
$D_\theta$ is nonreduced.

On the full scheme $D=D_\theta$, write $K,H$ for the core
cohomology bundles of ranks one and two. They exist over its
nilpotents because the residual complex $(g,0)$ becomes zero on
$g=0$. The residue pairing gives a surjection $H\to K^*$;
write $H'$ for its kernel line. The first infinitesimal neighborhood
of $D$ in $\mathcal F$ gives an intrinsic isomorphism

$$
 K\xrightarrow{\sim}H'\otimes I_D/I_D^2,
 \qquad \mathcal O(D)|_D\simeq K^*\otimes H'\simeq\det H.
                                                                  \tag{18b}
$$

In a residual chart the map sends a generator to $g$ times a
generator. This remains an isomorphism if $g=z^m$; it uses the
generator of $(z^m)/(z^{2m})$, not the ordinary derivative on
the reduced support. Consequently, with $\xi=-c_1(K)$,

$$
       c_1(W_{\mathrm{core}})-\xi=c_1(\mathcal O(D)|_D).       \tag{18c}
$$

On $\Gamma_\theta$, normalization at the elliptic attachment
cancels the constant matching term, giving
$W=W_{\mathrm{core}}+\operatorname{ev}^*\theta\otimes\mathbb E^*$,
where $\mathbb E$ is the elliptic Hodge line. Rational insertion
does not change this formula. The kernel line remains $K$.
Apply the actual constant-kernel Euler comparison (17). Products
of two positive-codimension core classes vanish on the base curve
$R\times_{\mathcal F}D$, so

$$
 j_\theta^*\Delta_\theta
       =-\rho^*c_1(\mathcal O(D)|_D)\,\lambda.
                                                                  \tag{18d}
$$

Finite-flat pushforward has degree two, and
$\int_{\overline{\mathcal M}_{1,1}}\lambda=1/24$. Therefore its
degree is $-2D^2/24=1$, proving Theorem 3. This proof evaluates
the specified fundamental cycle with all determinantal
multiplicities; it makes no assertion that it is already the
required component of $\Delta_\theta\cap[M]$.
