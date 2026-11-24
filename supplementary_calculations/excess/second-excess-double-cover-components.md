# B1. The components of the second-excess double-cover integral

Let $C$ be a smooth compact curve of genus $g\ge2$, let $a,b$ be theta
characteristics, and use the direct sum of their positive real spin
operators as the normal operator $D$. Write $V_a,V_b$ for the oriented
real cokernel bundles on stable maps, with determinant orientations
continued from the complex Dolbeault operators. The connected degree-two
local invariant in source genus $2g+1$ is

$$
 L_{2,2g+1}(D)=\int_{[\overline{\mathcal M}_{2g+1}(C,2)]^{\rm vir}}
                         e(V_a)e(V_b).
$$

The genus $2g-1$ is the unramified minimum, so this is excess two.
Let $M$ be the reduced closure of the locus of smooth double covers with
four simple branch points, and set

$$
 I_g(a,b)=\int_{[M]}e(V_a)e(V_b).
$$

This integral is over a proper stable-map component. It is not an integral
over the open finite-map locus. The spin obstruction bundle is the one
constructed by Lee--Parker, restricted from the full stable-map space.
[Lee--Parker, Theorem 8.1 and §§11--12](https://users.math.msu.edu/users/parker/LeeParkerObstructionBundle.pdf)

## Component reduction

Put $s(t)=2\sin(t/2)$, let $c_2(D)$ be the signed stack count of connected
étale double covers, and let $F_2(D)$ be the finite-map integral at first
excess, on the compact surface of source genus $2g$. Then

$$
 \boxed{L_{2,2g+1}(D)=I_g(a,b)-\frac{g-1}{6}F_2(D)
             +c_2(D)[t^{4g}]s(t)^{4g-4}.}
 \tag{1}
$$

For distinct spin lines, put $\epsilon=(-1)^{h^0(a)+h^0(b)}$.
The previously proved values are $F_2(D)=0$ and $c_2(D)=-1/2$.
The elementary exponent at second excess therefore satisfies

$$
 \boxed{n_{2,2g+1}=I_g(a,b)
            -\epsilon\,2^{2g-6}(g-1)(g-2).}
 \tag{2}
$$

In particular, for a genus-two positive distinct-spin core,

$$
                 L_{2,5}=I_2(a,b)-\frac1{160},\qquad
                 n_{2,5}=I_2(a,b).
 \tag{3}
$$

The remaining integer is an Euler integral on the compactified four-branch
component. Equations (1)--(3) do not assert its vanishing.

### The stable-map strata

Normalize the nonconstant components of a stable map. Let $2r$ be their
total ramification, $t$ the number of direct nodes between nonconstant
branches, and let $G_i$ be the maximal connected contracted subcurves.
Write $q_i=p_a(G_i)$ and $a_i$ for their numbers of attachments. Since
the total degree is two, $a_i\le2$: distinct attachments use distinct
points of the normalized cover over one target point. Stability gives
$2q_i-2+a_i>0$, hence $q_i\ge1$. Riemann--Hurwitz and the dual-graph
genus formula give

$$
                 r+t+\sum_i(q_i+a_i-1)=2.
 \tag{4}
$$

Thus finite maps have $(r,t)=(2,0),(1,1),(0,2)$. The other possibilities
are an elliptic tail on a two-branched cover; an elliptic tail on a
one-node finite cover; two elliptic tails on a connected étale cover;
a genus-two tail on a connected étale cover; or an elliptic bridge
joining the two sheets of an étale normalization. The last normalization
may be connected or two copies of $C$. Contracted blocks include their
nodal stable limits. This convention includes rational components carrying
several positive-genus branches.

The image of the source fundamental group in $\pi_1(C)$ is locally
constant in stable-map families. Consequently the stack is the disjoint
union of its open and closed étale-lift summands and its primitive
summand, where primitive means surjective on $\pi_1(C)$.
This is the decomposition of Bryan--Pandharipande, not primitivity in
homology or in the cycle algebra. In degree two, a nontrivial étale lift
has degree one upstairs. Only the two-elliptic-tail and genus-two-tail
types belong to these summands. A bridge joining conjugate sheets does
not lift: its constant value upstairs would have to be two different
points. [Bryan--Pandharipande, Definition 2.2 and Proposition 2.3](https://arxiv.org/html/math/0009025v3#S2.SS2)

Let $\mathcal F_2$ be the proper finite-map surface at first excess and
$\pi:\mathcal U\to\mathcal F_2$ its universal curve, viewed as the stack of
one-pointed stable maps. A marking at a node is carried by an inserted
rational component. Gluing a pointed elliptic curve defines a proper
map

$$
 \gamma:\mathcal U\times\overline{\mathcal M}_{1,1}
       \longrightarrow\overline{\mathcal M}_{2g+1}(C,2).
$$

Its image $G$ is the elliptic-tail component, and its generic stack
degree is one. Every primitive map lies in $M\cup G$. Pair collisions
produce the finite nodal maps. A three-branch collision replaces a cusp
by an elliptic tail at ramification, and a four-branch collision replaces
a tacnode by an elliptic bridge. These limits follow by introducing a
rational target bubble and taking its degree-two cover. Conversely
$|2p|$ on a pointed elliptic curve and $|p+q|$ on a twice-pointed one
produce the required admissible covers. Smoothing their nodes with
stabilized target $C$ gives every smooth bridge in $M$; nodal bridges
follow by degeneration. Tails on nodal finite cores lie in $G$.

Both components have dimension four, the expected dimension, and are
generically unobstructed. The smooth-cover deformation complex has
finite-length cohomology at the ramification points. At a generic
elliptic tail the attachment is unramified, and its smoothing line maps
isomorphically to the tail obstruction by $df\otimes\delta$, where
$\delta:T_pE\to H^1(E,\mathcal O_E)$ is the principal-part map.
The four remaining dimensions are those of the first-excess cover,
attachment point and pointed elliptic curve. Equation (4) leaves no
other component. A four-dimensional virtual cycle on this pure
four-dimensional stack is determined by its generic multiplicities;
these are one. Thus

$$
 [\overline{\mathcal M}_{2g+1}(C,2)^{\rm prim}]^{\rm vir}
       =[M]+\gamma_*[\mathcal U\times\overline{\mathcal M}_{1,1}].
 \tag{5}
$$

This top-cycle statement does not require global smoothness of $M$ or
separation of the components at their intersection. At a smooth
ramification-attached elliptic tail the local equation is $uv=0$:
the obstruction is divisible by the node-smoothing parameter, and its
quotient has nonzero derivative in the attachment coordinate transverse
to simple ramification.

### The two contributions outside the main component

On the source of $\gamma$, normalizing the attaching node cancels the
constant sections of the tail against the matching fiber. The normal
cokernel has the core cokernel and
$\mathbb E^\vee\otimes\operatorname{ev}^*(a\oplus b)$ as its two
successive factors, with complex orientation on the latter. This
continues to a marked node: the inserted rational component has no
$H^1$, and its matching sequence recovers the original nodal core.
With $\lambda=c_1(\mathbb E)$ and $\kappa=c_1(K_C)$,

$$
 c_2(\mathbb E^\vee\otimes\operatorname{ev}^*(a\oplus b))
       =-\lambda\operatorname{ev}^*\kappa.
$$

Here $\lambda^2=0$ on $\overline{\mathcal M}_{1,1}$ and both normal
Chern roots come from $C$. The two pushforwards are
$\int\lambda=1/24$ and $\pi_*\operatorname{ev}^*\kappa=4(g-1)$.
The ghost component therefore contributes $-(g-1)F_2(D)/6$.
All automorphism factors are retained by the universal curve and the
pointed elliptic stack; there is no additional uniform factor of two.

For a connected étale double cover $p:Y\to C$, the other summand is
$[\overline{\mathcal M}_{2g+1}(Y,1)/\operatorname{Aut}_C(p)]$.
The differential $dp$ identifies the tangential obstruction theories.
Every contracted subcurve of a degree-one map attaches just once, so
the Dolbeault kernels with coefficients $p^*a,p^*b$ are the constant
spaces $H^0(Y,p^*a),H^0(Y,p^*b)$ on the entire proper moduli stack.
Lee--Parker's constant-rank quotient formula therefore identifies the
product Euler class with the complex index Chern class, multiplied by
the determinant sign of $p^*D$. Pandharipande's degree-one formula
evaluates its series as that sign times $s(t)^{4g-4}$.
Summing with the deck weights gives the last term of (1).
[Pandharipande, Theorem 2](https://arxiv.org/abs/math/9811140)

The same argument works in every source genus and every covering degree
for the summands with degree one upstairs. Their full series is
$c_m(D)s(t)^{2m(g-1)}$. This extends the first-excess tail removal to
all contracted configurations on those summands; it is the signed spin
version of the classical étale-sector reduction.

### Extraction of the elementary exponent

Set $p=2g-2$ and write $B_{j,h,k}$ for the connected elementary
coefficients of Article III, §3. The degree-two elementary series is
$2s(t)^p s(2t)^p-s(t)^{2p}/2$. Its second-excess coefficient is

$$
 B_{2,g,2g+1}
 =\frac{2^{p+1}p(125p-34)-2p(5p-1)}{5760}.
$$

The needed degree-one coefficients are

$$
 B_{1,2g-1,2g+1}=\frac{4p(5p-1)}{5760},\qquad
 B_{1,2g,2g+1}=-\frac{p+1}{12}.
$$

For distinct spins, the earlier exponents are
$n_{2,2g-1}=-1/2-\epsilon(2^{2g}-1)/2$ and
$n_{2,2g}=\epsilon(g-1)2^{2g-3}$. Substitution in the triangular
recursion and cancellation of the étale term gives (2).
At $g=2$ the inherited coefficient is
$19/32-8/80-2/4=-1/160$, giving (3).

## A compact extension that does not compute the answer

In genus two, normalize a Poincaré line $M_0$ on
$C\times\operatorname{Pic}^2(C)$. Nonzero branch sections form
$B=\mathbb P(V)$, where $V=\pi_*M_0^2$ has rank three.
Write $P$ for the principal polarization and $\zeta=c_1(\mathcal O(1))$.
Grothendieck--Riemann--Roch gives

$$
 \operatorname{ch}(V)=3-4P,\qquad
 p_*\zeta^2=1,\quad p_*\zeta^3=4P,\quad p_*\zeta^4=8P^2.
$$

For any theta characteristic $\theta$, the line
$\theta M_0^{-1}$ has negative fiber degree. Thus
$W_\theta=R^1\pi_*(\theta M_0^{-1})$ is a rank-two bundle with
$\operatorname{ch}(W_\theta)=2+P$. On the root gerbe, where
$T^2=\mathcal O(1)$,

$$
 c_2(W_\theta\otimes T^{-1})
       =\frac{P^2}{2}-\frac{P\zeta}{2}+\frac{\zeta^2}{4}.
$$

On the finite nodal open locus this gives the spin Euler class, up to
the base spin sign. Its square integrates to zero on this compact
gerbe: the three pushforward terms are
$(1/2-1+1/2)P^2$. However, $B$ also contains nonnodal double covers.
It is not the stable-map component $M$.

There is an independent numerical check that this distinction matters.
For equal spins in genus two, the total balanced elementary coefficient
is $19/32$, the étale contribution is $3/32$, and (1) gives an
elliptic-tail contribution $1/3$. Hence the actual main integral is

$$
                         I_2(\theta,\theta)=\frac16.
$$

Therefore the zero integral of the displayed compact extension cannot
be used to prove the mixed-spin value. The difference is entirely in
the stable-map boundary and in the extension of the real obstruction
class across the cohomology-jumping loci.
