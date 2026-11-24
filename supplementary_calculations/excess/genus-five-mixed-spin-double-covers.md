# B2. The genus-five contribution of two distinct spin lines

**Theorem.** Let $C$ be a smooth curve of genus two, and let
$a,b$ be distinct even theta characteristics. For the direct sum of their
positive real spin operators, with the absolute determinant orientation,

$$
                         L_{2,5}=-\frac1{160},\qquad n_{2,5}=0.
 \tag{1}
$$

Here $L$ is the connected local Gromov--Witten degree and $n$ is the
exponent of the balanced elementary factor. The assertion concerns this
single source genus and covering degree. It neither identifies an
individual wall nor establishes the proposed all-genus partition series.

## The compact component and the common index class

Let $M$ be the reduced closure of smooth four-branched double covers in
$\overline{\mathcal M}_5(C,2)$. The two real spin cokernels restrict to
oriented rank-four bundles $V_a,V_b$ on $M$. Their existence on the full
stable-map space is an input from Lee--Parker. Let

$$
 c=c_2(-R\pi_*f^*\theta)\in H^4(M;\mathbb Q)
$$

for any theta characteristic $\theta$, where $\pi$ and $f$ are the
universal source and map. The class is independent of $\theta$: their
differences are torsion line bundles on the fixed curve, hence their
rational topological index classes agree. For even $\theta$, the
Dolbeault kernel vanishes away from the locus where an additional
section appears. On that complement, Lee--Parker's quotient and
orientation give $e(V_\theta)=c$.

We compute the difference on its compact support. No identification with
the Chern class of the spectral projective bundle is made. As shown in
Calculation B1, that extension has a different integral even for equal spins.

## The complete holomorphic-kernel locus

Write $\iota$ for the hyperelliptic involution and $W$ for the six
Weierstrass points. The kernel locus for an even $\theta$ has the
following surface components.

- $A_\theta$ consists of a finite genus-four double cover with root
  $\theta$ and an elliptic tail attached at a ramification point. Its
  branch divisor is $x+\iota(x)$. The chosen ramification point labels
  $x\in C$, and the tail is parameterized by
  $\overline{\mathcal M}_{1,1}$. The stack is a $\mu_2$-gerbe over
  this product, from the covering involution. At $x\in W$, marking
  the nodal cover inserts a rational bridge before the tail is attached.
- For each $p\in W$, $B_{\theta,p}$ consists of the étale double cover
  with character $\xi=\theta(-p)$, with an elliptic bridge joining its
  two points over $p$. Since $\theta$ is even, $\xi\ne0$.
  The bridge retains the automorphism combining sheet interchange with
  the elliptic involution exchanging its two attachments.

To see exhaustiveness, use the strata in Calculation B1. On a finite
arithmetic-genus-five double cover the anti-invariant spin line has
degree minus one. On a smooth genus-four normalization it has degree
zero, and a section forces its root to equal $\theta$. A direct node
on this normalization would require that section to vanish at two
distinct points of one unramified fiber, which it does not. On an
étale double normalization the anti-invariant line is $\theta\xi$,
of degree one. Matching across a direct node or a bridge requires its
section to vanish at the image of that node or bridge. This forces
$\theta\xi=\mathcal O_C(p)$ with $p\in W$. Two direct nodes at
different images are impossible. A tail on a one-node finite core lies
in the main component only in its marked-node limit; that limit belongs
to the intersection of the displayed surfaces. Finally, two
degree-one components have no sections since $H^0(C,\theta)=0$.

These arguments also show that the kernel is one-dimensional on the
entire locus and that the kernel loci for two distinct even theta
characteristics are disjoint. In particular an elliptic bridge supplies
one theta characteristic, determined by its character and target point.
Put

$$
                  Z_\theta=[A_\theta]+\sum_{p\in W}[B_{\theta,p}].
 \tag{2}
$$

The classes in (2) mean the Thom classes of the surface components in
their smooth ambient neighborhoods, extended to $M$.

## Deformation charts at their compact boundaries

The finite first-excess double-cover surface $F$ is the root gerbe over
the blow-up of $\operatorname{Pic}^1(C)$ at the sixteen theta
characteristics. Let $D_\theta\simeq\mathbb P^1$ be its exceptional
curve. Its coarse normal line is $\mathcal O(-1)$. The ramification-marked
cover of $D_\theta$ is $C\to\mathbb P^1$, the canonical double cover.
The exceptional curve meets each relevant nodal diagonal transversely:
it is the exceptional divisor of the blow-up, whereas the diagonal is
the strict transform of the smooth Abel curve through that theta point.

At a smooth ramification-attached tail, the tangential deformation space
has the four node-preserving parameters and one smoothing parameter
$\tau$, with one obstruction line. The obstruction is divisible by
$\tau$ and its quotient vanishes simply with $df$ at the attaching
point. Thus the two components have equation $\tau u=0$; the main
branch is smooth. Its tail divisor is Cartier and its normal is the
attaching-node smoothing line.

The only additional boundary relevant to the surface intersections is
a rational bridge carrying an elliptic tail. The one-pointed finite
cover has a smooth universal-curve chart in which its map near the
marked point $z=1$ of that rational bridge is

$$
                      f(z)=\nu+u z+v/z.
 \tag{3}
$$

The old finite node has smoothing parameter $uv$. The two bridge-node
smoothings are $u,v$, and $\nu$ moves their image in $C$. Together
with the elliptic parameter and $\tau$, these are five deformation
parameters. Normalization of the tangential deformation complex cancels
the unobstructed finite core and the constant matching terms; it leaves
exactly $H^1(E,\mathcal O_E)\otimes T_\nu C$ as obstruction. Thus a
versal analytic chart has one equation $F$ in these five parameters.

The family with $\tau=0$ exists identically. The family with $u=v=0$
also exists identically for every $\tau$: the whole bridge remains
contracted. The principal-part smoothing map is the derivative $u-v$
at $z=1$, multiplied by the nonzero elliptic connecting map. Therefore

$$
                     F=\tau(uU-vV),\qquad U(0)=V(0)\ne0,
 \tag{4}
$$

with analytic units. Replacing $u,v$ by $uU,vV$ gives a smooth main
branch $u=v$. The bridge divisor on it is $u=0$, and the tail divisor
is $\tau=0$. The root condition $D_\theta$ is transverse to the
finite nodal diagonal, so the surface equations take the form

$$
 A_\theta:\ \tau=0,\ \nu=h(u^2),\qquad
 B_{\theta,p}:\ u=0,\ \nu=0,
 \tag{5}
$$

after choosing the local target origin at $p$. Thus both surfaces are
regularly embedded and meet cleanly along
$\overline{\mathcal M}_{1,1}$. The common normal direction at this
intersection is the fixed target/root direction.

For a smooth elliptic bridge, the two smoothing lines map surjectively
onto its single elliptic obstruction; the main component is smooth
there as well. The same statement holds on the nonseparating boundary:
for an irreducible nodal elliptic bridge, or for a cycle of two rational
components with one attaching mark on each, the dualizing differential
is nonzero at each smooth attaching point. On each rational component
of the latter cycle the dualizing line has degree zero and the
differential has opposite nonzero residues at its two nodes. The
separating boundary is the rational bridge with an elliptic tail
already treated by (3)--(5), including a nodal elliptic tail. These
exhaust the stable genus-one dual graphs with two markings. The
calculations hold on orbifold charts
and retain all automorphisms. They establish the required smooth
neighborhoods of the compact kernel surfaces; global smoothness of
$M$ is not needed.

## The oriented Euler correction

We first specify the local family comparison at nodal maps. This is
stronger than the existence of the positive-mass cokernel bundle alone.
For a versal nodal family and $L=f^*\theta$, use Siebert's Banach
realization of the Dolbeault family: Corollary 4.7 identifies its kernel
and cokernel with sheaf cohomology, Corollary 4.11 identifies the family
complex with the holomorphic Cech complex, and Proposition 4.15(1),(3)
gives uniform sup-norm control on input sections and norm-continuity of
the operators in Banach-bundle trivializations.
[Siebert, version 2, §4](https://arxiv.org/abs/dg-ga/9608005v2)
Zinger explicitly records this stronger Fredholm-family statement for
flat stratified nodal families in §2.4 of his comparison theorem.
[Zinger, §2.4](https://arxiv.org/html/0807.0805v3#S2.SS4)

The real zeroth-order term $R=A(df)$ is continuous in these realizations
as well. Locally write the source plumbing as $xy=s$ and the holomorphic
map in a target coordinate as $f=F(x,y,b)$. On the half-neck
$|x|\ge\sqrt{|s|}$ its derivative is

$$
              df=(F_x-sx^{-2}F_y)\,dx.
$$

The coefficient is uniformly bounded, and the other half-neck has the
same estimate with $x,y$ exchanged. The weighted $L^p$ norm on collars
of shrinking radius tends uniformly to zero; away from those collars
the coefficients converge smoothly. Multiplication by $df$ therefore
converges in operator norm, using the uniform sup bound for sections.
Consequently $\bar\partial+tR$, including $t=0$, is a continuous local
family of Fredholm operators.

Choose a finite complex space of obstruction forms supported away from
the nodes spanning the central Dolbeault cokernel. Stabilization stays
surjective on a small product of the base and mass interval. It gives
a finite-dimensional family, whose complex at $t=0$ represents
$R\pi_*L$ by the cited Cech comparison. At positive $t$ its quotient
is the Lee--Parker bundle with its actual topology: the same
away-from-node forms give Lee--Parker's local cokernel charts
(Lemmas 8.2--8.3). The initial complex obstruction space has real
dimension six; at a positive-mass point a real four-dimensional
complement to its relation space gives a rank-four cokernel chart.
Relations among these forms in either realization
are exactly solutions of the same equation on the normalized
components, with matching values. Sobolev embedding, elliptic
regularity and removal of bounded punctures identify these solutions.
Thus the quotient charts and their overlap maps agree, rather than
merely having isomorphic individual fibers. Both determinant
orientations come from the same homotopy to $t=0$. Finite stabilizer
groups are retained by equivariant stabilization.

At a general point of either surface, a finite-dimensional reduction
of the Dolbeault operator has a one-dimensional kernel and
three-dimensional cokernel. Its variation has rank two transverse to
the kernel locus. On $A_\theta$, these directions are the root
variation transverse to $D_\theta$ and smoothing the elliptic node;
the anti-invariant section has a simple zero at the ramification point.
On $B_{\theta,p}$, moving the bridge image violates value matching
at its simple zero. For the other direction, the tangential node
smoothings satisfy $\tau_1\delta_1+\tau_2\delta_2=0$; the opposite
derivatives of the anti-invariant section give
$\tau_1\delta_1-\tau_2\delta_2=0$. This is nonzero on the surviving
tangential smoothing line.

The spin mass supplies the remaining conjugate-linear direction. Its
projection on the cokernel is nonzero by the spin norm identity. It
is orthogonal to the two complex variation directions: differentiate
the Stokes identity for $\sigma u\bar\partial u$, where
$\sigma=df$, at a holomorphic kernel section. Equivalently the Serre
pairing of either variation with $\sigma u$ vanishes. This argument
uses the dualizing sheaf and the normalization sequence also at a nodal
source. After discarding invertible factors, the transverse real model
is consequently homotopic through injective reductions to

$$
        u\longmapsto(z_1u,z_2u,\overline u):
                      \mathbb C\longrightarrow\mathbb C^3.
 \tag{6}
$$

Its quotient has the negative of the complex orientation on the first
two coordinates. The projected last coordinate has zero degree $-1$.
For the complex perturbation replacing $\overline u$ by $u$, it has
degree $+1$. These perturbations agree up to homotopy on the boundary
of a small transverse disk, where $(z_1,z_2)\ne0$. Their relative
Euler difference is $-2$ times the complex Thom class.

The difference $e(V_\theta)-c$ vanishes off the compact kernel locus.
The Thom sequence determines its degree-four supported class by the
generic coefficients of its codimension-two components. Their
intersections and further boundaries have complex codimension at
least three and support no additional degree-four class in the smooth
ambient neighborhoods just constructed. Thus

$$
                             e(V_\theta)=c-2Z_\theta.
 \tag{7}
$$

## Self-intersections and the mixed integral

For $A_\theta$ inside the ramified-tail divisor, the first normal line
is the pullback of $N_{D_\theta/F}=\mathcal O(-1)$ along the degree-two
canonical map. Its degree on $C$ is $-2$. The second normal line is
the node-smoothing line. Its first Chern class has the form
$\eta-\lambda$, where $\eta$ comes from the $C$ factor and
$\lambda=c_1(\mathbb E)$ from $\overline{\mathcal M}_{1,1}$.
At the six nodal markings the core cotangent line changes, but remains
pulled back from $C$. These changes do not contribute to the product
with the first normal class. The extra core gerbe and the elliptic
Hodge integral give

$$
              \int_M[A_\theta]^2
                 =\frac12(-2)\left(-\frac1{24}\right)=\frac1{24}.
 \tag{8}
$$

The bridge surface is the fiber at $p$ of the bridge divisor with fixed
étale character and varying image in $C$. Its normal within that
divisor is the constant line $T_pC$, so its self-intersection is zero.
At its intersection with $A_\theta$, the normal of the common curve
inside the bridge surface is the elliptic-tail smoothing line. The
excess intersection sequence cancels that line against the second
normal of $A_\theta$. The remaining line is the constant first normal
at the chosen point of $C$, and its Chern class vanishes rationally.
Hence $[A_\theta][B_{\theta,p}]=0$. Different bridge surfaces have
different images in $C$ and are disjoint. In particular

$$
                \int_M Z_\theta^2=\frac1{24},\qquad
                Z_aZ_b=0\quad(a\ne b).
 \tag{9}
$$

For each equal-spin pair, Calculation B1 gives
$\int_M e(V_\theta)^2=1/6$ from the elementary total, the complete
étale contribution, and the separate elliptic-tail component.
Set $I=\int_M c^2$ and $J_\theta=\int_M cZ_\theta$. Equations
(7)--(9) imply

$$
                 \frac16=I-4J_\theta+\frac16,
                 \qquad J_\theta=I/4.
$$

For distinct even $a,b$ it follows that

$$
       \int_M e(V_a)e(V_b)=I-2J_a-2J_b+4\int_M Z_aZ_b=0.
$$

Equation (3) of Calculation B1 now gives (1). This cancellation does not
require the separate value of the common complex-index integral $I$.

There is a separate restriction check on the intersections. On
$A_\theta$, let $K=T^{-1}$ be the anti-invariant kernel line of the
finite core, where $T^2=\mathcal O_{D_\theta}(1)$. Its pullback to
the ramification parameter curve $C$ has degree $-1$. The core's
complex negative index is $2K-K=K$, and the tail supplies the line
$\operatorname{ev}^*\theta\otimes\mathbb E^\vee$. Hence

$$
            \int_{A_\theta}c
              =\frac12(-1)\left(-\frac1{24}\right)=\frac1{48}.
$$

On $B_{\theta,p}$ the complex index is the matching line $\theta_p$
plus $\theta_p\otimes\mathbb E^\vee$, modulo constant core index
spaces of equal dimension. Its second Chern class is zero. Thus
$J_\theta=1/48$ and $I=1/12$. Lee--Parker's quotient gives the real
first-excess core Euler class $-c_1(2K-\overline K)=-3c_1(K)$:
the mass conjugates the one-dimensional kernel and reverses its
determinant orientation. This class has degree $3$ on the ramification
parameter curve, so its product
with the tail line integrates to $-1/16$ on $A_\theta$. This agrees
with $1/48-2/24$ from (7)--(9). It checks both the orientation and
the extra core stack factor independently of the mixed-product
cancellation.

These restrictions also give a second proof of the mixed vanishing
which does not use the transverse model (6). The support argument
alone writes $e(V_\theta)-c$ as
$\alpha_\theta[A_\theta]+\sum_p\beta_{\theta,p}[B_{\theta,p}]$.
Restriction to $A_\theta$ and (8)--(9) give
$\alpha_\theta/24=-1/16-1/48=-1/12$, hence
$\alpha_\theta=-2$. All terms involving a bridge coefficient
$\beta_{\theta,p}$ vanish in the required pairings: $c$ integrates
to zero there, and all its self-intersections and intersections with
the other kernel surfaces are zero. The equal-spin calibration then
gives $I=1/12$, and the distinct-spin product is
$1/12-2/48-2/48=0$. This proof uses only the constant-rank
Lee--Parker formula on the compact kernel surfaces, their regular
neighborhoods, and the intersection calculation.
