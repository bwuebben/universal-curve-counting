# B4. Four branch points and the spin obstruction in genus two

test. The target curve $C$ has genus two. Unless specified otherwise,
$\theta$ is an even theta characteristic; thus $H^0(C,\theta)=0$.
The two spin lines in the mixed calculation can both be chosen even and
distinct. The results below distinguish a compact spectral parameter space
from the stable-map compactification. In particular, the first zero computed
below is not the local Gromov--Witten invariant.

## The finite-cover calculation

Let $J=\operatorname{Pic}^2(C)$, let $\mathcal M$ be the Poincaré line on
$C\times J$ normalized at a fixed point of $C$, and let $P$ be the principal
polarization class on $J$. Put

$$
 V=\pi_*\mathcal M^2,\qquad
 E_\theta=R^1\pi_*(\theta\mathcal M^{-1}),\qquad
 B=\mathbb P_J(V),\qquad \zeta=c_1(\mathcal O_B(1)).
$$

Here projectivization parametrizes lines. The ranks are three and two:
every degree-four line on $C$ has three sections and no first cohomology,
whereas $\theta M^{-1}$ has degree minus one and no sections. The point
$(M,[s])$ of $B$ gives the double-cover algebra
$\mathcal O_C\oplus M^{-1}$ with equation $w^2=s$.

The stack of these spectral curves is the gerbe on which a line $T$ with
$T^2=\mathcal O_B(1)$ exists. The universal root is $\mathcal M\otimes T$.
On the open subset $B^{\mathrm{nod}}$ where the zeros of $s$ have
multiplicity at most two, the source is nodal. The finite-cover normal
cokernel is, as an oriented real bundle, the underlying real bundle of
$E_\theta\otimes T^{-1}$. Indeed, the invariant Dolbeault cohomology is
zero, the anti-invariant kernel vanishes, and the positive spin mass gives
the complex cokernel projection of Lee--Parker, Proposition 12.2.
This statement includes finite nodal domains by the normalization sequence.

Grothendieck--Riemann--Roch gives

$$
 \operatorname{ch}(V)=3-4P,\qquad
 \operatorname{ch}(E_\theta)=2+P.
$$

Consequently the Euler class on this open set is the restriction of

$$
 c=\frac{P^2}{2}-\frac{P\zeta}{2}+\frac{\zeta^2}{4}.
 \tag{1}
$$

This class is independent of $\theta$ in rational cohomology. The
projective-bundle formula gives

$$
 p_*\zeta^2=1,\qquad p_*\zeta^3=4P,\qquad
 p_*\zeta^4=8P^2.
$$

Therefore

$$
 \int_B c^2=
 \int_J\left(\frac12P^2-P^2+\frac12P^2\right)=0.
 \tag{2}
$$

The gerbe multiplies this zero by $1/2$. Equation (2) is an intersection
calculation on the spectral compactification. At a triple zero the source
has a cusp, and at a quadruple zero it has a tacnode. These are not stable
nodal domains. The corresponding stable limits contain contracted elliptic
tails or bridges. Thus (2) cannot be substituted for an Euler integral on
$\overline{\mathcal M}_5(C,2)$.

## The holomorphic kernel on the main component

Consider the closure of the smooth-source double covers inside the
stable-map stack, called the main component here. For even $\theta$, the
holomorphic pullback acquires a section only in the following two kinds of
codimension-two loci.

First, a degree-two genus-four finite cover with root $M_1=\theta$ has one
anti-invariant section. Attach an elliptic tail at a ramification point.
The section extends over the tail by its attaching value, which is zero.
The resulting surface is denoted $A_\theta$. Its compact parameter curve
for the pointed finite cover is $C$: the branch divisor is $x+\iota(x)$,
where $\iota$ is the hyperelliptic involution, and the chosen ramification
point lies above $x$. The tail parameter is
$\overline{\mathcal M}_{1,1}$. There is an additional generic stack factor
$1/2$ from the double-cover involution, which fixes the ramification point.
At a Weierstrass point the pointed finite cover is stabilized by inserting
a rational bridge at its node before attaching the elliptic tail.

Second, let $p$ be one of the six Weierstrass points and put
$\xi=\theta(-p)\in\operatorname{Pic}^0(C)[2]$. This is nontrivial because
$\theta$ is even. Take the étale double cover with character $\xi$ and
join its two points over $p$ by a contracted two-pointed elliptic curve.
The anti-invariant section corresponds to the unique section of
$\theta\xi=\mathcal O_C(p)$; it vanishes at both attaching points.
These maps form a surface $B_{\theta,p}$. It meets $A_\theta$ along the
boundary where the two attaching points lie on a rational component with
an elliptic tail. This common curve has parameter
$\overline{\mathcal M}_{1,1}$.
The bridge retains the simultaneous core deck transformation and
elliptic involution exchanging the two attaching points. An elementary
quotient of a sheet-labelled product can lose this stabilizer, so it
should not be substituted for the bridge stack without a separate
comparison. The vanishing intersections below use its constant target
normal line and do not require such a global presentation.

The elementary kernel argument is useful independently of intersection
theory. With four simple branch points the anti-invariant line has degree
minus one. With two branch points it has degree zero, and a section forces
$M_1=\theta$. On an étale normalization it is a theta characteristic;
matching the two sheets across a bridge forces its zero to be the bridge
image. A genus-two effective theta characteristic is $\mathcal O_C(p)$
for a Weierstrass point $p$. Two degree-one positive components have no
sections because $H^0(C,\theta)=0$. These observations also show that the
kernel loci for two distinct even theta characteristics are disjoint.
The bridge surfaces cannot be omitted from the compact kernel locus.

## The tail--bridge corner

The following calculation supplies a local model for the intersection
$A_\theta\cap B_{\theta,p}$. It uses the usual analytic Kuranishi
presentation of the deformation complex of a stable map. At the central
map the nonconstant normalization is étale. The contracted subcurve is a
rational bridge with an elliptic tail. After matching evaluations, its
only tangential obstruction is the one-dimensional space
$H^1(E,\mathcal O_E)\otimes T_pC$.

In a chart for a one-pointed finite nodal cover, write the map near the
inserted rational bridge as

$$
 f(\xi)=\nu+a\xi+b/\xi,
 \qquad \text{marked point }\xi=1.
 \tag{3}
$$

The two attaching-node smoothing parameters are $a,b$, and the smoothing
parameter of the finite node before inserting the marked bridge is $ab$.
The parameter $\nu$ moves its image. The derivative at the marked point
in (3) is $a-b$. Let $\tau$ smooth the node joining the elliptic tail.
Together with an orbifold parameter $j$ for the pointed elliptic curve,
these give five tangent coordinates and one obstruction equation $F$.

The whole elliptic-tail component exists for $\tau=0$, so $F$ is divisible
by $\tau$. When $a=b=0$, the rational bridge and elliptic tail remain
contracted for every $\tau$; therefore $F/\tau$ belongs to $(a,b)$.
The principal-part connecting homomorphism on an elliptic curve is an
isomorphism, also on its irreducible nodal limit. Its leading term is
$a-b$. Thus

$$
 F=\tau(aU-bV),\qquad U(0)=V(0)\ne0,
 \tag{4}
$$

with $U,V$ analytic units. Replacing $a,b$ by $aU,bV$ puts the equation
in the form $\tau(a-b)=0$. In particular the main branch is smooth in
this chart; its elliptic-tail and elliptic-bridge divisors have equations
$\tau=0$ and $a=b=0$, respectively.

The exceptional curve $D_\theta$ on the first-excess finite-cover surface
meets the nodal diagonal transversely. In branch-center and node-smoothing
coordinates it has equation $\nu=h(ab)$ for an analytic $h$ with
$h(0)=0$. Hence, on the main branch, the two kernel surfaces have local
equations

$$
 A_\theta:\quad \tau=0,\quad\nu=h(a^2),
 \qquad
 B_{\theta,p}:\quad a=0,\quad\nu=0.
 \tag{5}
$$

They are smooth and meet cleanly along the elliptic-modulus curve. The
argument for (4) needs both exact vanishing loci, not only a first-order
smoothing calculation: the latter alone would not justify the
normal-crossing equation. A complete assembled proof should spell out
the Kuranishi construction and verify its compatibility with the
stack charts at every boundary point.

## The local Euler correction

There is a simple finite-dimensional model for a transverse
codimension-two jump of a one-dimensional Dolbeault kernel. After
discarding invertible summands it is

$$
 D_z(u)=(z_1u,z_2u,\overline u):\mathbb C\longrightarrow\mathbb C^3.
 \tag{6}
$$

The real cokernel has rank four. Its quotient identification with the
first two complex coordinates has the negative of their complex
orientation, since conjugation reverses the orientation of the
one-dimensional kernel. The section induced by the third coordinate has
local zero degree minus one. Replacing the last coordinate of (6) by
$\varepsilon u$, $\varepsilon\ne0$, gives the corresponding complex
index representative, for which this zero degree is plus one. Their
difference is therefore minus two times the complex orientation of the
kernel locus.

The two required transverse directions are visible geometrically. On
$A_\theta$ they are variation of the root away from $\theta$ and smoothing
the elliptic attachment: the anti-invariant section has a simple zero at
the ramification point. On a generic bridge surface, moving its image
away from $p$ violates matching to first order. If $\tau_1,\tau_2$ are
the node smoothings and $\delta_1,\delta_2$ the elliptic principal-part
maps, the tangential equation is
$\tau_1\delta_1+\tau_2\delta_2=0$. The anti-invariant section has opposite
derivatives on the two sheets, so its equation is
$\tau_1\delta_1-\tau_2\delta_2=0$. Restricted to the tangential smoothing
line, the latter is nonzero. Both surface components therefore have
generic multiplicity one in the kernel degeneracy cycle.

If these charts are used to identify the global real Euler correction,
the resulting formula is

$$
 e(V_\theta)=c_2(-R\pi_*f^*\theta)-2Z_\theta,
 \qquad
 Z_\theta=[A_\theta]+\sum_{p\ {
m Weierstrass}}[B_{\theta,p}].
 \tag{7}
$$

For a full proof one must use the compact-support/Thom sequence to
exclude any additional codimension-two contribution. Equation (6)
determines the generic coefficients; strata of complex codimension at
least three cannot support a degree-four correction on a smooth ambient
chart. The common rational Chern class in (7) is independent of
$\theta$, by Grothendieck--Riemann--Roch and the torsion difference
between theta characteristics. This Chern class is on the stable-map
main component, and should not be confused with the particular
extension (1) on $B$.

## Intersections of the compact kernel surfaces

The expected intersection calculation has a useful robustness property.
The normal line to $A_\theta$ inside the ramified-tail divisor is the
pullback of $N_{D_\theta/B_2}=\mathcal O_{\mathbb P^1}(-1)$. The
ramification-marking map $C\to\mathbb P^1$ has degree two, so this normal
line has degree minus two. The other normal line is the attaching-node
smoothing line. Its dependence on the elliptic modulus is
$-\lambda$, where $\lambda=c_1(\mathbb E)$ on
$\overline{\mathcal M}_{1,1}$. Any corrections at the six nodal marked
covers are pulled back from the $C$ factor and disappear when multiplied
by the first normal class. Thus, provided the regular embeddings and
normal identifications above extend on the compact stack,

$$
 A_\theta^2=\frac12\cdot(-2)\cdot
       \left(-\frac1{24}\right)=\frac1{24}.
 \tag{8}
$$

The factor $1/2$ is the core involution; $1/24$ already includes the
elliptic stack stabilizer.

Each $B_{\theta,p}$ is a fiber at the fixed image $p$ of the
elliptic-bridge divisor. Its normal line inside that divisor is the
constant line $T_pC$, so $B_{\theta,p}^2=0$. For the clean intersection
with $A_\theta$, the normal of the common curve inside the bridge surface
is the elliptic-tail smoothing line. It cancels the same line in
$N_{A_\theta}$, leaving a constant excess line. Hence
$A_\theta B_{\theta,p}=0$. Distinct bridge surfaces have distinct bridge
images and are disjoint. Consequently (8) gives

$$
 Z_\theta^2=\frac1{24},\qquad
 Z_aZ_b=0\quad(a\ne b).
 \tag{9}
$$

The zero mixed intersections are substantive: merely dropping the six
bridge surfaces from $Z_\theta$ would leave the degeneracy analysis
incorrect, even if it accidentally produced the same number.

Equations (7)--(9), together with a justified decomposition of the
primitive virtual class and its equal-spin calibration, would imply
vanishing of the mixed main-component integral. At this stage the
assembled virtual-class argument and a fresh independent check of the
compact charts remain necessary. No value of $n_{2,5}$ is asserted by
this investigation alone.

For reference, the numerical calibrations furnished by the elementary
hook formula and degree-one contribution are

$$
 B_{2,2,5}=\frac{19}{32},\qquad
 B_{1,3,5}=\frac1{80},\qquad
 L^{\rm elem,prim}_{2,5}=\frac{19}{32}
       -\frac{15}{2}\frac1{80}=\frac12.
$$

On the elliptic-tail component over the first-excess finite-cover
surface, normalization identifies the normal Euler product with

$$
 e(V_a^{\rm core})e(V_b^{\rm core})
    (\operatorname{ev}^*c_1(a)-\lambda)
    (\operatorname{ev}^*c_1(b)-\lambda).
$$

The square of a target class and $\lambda^2$ vanish. Pushforward along
the degree-two universal cover, then along
$\overline{\mathcal M}_{1,1}$, gives $-F_2(D_{ab})/6$ in genus two.
Thus the proposed ordinary-fundamental-class contribution is $1/3$ for
equal spins and zero for distinct spins. The equal-spin main integral
would then be $1/2-1/3=1/6$. With (7)--(9), comparing equal and distinct
spin products subtracts $4Z_\theta^2=1/6$, giving a zero mixed main
integral. Finally the signed étale component has coefficient
$(-1/2)(1/80)=-1/160$, and triangular extraction gives
$n_{2,5}=L_{2,5}+1/160$. These identities explain why the proposed
intersection calculation addresses precisely the desired integer;
their geometric use still requires the component and virtual-class
arguments stated above.

## Primary-source comparison

- Lee--Parker, *An obstruction bundle relating Gromov--Witten invariants
  of curves and Kähler surfaces*, §12, gives the constant-rank quotient
  and orientation used above. Its construction gives the stable real
  obstruction bundle; it does not give the mixed-spin integral here.
  [Primary PDF](https://users.math.msu.edu/users/parker/LeeParkerObstructionBundle.pdf).
- Kiem--Li, *Low degree GW invariants of spin surfaces*, Theorem 1.1
  and §4, treats the analogous obstruction-theoretic comparison and
  analyzes corrections where Dolbeault cohomology fails to be locally
  free. It concerns one spin class. It is not a theorem about products
  of two distinct spin Euler classes. An earlier manuscript contains
  the related calculation in §4:
  [arXiv:0707.2986](https://arxiv.org/abs/0707.2986).
- Nesterov--Schimpf, *Admissible covers and stable maps*, Theorem 2.7
  and §7, gives a cycle-valued comparison with explicit contracted-curve
  corrections, including a fixed-target formulation. The spin Euler
  classes would require a compatible pullback comparison before applying
  that theorem here.
  [arXiv:2505.03487](https://arxiv.org/abs/2505.03487).
