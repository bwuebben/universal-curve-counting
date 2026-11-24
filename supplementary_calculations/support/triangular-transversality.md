# A9. Triangular graph equations and degree-two transversality

**Conclusion.** A triangular family has more perturbation freedom than the
split family in the preceding investigation, while retaining its genus
bound. In base genus two, this freedom suffices for the local, degree-at-most-two
version of the four conditions in Doan--Ionel--Walpuski Definition 3.13.
The argument uses a scalar normal line, rather than an assertion that a
generic ambient perturbation preserves the triangular form. It also gives
the smooth-domain double-cover transversality tests used by
Bai--Swaminathan. These are local analytic assertions. Applying the
cluster comparison still requires its ambient compactness and localization
hypotheses and the identification of the endpoint contributions.

Throughout, all dimensions and codimensions are real. Claims of residual
genericity refer to sufficiently differentiable Banach completions and
then to smooth parameters by the usual countable local perturbation
argument. The proofs below specify the surjectivity and finite-jet tests
needed for that argument.

## 1. Parameter space and genus bound

Let $C$ be a compact Riemann surface, let $L_1,L_2$ be complex line
bundles of degree $g(C)-1$, and fix a real Cauchy--Riemann operator $D_2$
on $L_2$ such that

$$
                  \ker h^*D_2=0                         \tag{1}
$$

for every nonconstant holomorphic map $h$ from a smooth compact connected
Riemann surface to $C$. A positive-mass spin operator has this property;
the spin-square/Stokes proof is recalled in
Article III, Lemma 14.1.
The same conclusion holds for connected stable sources: a kernel
section vanishes on every positive-degree component, is constant on
every contracted component, and vanishes everywhere by matching and
connectedness.

Fix a complex-linear reference $\bar\partial_{L_1}$. Consider smooth maps
over $C$

$$
 A:L_1\longrightarrow\Lambda^{0,1}T^*C\otimes L_1,
 \qquad F:L_1\longrightarrow\Lambda^{0,1}T^*C\otimes L_2,
$$

and the graph equations

$$
 \bar\partial_{f^*L_1}u_1+f^*A(u_1)=0,
 \qquad f^*D_2u_2+f^*F(u_1)=0.                          \tag{2}
$$

The form factors in $A,F,D_2$ are pulled back using $df$. The associated
almost-complex structure has holomorphic projections to $C$ and to the
almost-complex surface $L_1$ defined by the first equation. All fibers
are standard complex vector spaces. The parameters may be supported in
a prescribed disk bundle; the assertions concern curves contained in
the region where (2) holds.

**Proposition 1.** The conclusion of Calculation A1 holds for (2), with
$d=\deg L_1$. Moreover, every reduced simple graph in (2) has a simple
projection to $L_1$; distinct full graph components have distinct scalar
images.

**Proof.** Use the normalized off-diagonal fiber product $Z$ from
[Calculation A1](genus-bound-for-split-multisections.md). On a connected
component of $Z$, the difference $v_1$ of the first coordinates satisfies
a homogeneous scalar real Cauchy--Riemann equation, by the integral
formula for the fiber derivative of $A$. If $v_1$ vanishes identically,
the difference of the two second equations is $h^*D_2v_2=0$.
Equation (1) gives $v_2=0$. This contradicts reducedness and generic
injectivity of the full graph. Thus $v_1$ is nonzero on every component
of $Z$. Its zero count is $(\deg h)\deg L_1$.

The forced ramification and nodal zeros are exactly those counted in
Calculation A1, so

$$
 R(f)+2b\le m(m-1)\deg L_1,
 \qquad
 p_a\le1+m(g-1)+\frac{m(m-1)}2\deg L_1.                 \tag{3}
$$

The same argument excludes a generic identification in the scalar
projection, both within a component and between components. $\square$

The scalar zero theorem used here is the similarity principle and its
first-Chern-class consequence; see
[Wendl, *Lectures on Symplectic Field Theory*, Theorem 2.5.3](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf).
The coupling allowed in (2) is essential: arbitrary dependence of the
first equation on $u_2$ does not satisfy this proof.

### The special geometry in genus two

Now fix $g=2$. For a simple degree-two graph with smooth domain, let
$\tau$ be the deck involution and put $v_1=u_1-u_1\tau$. It is nonzero
and has precisely two zeros, counted with multiplicity. Consequently:

* A branched graph has source genus four and exactly two ramification
  points. These exhaust the zeros of $v_1$, both simply. Thus its scalar
  projection is injective, and its derivative at each ramification point
  is nonzero. It is an embedding.
* An étale graph has source genus three. The zeros of $v_1$ are one
  exchanged pair and are simple. Its scalar projection is an immersion
  with exactly one transverse double point. Its full graph is embedded
  precisely when the two second-coordinate values at that pair differ.
* No simple degree-two graph with smooth domain has genus at least five.

At a ramification point choose a coordinate with $f(w)=w^2$. The equation
at $w=0$ gives $\partial_{\bar w}u_1(0)=0$, and the simple zero of
$u_1(w)-u_1(-w)$ gives $\partial_wu_1(0)\ne0$. This proves the immersion
claim without a transversality assumption. At an étale coincidence the
same calculation in base coordinates gives distinct tangent lines.
Degree-one graphs are sections and hence embedded.

## 2. Linear double-cover kernels

Suppose $A(z,0)=F(z,0)=0$. The normal operator of the zero section has
the form

$$
 T=\begin{pmatrix}D_1&0\\ B&D_2\end{pmatrix},\qquad
 D_1=\bar\partial_{L_1}+d_{\rm fib}A|_0,
 \quad B=d_{\rm fib}F|_0.                              \tag{4}
$$

Let $f:Y\to C$ be a connected smooth double cover of genus
$h=2g-1+r$, so there are $2r$ branch points. Write $T_f^-$ for the
operator on the anti-invariant sections. Its index is $-4r$.
The invariant summand has index zero and the index of the full pullback
is $-4r$, as follows by deformation to a complex-linear operator and
Riemann--Roch.

If $\sigma=(s,v)\in\ker T_f^-$ is nonzero, then $s\ne0$, by (1).
For arbitrary first-column variations,

$$
      \dot T_f\sigma=(f^*\dot A\,s,f^*\dot B\,s).       \tag{5}
$$

On an unbranched disk where $s\ne0$, the two real bundle maps
$\dot A,\dot B$ can prescribe both output coordinates arbitrarily.
The prescription on the other sheet is forced by oddness and is
consistent with descent to the base.

**Proposition 2.** In the triangular operator space times the smooth
double-cover space:

1. The locus $\dim\ker T_f^-=1$ is a submanifold of codimension
   $4r+1$.
2. Near a point of nullity $k\ge2$, the nullity-$k$ locus is contained
   in a submanifold of codimension at least $2(k+4r)$.

**Proof.** At nullity one, a cokernel functional annihilating (5) is
represented by an anti-invariant adjoint solution. Arbitrary local
prescriptions in (5) force that solution to vanish on a disk. Unique
continuation makes it zero globally. Thus the normal derivative onto
$\operatorname{Hom}(\ker T_f^-,\operatorname{coker}T_f^-)$ is onto.

For the second assertion, choose any independent pair
$\sigma_1,\sigma_2$ in the kernel and write their first coordinates
as $s_1,s_2$. They remain linearly independent, because the second
coordinate is uniquely determined by the first. They are pointwise
linearly independent over $\mathbb R$ on some disk. Indeed, if
$s_2=a s_1$ with real $a$ on an open set where $s_1\ne0$, the common
scalar equation gives $\bar\partial a=0$; hence $a$ is constant, and
unique continuation would make the sections globally dependent.
If they were pointwise dependent everywhere this argument would apply
on such an open set.

The dual of the first-column variation map on this pair is

$$
 (\eta_1,\eta_2)\longmapsto
                 s_1\otimes_{\mathbb R}\eta_1+
                 s_2\otimes_{\mathbb R}\eta_2.           \tag{6}
$$

If (6) vanishes, pointwise independence forces both adjoint solutions
to vanish on the disk, hence everywhere. The derivative restricted to
the two chosen kernel columns is therefore onto
$\operatorname{Hom}_{\mathbb R}(\mathbb R^2,
\operatorname{coker}T_f^-)$, of dimension $2(k+4r)$.
Finite-dimensional Fredholm reduction and the implicit function theorem
give the stated local codimension bound. This does not claim that every
higher-kernel stratum is a smooth manifold of its maximal expected
codimension. $\square$

The cover space has dimension $4r$. Thus generic paths in the triangular
operator space meet only simple kernels, and the map

$$
 T_f\mathcal M_h(C,2)\oplus\mathbb R_t
 \longrightarrow
 \operatorname{Hom}(\ker T_f^-,\operatorname{coker}T_f^-)
                                                               \tag{7}
$$

is an isomorphism at each intersection. Higher-kernel loci have
codimension at least $4+4r$ after allowing the cover to vary, so a
one-parameter family avoids them. At a base-invertible operator the
anti-invariant spaces in (7) are the full kernel and cokernel spaces.

These statements verify the infinitesimal condition in
[Bai--Swaminathan, Definition 6.5](https://arxiv.org/html/2106.01206v4).
They do not by themselves verify that definition's separate condition
on the compactified stable-map space. The distinction is material at
an étale wall with attached contracted components.

### Removing the étale double point

For $g=2$ an étale odd scalar kernel $s$ descends to a section of
$L_1\otimes\alpha$ with one simple zero $p$, where $\alpha$ is the
nontrivial sign line of the cover. The operator $D_2\otimes\alpha$ is
invertible. The second coordinate is

$$
                 v=-(D_2\otimes\alpha)^{-1}(B s).        \tag{8}
$$

Varying $B$ can prescribe $v(p)$ arbitrarily while leaving $D_1,s$
fixed. Here is a direct proof, including the constraint $s(p)=0$.
Choose a local solution $\xi$ of $(D_2\otimes\alpha)\xi=0$ with any
specified value at $p$; local solutions with arbitrary value follow
from [Wendl, Lemma 2.5.1](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf#page=63)
and elliptic regularity. Multiply by a cutoff equal to one near $p$. The result
has $(D_2\otimes\alpha)\xi$ supported in an annulus where $s\ne0$.
There choose a real-linear map $\dot B$ with

$$
            \dot B s=-(D_2\otimes\alpha)\xi,
$$

and set it equal to zero elsewhere. Then $(s,v+\epsilon\xi)$ is an
exact kernel vector for $(D_1,B+\epsilon\dot B,D_2)$. In particular,
the condition $v(p)=0$ has two additional independent real equations
inside the simple-kernel wall. Its projection to operator space has
codimension three, so generic paths avoid it. Branched kernel graphs
are already embedded by §1.

In genus two there are only fifteen connected étale covers. The
two-branch Hurwitz space compactifies by the spectral double covers
over branch divisors in $\operatorname{Sym}^2C$. When the branch points
coincide, the finite curve has one node. Its normalization is either
an étale double cover or two copies of $C$. A nonzero odd kernel on
the nodal limit would therefore give either an étale kernel vanishing
at its exchanged pair or a base kernel. Along a base-invertible path
whose étale kernel graphs are embedded, neither occurs. Normalized
kernel graphs, scaled into a uniformly tame small tube, satisfy
fixed-genus Gromov compactness. Thus simple branched kernel crossings
cannot escape through this Hurwitz boundary; the smooth kernel
crossings are finite along a compact generic path.

At a genus-four branched crossing with no simultaneous étale kernel,
the remaining stable degree-two genus-four strata have no kernel:
their positive-degree normalizations have either degree one or are
étale. A genus-three étale crossing satisfies the simple-wall
condition in genus three. It does not satisfy it in every larger
source genus, because contracted tails inherit its kernel.

## 3. The normal line of a scalar image

The following observation is the additional ingredient needed to
control simple curves away from the zero section.

Let $u:Y\to L_1\oplus L_2$ be a simple solution of (2) of degree one
or two over a genus-two base. Let $\bar u:Y\to L_1$ be its scalar
projection. By §1 both maps are immersions, even before generic
embeddedness is imposed. Let $N_u$ and $N_{\bar u}$ be their normal
bundles on $Y$; $N_u$ has complex rank two and $N_{\bar u}$ has
rank one. The holomorphic projection induces an exact sequence of
normal Cauchy--Riemann operators

$$
 0\longrightarrow(f^*L_2,f^*D_2)
 \longrightarrow(N_u,D_u^N)
 \longrightarrow(N_{\bar u},Q_{\bar u})
 \longrightarrow0.                                    \tag{9}
$$

The operator $Q_{\bar u}$ is the normal real Cauchy--Riemann operator
of the immersed curve in the almost-complex surface $L_1$.
The exactness of the bundles follows because $d\bar u$ is injective;
compatibility of the operators follows by differentiating the
holomorphic projection. The operator on the vertical subbundle is
exactly $f^*D_2$, since the second equation is linear in $u_2$.
After a smooth splitting, (9) is a triangular system; its off-diagonal
term has order zero.

In particular,

$$
               \ker D_u^N\hookrightarrow\ker Q_{\bar u}. \tag{10}
$$

The index of $D_u^N$ is zero: $c_1(T(L_1\oplus L_2))[C]
=2-2g+\deg L_1+\deg L_2=0$, and the unmarked curve index in
complex dimension three is twice this Chern number.

Every nonzero kernel element therefore has a nonzero scalar normal
component. Any two independent kernel elements project to two
independent solutions of the same scalar operator. The independence
argument used in (6) applies unchanged to these normal sections.

## 4. Universal simple-curve transversality with prescribed core jets

Fix a smooth path of operators (4) with $D_{1,t}$ invertible on $C$.
Require the nonlinear perturbations to preserve the zero section and
its first fiber jet:

$$
 A_t(z,0)=F_t(z,0)=0,
 \quad d_{\rm fib}A_t|_0=D_{1,t}-\bar\partial_{L_1},
 \quad d_{\rm fib}F_t|_0=B_t.                           \tag{11}
$$

One may keep the complete germs fixed in collars of the two endpoint
parameters. Perturbations in the interior may have arbitrary higher
fiber jets and arbitrary support away from the zero section.

**Proposition 3 (local isolation at degree two).** Suppose the graph
structures are realized in a common tame tube in a closed symplectic
six-manifold. Assume the fixed endpoint collars contain only the
regular zero-section support; constant positive-spin endpoint germs
satisfy this assumption. There is a residual set of arbitrarily small smooth
triangular perturbations satisfying (11) and the endpoint constraints
such that, among simple closed curves contained in the tube and of
degree at most two:

1. every curve is embedded;
2. distinct curves are disjoint, including disjointness from the zero
   section;
3. for each fixed parameter the moduli space of simple curves is
   discrete;
4. away from the zero section the moduli over the parameter interval
   is a smooth one-dimensional manifold, and its projection to the
   interval has only nondegenerate critical points.

The zero section is regular throughout. All the curves in this
statement have index zero, so this is precisely the four local
conditions of Doan--Ionel--Walpuski Definition 3.13 at this degree
cutoff.

**Proof.** The relevant perturbation tests are as follows.

**(a) Free values and Fredholm surjectivity.** Every noncore curve has
an unbranched disk whose scalar image is injective and disjoint from
the scalar images of any prescribed finitely many other distinct
curves and from the zero section. To see the last assertions, scalar
images are simple, and distinct scalar images cannot agree on an open
set: if they agreed globally, their second-coordinate difference
would be a kernel of $f^*D_2$. Distinct scalar images have isolated
intersections. The same applies to intersection with the zero section.
Choose the disk away from these intersections and from branch points.

On this disk the vertical directions identify with the full normal
bundle of the graph, and arbitrary values of $(\dot A,\dot F)$
prescribe arbitrary normal $(0,1)$-forms. These values extend to smooth
perturbations supported in a neighborhood of the scalar disk in $L_1$.
They satisfy all of (11). An adjoint cokernel element annihilating
these perturbations vanishes on the disk and hence everywhere by
unique continuation. The universal normal Fredholm section is
therefore transverse. At a fixed $t$, the projection of the universal
curve space to the space of static perturbations has Fredholm index
zero. After including $t\in[0,1]$, the projection of the universal
path-curve space to the space of perturbation paths has Fredholm index
one. A regular chosen path consequently has one-dimensional curve
moduli; the subsequent projection of this one-manifold to $[0,1]$ is
the map whose critical points are considered below.

**(b) Kernel strata.** Hold the current curve fixed by requiring the
perturbation values to vanish along its scalar image. The first
derivative transverse to that scalar image remains an arbitrary real
map from the scalar normal line to the two normal output coordinates,
supported in the chosen disk. It acts on $\ker D_u^N$ through the
injection (10). The proof of Proposition 2 now shows that the
one-dimensional-kernel locus has codimension one in the universal
simple-curve space. At nullity $k\ge2$, the derivative on any two
kernel columns is onto a space of dimension $2k$, so the nullity-$k$
locus is locally contained in a submanifold of codimension at least
four. These estimates suffice for avoidance by a one-parameter family.

**(c) Nondegenerate quadratic term.** At nullity one choose a kernel
generator $\xi$ and an adjoint cokernel generator $\eta$. Its scalar
normal component $\bar\xi$ is nonzero on an open set. Vary the second
normal derivative of $(A,F)$ on a smaller disk there, while holding
their values and first derivatives fixed along the scalar curve.
Such a perturbation changes the scalar quadratic Kuranishi coefficient
by a nonzero constant multiple of

$$
   \int_Y\left\langle\eta,
          \dot P^{(2)}(\bar\xi,\bar\xi)\right\rangle.    \tag{12}
$$

Because $\eta$ satisfies unique continuation, the disk can be chosen
where both factors are nonzero. A real quadratic map evaluated on a
fixed nonzero real vector can take any prescribed output value.
Thus (12) can be varied arbitrarily, without changing the curve or
its linearized operator. Vanishing of this coefficient imposes one
additional independent equation inside the nullity-one locus.

Relative parametric transversality therefore makes the curve moduli
over the path a one-manifold, avoids nullity at least two and the
vanishing-quadratic locus, and meets the nullity-one locus transversely.
The local equation for the projection is $t-t_0=\pm x^2$ at a
critical point. In particular its fiber at every fixed $t$ is locally
discrete. At the zero section regularity follows directly from (11)
and the invertibility of both diagonal base operators.

**(d) Evaluation constraints.** For an immersed simple graph, the
universal normal equation together with evaluation at two distinct
marked points is transverse to the coincidence condition. The same
holds for a pair of distinct curves, and for intersection with the
fixed zero section. We first prove a bounded linear statement at
fixed smooth maps and fixed marked points, and only then allow the
marks to move in higher-regularity charts.

Fix $p>2$ and distinct marks $z_1,z_2$. On $W^{1,p}(N_u)$, evaluation
at these fixed points is bounded. Consider the augmented linear map

$$
 (\xi,\dot P)\longmapsto
 \left(D_u^N\xi+\dot P_u,\xi(z_1),\xi(z_2)\right)
 \in L^p(\Lambda^{0,1}\otimes N_u)
       \oplus N_{u,z_1}\oplus N_{u,z_2},
$$

where $\dot P=(\dot A,\dot F)$ ranges over smooth permitted
perturbations. The operator on $\xi$ alone is Fredholm, since the
evaluation terms are finite-dimensional. Its image has finite
codimension, so adding parameter directions still gives a closed
image. An annihilator is therefore an $L^q$ adjoint section, with
$1/p+1/q=1$, together with covectors at the fixed marks. The adjoint
equation has point sources at the marks and is homogeneous elsewhere.
Perturbations supported on projected injective disks away from the
marks force the adjoint section to vanish on those disks. Elliptic
regularity and unique continuation make it zero on each punctured
domain. An $L^q$ section supported at finitely many points is zero.
Arbitrary fixed-point normal values of smooth $\xi$ then force the
remaining point covectors to vanish. The augmented linear map is
surjective. For distinct curves one chooses disjoint projected
disks and applies the same argument to their direct-sum normal
operator, killing the adjoint sections independently.

Because the cokernel of the augmented operator on $\xi$ is finite
dimensional, finitely many smooth parameter variations suffice for
this surjectivity. Keep this finite-dimensional parameter space and
replace $W^{1,p}\to L^p$ by $W^{k,p}\to W^{k-1,p}$. Given a
$W^{k-1,p}$ right-hand side and prescribed evaluation values, solve
using the $W^{1,p}$ result. The parameter term is smooth, and elliptic
regularity upgrades the resulting section to $W^{k,p}$, with the
same values at the marks. This is the regularity upgrade used in
[Wendl, *Lectures on Symplectic Field Theory*, Lemma 9.2.9](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf#page=369),
now with the fixed-point evaluation factors retained.

Take $k\ge3$ sufficiently large and parameter regularity $\ell$
sufficiently larger than $k$ for the required differentiability.
The nonlinear evaluation maps with moving marks are then defined
in these $W^{k,p}$ and $C^\ell$ charts with the differentiability
needed for Sard--Smale. Their derivative adds
$du(z_i)\dot z_i$ to the fixed-point variations. Since each map is
immersed, these tangent directions together with the arbitrary
normal values span the full target tangent spaces. Thus the
coincidence constraint is transverse in the moving-mark problem.
For the fixed core, the variable curve already supplies arbitrary
full evaluation variations. The multi-evaluation setup is the one
in [Wendl, *Lectures on Holomorphic Curves*, version 3.3, Section 4.6](https://www.mathematik.hu-berlin.de/~wendl/pub/jhol_bookv33.pdf#page=200).
No differentiability of moving-point evaluation on $W^{1,p}$,
or $L^q$ description of a higher-Sobolev dual, is used.

In a path the two-marked self-intersection locus, and likewise the
intersection locus for two distinct curves, has expected dimension

$$
                         1+4-6=-1.                     \tag{13}
$$

The index for intersection with the fixed core is the same: the two
marks are one on the moving curve and one on the core. Hence all
these loci are empty for generic paths. Immersion was already proved
in §1, so no missing jet condition for critical points of the maps
is being suppressed.

Apply these statements in countably many normal/Fredholm charts and
in an exhaustion of the punctured tube. On each chart the required
perturbations are supported a positive distance from the core and
from the endpoint collars. The conditions (11) therefore do not
restrict any of the value or jet variations used above. Parametric
Sard--Smale and the Baire theorem give a simultaneous residual choice.
Elliptic regularity and smooth local perturbations give the smooth
version. $\square$

There is a dense set of parameters at which all these simple curves
are regular. Indeed, for each of the finitely many degree/genus
possibilities, the path moduli is a second-countable one-manifold and
the projection has only isolated critical points. The union of their
critical values is countable and has measure zero; its complement is
dense. This statement does not assert that the entire path moduli is
compact or that there are only finitely many simple-curve folds.

For the underlying Banach/Floer-space and smooth Baire argument, see
[Wendl, Sections 9.1--9.2 and Remark 9.2.6](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf#page=359).
That argument is used only after verifying (a)--(d) for the restricted
parameter space; ambient genericity alone would not imply them.

The normal-bundle formulation also explains why prescribing the core's
first jet is compatible with this assertion. At a noncore curve, its
projected injective disk lies away from the locus on which the jet was
prescribed. Keeping the zero section fixed would obstruct a generic
fold argument *at a singular zero section*; that issue is absent here
because its operator is invertible for the entire path.

## 5. Endpoint paths and optional cover bifurcation jets

For distinct even theta characteristics on a genus-two curve, both
holomorphic line operators are invertible: a degree-one line has at
most one independent holomorphic section, so even parity means none.
One can therefore connect their positive spin operators through
base-invertible scalar real operators. First lower the positive spin
mass to zero. The spin norm identity preserves injectivity at every
positive mass, and the zero-mass endpoint is invertible. Next join
the holomorphic lines through
$\operatorname{Pic}^1(C)\setminus\Theta$, lifting the path to Dolbeault
operators on a fixed smooth complex line. This complement is
path-connected: the theta divisor is a complex curve in the connected
complex surface $\operatorname{Pic}^1(C)$, and a general path avoids
its real codimension-two locus. Finally restore the other spin mass.

The resulting compact path has a uniform bounded inverse on the base.
It can be perturbed slightly to satisfy the double-cover conditions
of §2 without losing base invertibility. The lower term $B$ does not
affect base invertibility, since $D_2$ is invertible there. Choose
$B=0$ in endpoint collars. The higher-order perturbations of Proposition
3 can then be made with exactly these prescribed first jets.

If a detailed double-cover bifurcation model is also desired, its
nonlinear coefficient is accessible inside the same triangular
family. At a simple cover kernel $(s,v)$, choose a disk where $s\ne0$.
A real cubic fiber map from $L_1$ to either output line can take an
arbitrary prescribed value at $s$; on the other sheet it takes the
negative value. Thus cubic perturbations, vanishing through second
order at the core, can vary the anti-invariant cubic obstruction in
any cokernel direction. Combining this with the isomorphism (7)
gives a local reduced equation

$$
 x\bigl(L(\delta f,t-t_0)+x^2c+\text{higher terms}\bigr)=0.
                                                               \tag{14}
$$

The deck action is $x\mapsto-x$. Solve the bracket using the invertible
map $L$; the coefficient of $x^2$ in $t-t_0$ can be made nonzero.
This is the required nondegenerate pitchfork coefficient. It is an
additional finite-jet condition, not a consequence of (7) alone.
Multiplying all perturbations by an interior time cutoff preserves
every endpoint jet. Proposition 3 does not require this extra
condition for local isolation of the *simple* curves.

## 6. Precise interface with the cluster argument

The verified primary-source conditions are
[Doan--Ionel--Walpuski, Definition 3.13 and Proposition 3.30](https://arxiv.org/html/2103.08221v2):
their simple-isotopy argument uses isolated embedded disjoint curves
and an unobstructed core; its Hausdorff neighborhoods may change along
the path. It does not require one fixed disk tube with no curves on
its boundary for every parameter. This removes the particular
objection based on scaling solutions of a linear equation through a
fixed tube boundary.

Proposition 3 supplies the local isolation conditions at degree at
most two for a nonlinear triangular path with the specified
base-invertible core. The genus bound supplies the only possible
genera of simple degree-two correction cores, namely three and four.
In a compact realization, one must still justify that the localized
cycle neighborhoods and all the correction cores in the chosen
application remain inside the triangular region. Curves outside that
region are not covered by Proposition 3. One may perturb the ambient
structure away from a smaller closed tube, but no global
$\mathcal J_{\rm isol}$ assertion in arbitrary degrees is proved here.

The compactness statement needed locally is the bounded-degree,
bounded-region version of Calculation A1, now with (3). It bounds reduced
supports before applying compactness. It does not assert an arithmetic
genus bound for all stable maps with contracted components. The
cluster series includes those maps; they must be retained in the
cluster comparison and in the multiplicity-one elementary factor.

Accordingly, the local transversality obstacle identified for split
operators has been removed in degree two. If the separate localized
cluster comparison applies along this path, its correction cores
have degree two and genus at most four. At a degree-two cutoff those
cores occur only with multiplicity one. This is the precise input
needed to infer an elementary-exponent cutoff after identifying their
complete local series. The present note does not replace those
localization and series-identification arguments with a count of
smooth kernels.
