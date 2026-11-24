# A2. An all-genus double-cover comparison for distinct even spin lines

## The statement and the existing inputs

Let $C$ be a smooth compact complex curve of genus two, and let $a,b$ be
distinct even theta characteristics. On $N=a\oplus b$ take the direct
sum of their positive-mass real spin operators, with the absolute
determinant orientation. Write $L_{m,k}$ for the connected local
Gromov--Witten degree of covering degree $m$ and source genus $k$.
Write $n_{m,h}$ for the integral exponent of the balanced elementary
factor of degree $m$ and support genus $h$, as in
Article III, Proposition 3.1.
These are different invariants.

**Theorem.** For the indicated germ,

$$
              n_{2,3}=-8,\qquad n_{2,4}=2,
              \qquad n_{2,h}=0\quad(h\ge5).              \tag{1}
$$

Equivalently, with $s(t)=2\sin(t/2)$,

$$
           \sum_{k\ge3}L_{2,k}t^{2k-2}=-\frac12s(t)^4.   \tag{2}
$$

Consequently, in every source genus $k\ge3$,

$$
         L_{2,k}=\frac{(-1)^k(2^{2k-2}-4)}{(2k-2)!}.     \tag{3}
$$

If $F_h(q)$ denotes the balanced elementary series in the complex
cycle algebra, the unique balanced lift of the germ satisfies

$$
       \widetilde P_{a,b}(q)
          =F_2(q)F_3(q^2)^{-8}F_4(q^2)^2\pmod{q^3}.     \tag{4}
$$

The entries $-8$ and $2$ were already proved by the étale and
first-excess calculations in Papers 7--8. The independent genus-five
calculation in [Calculation B2](../excess/genus-five-mixed-spin-double-covers.md)
also gave $n_{2,5}=0$. The new assertion here is the uniform vanishing
in (1), established by controlling the genera of correction supports.

The proof uses the existing geometric factorization and identification
of its exponents with the connected elementary expansion from
Article III, Proposition 3.1 and Corollary 3.2. It uses
the published cycle compactness and cluster results of
[Doan--Ionel--Walpuski](https://arxiv.org/html/2103.08221v2), abbreviated
DIW. It does not assume that Bai--Swaminathan's all-genus wall invariant
has already been identified with elementary exponents, and does not
construct a numerical evaluation on arbitrary analytic families.

## 1. Triangular equations preserve the support bound

Work on fixed smooth degree-one lines $L_1,L_2$ over $C$. Fix a
positive spin operator $D_2$ on $L_2$ and a reference complex-linear
operator $\bar\partial_{L_1}$. Consider the equations

$$
 \bar\partial_{f^*L_1}u_1+f^*A_t(u_1)=0,\qquad
 f^*D_2u_2+f^*F_t(u_1)=0.                              \tag{5}
$$

Here $A_t$ and $F_t$ are smooth fiber maps from $L_1$ to the
$(0,1)$-form bundles with coefficients in $L_1$ and $L_2$,
respectively. The pullbacks include the differential of $f$ on the
form factor. We require $A_t(0)=F_t(0)=0$. These equations define a
graph almost-complex structure whose projections to $C$ and to the
almost-complex surface $L_1$ are holomorphic. Its fiber complex
structure is standard. A compact nonconstant curve cannot lie in
one affine fiber.

The spin/Stokes argument in Article III, Lemma 14.1, shows that $h^*D_2$ is injective
for every smooth connected finite cover $h$ of $C$. It is also
injective on connected stable sources: sections vanish on the
nonconstant normalized components, and matching propagates their
vanishing over the contracted components.

For a reduced graph of total degree $m$, use the normalized
off-diagonal fiber product of its normalization with itself over
$C$, as in [Calculation A1](genus-bound-for-split-multisections.md).
On every component, the first-coordinate difference $v_1$ solves a
homogeneous scalar real Cauchy--Riemann equation. If it were identically
zero, the second-coordinate difference would solve $h^*D_2v_2=0$,
and hence vanish. That would identify two generic sheets of the full
graph, contrary to reducedness and generic injectivity. Thus $v_1$
is nonzero on every retained component.

The scalar zero count, followed by the ramification and node count
of Calculation A1, gives

$$
       R(f)+2b\le m(m-1),\qquad
       p_a\le1+\frac{m(m+1)}2.                          \tag{6}
$$

It also proves that every simple full graph projects to a simple
scalar graph in $L_1$, and that distinct full graphs have distinct
scalar images. No coordinate is discarded in this argument.

In degree two the scalar image has a particularly rigid description.
For the smooth normalization, let $\iota$ be the deck involution.
The section $u_1-u_1\iota$ has precisely two zeros. If the cover is
branched, its two ramification points exhaust those zeros simply.
The scalar graph is then embedded and has genus four. At a branch
point with $f(w)=w^2$, the equation gives
$\partial_{\bar w}u_1(0)=0$, and the simple odd zero gives
$\partial_wu_1(0)\ne0$. If the cover is étale, the scalar graph is
immersed, has genus three, and has exactly one transverse double
point. Degree-one scalar graphs are embedded sections. In particular
every simple full graph of degree at most two is immersed before
any perturbation is made.

## 2. A path with regular core

Both even theta characteristics are non-effective. Indeed a line of
degree one on a genus-two curve has at most one independent section,
and even parity therefore forces $h^0=0$.

Keep $D_2=D_b$ fixed and join $D_b$ to $D_a$ on the first smooth line.
First lower the positive spin mass on $b$ to zero. At positive mass
the spin norm identity gives injectivity, and at zero mass it follows
from $h^0(b)=0$. Join the Dolbeault structures of $b$ and $a$ inside
$\operatorname{Pic}^1(C)\setminus\Theta$, and then restore the positive
spin mass on $a$. The complement of the theta divisor is path
connected: a smooth path in the complex surface $\operatorname{Pic}^1(C)$
can be perturbed to miss its real-codimension-two divisor. Smooth
Dolbeault representatives can be chosen on the fixed underlying line
over the interval. All these operators have index zero and zero
kernel, so they are invertible.

Denote this smooth base-invertible path by $D_{1,t}$. Initially set
$A_t$ equal to its zeroth-order term and $F_t=0$. Allow arbitrarily
small nonlinear perturbations subject to

$$
 \begin{aligned}
 A_t(0)&=F_t(0)=0,\\
 d_{\rm fib}A_t|_0&=D_{1,t}-\bar\partial_{L_1},\qquad
 d_{\rm fib}F_t|_0=0.                                  \tag{7}
 \end{aligned}
$$

Take the path constant in endpoint collars and keep the complete
positive-spin equations fixed there. The spin/Stokes argument ensures
that only the zero-section support occurs in those collars. The
normal operator of the zero section remains $D_{1,t}\oplus D_2$
and is invertible throughout. Its determinant sign is positive.
Allowing a lower off-diagonal first jet would also preserve this
invertibility, but is unnecessary here. The nonlinear freedom away
from the zero section will suffice.

The compact-realization construction in the proof of Article III, Proposition 3.1, realizes the entire
compact parameter family on one tube $V$ in a closed symplectic
six-torus, with common taming form and primitive class $A=[C]$.
The family agrees along the zero section, and compactness of the
interval permits one sufficiently small tube for the construction.
All local classes have the positive integral degree label $mA$.
Take an energy cutoff $2\int_C\omega\le\Lambda<3\int_C\omega$.

## 3. Actual-curve genericity in the triangular class

We require genericity for actual simple curves of degree one or two
contained in $V$. We do not require a split path to be transverse
to every smooth-cover wall, or to satisfy a fixed-genus wall theorem
on every compactified stable-map space.

**Lemma.** There are arbitrarily small smooth triangular perturbations
satisfying (7), fixed in the endpoint collars, for which all simple
curves of degree at most two in $V$ are embedded and pairwise disjoint
at every parameter. Their fixed-parameter moduli are discrete. A dense
set of parameters has every such curve unobstructed.

**Proof.** The zero section is regular by (7). For a noncore simple
graph $u:Y\to N$, write $\bar u:Y\to L_1$ for its scalar projection.
Both are immersions by Section 1, so their normal bundles fit into
an exact sequence of real Cauchy--Riemann operators

$$
 0\longrightarrow(f^*L_2,f^*D_2)
 \longrightarrow(N_u,D_u^N)
 \longrightarrow(N_{\bar u},Q_{\bar u})
 \longrightarrow0.                                    \tag{8}
$$

The vertical operator is $f^*D_2$ because the second equation in (5)
is linear in $u_2$; the quotient is the scalar normal operator of
the immersed graph in $L_1$. Thus
$\ker D_u^N\hookrightarrow\ker Q_{\bar u}$.
The full normal operator has real index zero, since the ambient
first Chern number on every local class is zero.

There is an injective unbranched disk of the scalar graph away from
the zero section and from the scalar images of any finitely many
other distinct curves. Distinct simple scalar curves have isolated
intersections; identical scalar images would contradict injectivity
of $D_2$ applied to their second-coordinate difference. On this disk,
arbitrary values of $(\dot A,\dot F)$ supply arbitrary normal
$(0,1)$-forms. They extend to parameter perturbations supported away
from the zero section. The adjoint unique-continuation argument
therefore makes the universal normal Fredholm section transverse.
At a fixed parameter its projection has Fredholm index zero; for
paths the index is one.

Next hold the curve fixed by making the perturbation values zero
along its scalar image. First normal derivatives remain arbitrary
real maps from the scalar normal line to the rank-two normal output.
At a one-dimensional kernel this makes the derivative to
$\operatorname{Hom}(\ker D_u^N,\operatorname{coker}D_u^N)$ surjective.
At nullity $k\ge2$, choose two independent kernel elements and project
them to scalar normal sections $q_1,q_2$. They are real-linearly
independent on some disk: otherwise their real ratio would satisfy
$\bar\partial a=0$, forcing a constant ratio and global dependence
by unique continuation. An annihilator for the variation of these
two kernel columns has pointwise expression

$$
                       q_1\otimes\eta_1+q_2\otimes\eta_2=0.
$$

On the disk of independence this forces both adjoint sections to
vanish, and unique continuation makes them zero globally. Hence
the derivative on the two columns is onto a space of dimension $2k$.
Finite-dimensional Fredholm reduction puts the nullity-$k$ locus in
a local submanifold of codimension at least $2k\ge4$. This is enough
for generic paths to avoid it; no maximal-codimension claim for
all kernel strata is needed.

At nullity one, second normal derivatives of $(A,F)$ can be varied
while their values and first derivatives remain fixed along the
scalar curve. A nonzero kernel has a nonzero scalar normal component
$q$; an adjoint generator $\eta$ is nonzero on an overlapping disk.
The change of the quadratic Kuranishi coefficient is a nonzero
constant times $\int\langle\eta,\dot P^{(2)}(q,q)\rangle$.
A real quadratic fiber map on a fixed nonzero vector has arbitrary
output. Thus vanishing of this coefficient imposes one further
independent equation. Parametric transversality gives a smooth
one-dimensional moduli over the path, with only Morse critical
points of its projection: locally $t-t_0=\pm x^2$. Every slice is
therefore locally discrete.

The two-point evaluation conditions are transverse with these same
restricted parameters. First check surjectivity of the augmented
linearization at fixed marked points on $W^{1,p}$, $p>2$, with target
$L^p$ plus the finite-dimensional evaluation spaces. An annihilator
consists of an $L^q$ adjoint section and covectors at the marked
points. Parameter values on projected injective disks away from the
marks kill the adjoint section there; adjoint unique continuation
kills it on the punctured domain. An $L^q$ section supported at
finitely many points is zero. Varying normal values and marked-point
positions then kills the point covectors. For distinct curves,
disjoint projected disks kill the adjoint sections separately.
The same argument works against the fixed core: the other curve
has projected injective disks off that core. For applying parametric
transversality with moving marks, use sufficiently high $W^{k,p}$
and $C^\ell$ charts, with $k\ge3$ and $\ell$ large, so that variable-point
evaluation has the required differentiability. The low-regularity
surjectivity just proved implies surjectivity there: a finite collection
of smooth parameter variations spans the finite-dimensional augmented
cokernel, and elliptic regularity upgrades solutions with
$W^{k-1,p}$ right-hand side to $W^{k,p}$. No differentiability of moving
evaluation on $W^{1,p}$ is assumed. The expected dimension of a
collision locus in a path is $1+4-6=-1$. Thus there are no
self-intersections, intersections between distinct curves, or
intersections with the core. Immersion was already known, so there
is no unstated critical-point transversality requirement.

These tests apply in countably many Fredholm charts. Each uses
perturbations supported off the core and away from endpoint collars;
the prescribed jets in (7) do not restrict them. Relative
Sard--Smale and a countable smooth local perturbation argument give
the asserted path. The critical points of a Morse function on a
second-countable one-manifold form a countable set; their values
therefore have dense complement. At those parameters every noncore
simple curve is regular, and the core is always regular. This proves
the final assertion. $\square$

The scalar zero and unique-continuation inputs are the usual
similarity principle; see
[Wendl, *Lectures on Symplectic Field Theory*, Theorem 2.5.3](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf).
The relevant Fredholm, evaluation, and parametric-transversality
framework is described in
[Wendl, *Lectures on Holomorphic Curves*, Chapter 4](https://www2.mathematik.hu-berlin.de/~wendl/pub/jhol_bookv33_2.pdf).
For the low-to-high Sobolev regularity step, see also Lemma 9.2.9
of the cited SFT lectures.
The first-column and second-jet tests above verify the restricted
parameter freedom instead of importing genericity from the full
ambient almost-complex space. Further detail is recorded in
[the transversality investigation](triangular-transversality.md).

## 4. Localized cluster transport at degree two

The preceding lemma supplies the local version of DIW Definition
3.13 in $V$ at energy cutoff $\Lambda$. We spell out the localization
of their Proposition 3.30, so that curves outside $V$ and arbitrary
coupled perturbations do not silently enter the genus estimate.

At each parameter choose Hausdorff neighborhoods of the core whose
closures consist of compact sets in a smaller tube $V'\Subset V$.
Every connected cycle in them of degree at most two is a
multiplicity-one embedded curve or twice a degree-one curve.
Two distinct degree-one components cannot form a connected support,
by disjointness. The embedded support spaces of fixed degree are
discrete, by the lemma and DIW Propositions 2.21 and 2.24. They are
countable. Apply bounded-energy cycle compactness, DIW Theorem 2.11,
to all connected cycles supported in a small closed Hausdorff ball,
including doubled degree-one curves. Their distances from $C$ form
a countable compact set of radii, including radius zero. Its complement
contains arbitrarily small radii with no cycle support on the boundary,
as in their Proposition 3.17. All limiting cycles under consideration
remain in $V$. In particular a sequence of degree-two curves tending
to $2C$ has not been omitted from the compact set. No genericity
outside $V$ is used.

Invertibility of the core operator makes $C$ the unique degree-one
support in such a sufficiently small neighborhood. DIW Proposition
3.10 transports its compact open-and-closed cycle insertion over a
short interval. The parametrized implicit function theorem retains
uniqueness in degree one after shortening the interval. Finitely
many intervals cover the path.

On every overlap choose a regular slice from the dense set in the
lemma. Refine the two cluster neighborhoods to a common smaller one,
as in DIW Proposition 3.19. The difference of the selected cycle
sets is compact and open and closed; it contains no degree-one
support or double of such a support. Its connected cycles are
therefore simple degree-two supports. They form a compact discrete
set, hence a finite set. Their cores are regular on the chosen slice.
Every one is an actual triangular curve in $V$, so its genus is
three or four by (6).

Adding the refinement differences along the finite interval cover
expresses the difference between the two complete connected local
series as a finite signed sum of these degree-two cluster series.
Each correction core is primitive at this cutoff: its proper
multiples have degree at least four. Its contribution is therefore
its sign times $G_h(q^2,t)$ modulo $q^3$. To apply the degree-one base
case of DIW Proposition 3.54 with its global genericity hypothesis,
first transport a sufficiently small isolating insertion around this
regular core to a nearby globally generic ambient structure. DIW
Propositions 3.20--3.21 preserve the insertion; the implicit function
theorem and openness of the embedded multiplicity-one locus preserve
its unique core, genus, and sign. One may also choose this perturbation
fixing the core, as in the proof of Proposition 3.54. The primitive
base case then applies. This final deformation may leave the
triangular class, since all its further proper-multiple corrections
are above the cutoff and the genus of the already selected core is
fixed.

This argument uses complete support insertions, not only maps with
smooth domains. All contracted components and stable-map boundaries
over each support remain included. In particular the degree-one
part of $G_h$ is $s(t)^{2h-2}$ and has arbitrarily high source-genus
terms. Bounding the genus of the correction core does not remove
these terms.

Consequently the endpoint difference has the form

$$
          \mathsf L_{a,b}-\mathsf L_{b,b}
             =c_3G_3(q^2,t)+c_4G_4(q^2,t)\pmod{q^3},   \tag{9}
$$

with integers $c_3,c_4$. This is the localized numerical cluster
identity. The endpoint factorization and its identification with the
same numerical elementary exponents, given by Article III, Proposition 3.1 and Corollary 3.2, imply that (9) compares those universal exponents. This
step requires no new universal analytic-family GW evaluation. The
full source-by-source localization is also recorded in
[the cluster investigation](triangular-cluster-transport.md).

## 5. Evaluation and scope

The equal-spin endpoint is elementary and has no degree-two
elementary exponents. Equation (9) therefore proves that the mixed
endpoint has no such exponents except in genera three and four.
The earlier étale and first-excess calculations give $c_3=-8$ and
$c_4=2$, proving (1) and (4).

The two degree-two partitions have hook lengths $1,2$, so

$$
 \begin{aligned}
 [q^2]G_2(q,t)
   &=2s(t)^2s(2t)^2-\tfrac12s(t)^4\\
   &=\tfrac{15}{2}s(t)^4-2s(t)^6.
 \end{aligned}
$$

Adding $-8s(t)^4+2s(t)^6$ gives (2). Expanding the elementary
trigonometric identity gives (3), with the independently known values
$L_{2,3}=-1/2$, $L_{2,4}=1/12$ and $L_{2,5}=-1/160$. The next values
are $17/60480$ and $-31/3628800$ in source genera six and seven.
These last values are consequences of the all-genus proof, not new
direct integrations on the higher-genus cover spaces.

The exact check
`code/all_genus_double_covers.py`
independently expands the sine products and the cosine expression,
compares the previously computed source genera three through five,
and performs triangular exponent extraction through genus twenty-four.
It checks the formal arithmetic, not the analytic or cluster proof.

The result establishes the degree-two part of the proposed
self-conjugate-partition formula for distinct even spins in base
genus two. It does not prove that formula in degrees three and
higher. The automatic immersion and scalar-image classification
used in Section 3 are specific to the present low degree and base
genus. Extending the transversality argument beyond this range
requires controlling critical points of the scalar projection.

Nor does (9) assign a correction to each singular parameter of the
original split linear path. It obtains a finite cluster comparison
through a nonlinear triangular path. The precise individual-wall
law, general source/base genera, arbitrary multiplicity, and the
assembled ten-objective programme remain open.
