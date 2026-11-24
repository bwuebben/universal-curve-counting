# A6. Triangular graph equations and degree-two cluster transport

**Conclusion.** Once a triangular path has the local isolation properties
below, no new all-genus wall formula is needed to obtain the degree-two
cutoff. A localized form of DIW Proposition 3.30 gives finitely many
correction cores of degree two. Each such core is primitive at this
cutoff, so its entire contribution, including contracted components, is
its signed elementary series. In genus two the support estimate bounds
their genera by four. For distinct even spin endpoints the degree-one
core can remain unobstructed throughout the path, avoiding the
two-parameter argument of DIW Proposition 3.37.

The substantive outstanding hypothesis is a genericity theorem **for
actual simple triangular curves**, including embeddedness and discreteness
at every parameter. Surjectivity for a smooth-cover linear wall alone
does not establish this hypothesis. A single tube isolating the entire
path is unnecessary for the argument below.

## 1. The geometric setting

Let $C$ have genus $g\geq2$, let $N=L_1\oplus L_2$ with
$\deg L_j=g-1$, and consider graph equations

$$
 \bar\partial_{f^*L_1}u_1+f^*A_t(u_1)=0,
 \qquad f^*D_2u_2+f^*F_t(u_1)=0.                 \tag{1}
$$

Here $A_t$ and $F_t$ are smooth fiberwise maps with values in the
appropriate $(0,1)$-form bundles; $A_t(0)=F_t(0)=0$; and $D_2$ is a
fixed positive-mass spin operator. Pullback includes the differential
of $f$ on the form factor. Thus the zero section remains holomorphic.
All compact curves have nonconstant projection: a curve lying in an
affine fiber is constant because the vertical complex structure is the
standard one.

There is a useful strengthening of the support argument in
[Calculation A1](genus-bound-for-split-multisections.md). On a component
of the normalized off-diagonal fiber product, write $v_j$ for the
difference of the two $j$-th coordinates. If $v_1$ vanishes identically,
then $v_2$ solves the pullback homogeneous equation for $D_2$. Spin
super-rigidity forces $v_2=0$. For a reduced simple graph this is
impossible. Therefore $v_1$ is nonzero on **every** off-diagonal
component. Its scalar difference equation supplies the complete zero
count without changing coordinates from component to component:

$$
 R(f)+2b\leq m(m-1)(g-1),\qquad
 p_a\leq1+\frac{m(m+1)}2(g-1).                  \tag{2}
$$

In particular, the projection of a simple full graph to
$\operatorname{Tot}(L_1)$ is also simple. In degree two every embedded
full graph has genus at most $3g-2$; the bound is four when $g=2$.
The first-coordinate image need not be embedded. The second coordinate
can separate its isolated double points.

Use the compact realization in the proof of
Article III, Proposition 3.1:
a compact smooth family fixing the zero section is realized on one
common tube $V$ in a closed symplectic six-manifold, with a common tame
form and a primitive core class $A=[C]$. All the graph structures have
the same value along $C$, so compactness of the parameter interval gives
one tube for the cutoff construction. The six-torus realization makes
$c_1=0$ and retains the positive integral degree label.

Choose $a=\int_C\omega$ and $2a\leq\Lambda<3a$. Every curve contained
in $V$ and of energy at most $\Lambda$ has degree one or two. In this
note “degree two” always refers to the projection/homology label, not to
the multiplicity of the underlying simple support.

## 2. The precise local hypotheses

DIW's space $\mathcal J_{\rm isol}$ is defined by four conditions:
nonnegative index for simple curves, embeddedness of index-zero simple
curves, disjointness of distinct such images, and discreteness of their
moduli at a fixed structure. Their $\mathcal J_*$ replaces discreteness
by unobstructedness. Proposition 3.30 requires the core path to be
unobstructed and the path to lie in $\mathcal J_{\rm isol}$; it does not
require every simple curve along the path to be unobstructed. These
are [DIW, Definition 3.13 and Proposition 3.30](https://arxiv.org/html/2103.08221v2#S3).

The following localized hypotheses suffice here. They concern every
curve contained in $V$, and not only graphs arbitrarily close to zero.

1. The normal operator of $C$ is invertible for every $t$.
2. At every $t$, all simple curves in $V$ of degree one or two are
   embedded, and the moduli of each such degree is discrete in the
   Gromov topology.
3. Distinct degree-one supports in $V$ are disjoint. It is sufficient,
   and convenient for genericity, to require all distinct supports of
   degree at most two to be disjoint. Only pairs whose degrees sum to
   at most two enter the present connected-cycle argument.
4. A dense set of parameter values has all these simple curves
   unobstructed. Equivalently for this argument, every nonempty overlap
   of the finitely many parameter intervals chosen below contains a
   slice on which its correction cores are unobstructed.

The index requirement is automatic locally: every class is $mA$ and
$c_1(A)=0$, so the unmarked simple-map index in dimension six is zero.
The hypotheses allow ordinary folds of degree-two curves. They do
not allow a positive-dimensional family of simple degree-two curves
on a fixed slice. An exact linear kernel equation therefore fails
condition 2 because nonzero solutions can be scaled.

Embeddedness alone is insufficient; regularity on a dense set alone is
also insufficient. The all-parameter discreteness condition is what
permits the cluster neighborhoods at the singular parameters. A
statement that a path is “generic in degree two” must explicitly imply
these conditions or prove a substitute for them.

## 3. Local degree-two simple isotopy

**Proposition.** Suppose a path as in (1) has a common compact tame
realization and satisfies conditions 1--4. Assume its endpoint germs
are super-rigid and contain only the zero-section support in a smaller
tube. Let $\mathsf L_0(q,t)$ and $\mathsf L_1(q,t)$ be their connected
local Gromov--Witten series. There are finitely many embedded triangular
curves $C_i$, of degree two and genera $h_i\leq3g-2$, and integers
$c_i$, such that

$$
 \mathsf L_1(q,t)-\mathsf L_0(q,t)
       =\sum_i c_iG_{h_i}(q^2,t)\pmod{q^3}.       \tag{3}
$$

Here $G_h$ is the connected elementary series of DIW Definition 3.45.
There is no sign change multiplying the first endpoint series: the
core is unobstructed throughout.

**Proof.** The proof of DIW Proposition 3.30 localizes to these
hypotheses as follows.

First choose every Hausdorff neighborhood of $C$ so small that its
closure consists of compact subsets of a smaller tube
$V'\Subset V$. Bounded-energy cycle compactness takes place in the
closed ambient manifold. Any limiting cycle belonging to one of these
closures is still contained in $V$, so its local degree, index,
embeddedness and discreteness are controlled. Nothing is assumed
about curves elsewhere in the closed manifold.

At a fixed $t$, each connected cycle in these neighborhoods is either
a multiplicity-one embedded curve of degree one or two, or twice an
embedded degree-one curve. A connected sum of two distinct degree-one
supports is excluded by disjointness. The support spaces of embedded
curves of fixed degree are discrete: use discreteness in the Gromov
topology together with DIW Propositions 2.21 and 2.24, which identify
the relevant embedded-curve, cycle and support topologies.

The argument of DIW Proposition 3.17 now gives arbitrarily small
Hausdorff radii with no cycle on the boundary. The exceptional radii
are a countable compact set, obtained from bounded-energy cycle
compactness and the discrete support spaces. Invertibility of the
degree-one core makes $C$ the unique nearby degree-one curve. Thus
the chosen neighborhoods are clusters at the fixed parameter.

DIW Proposition 3.10 makes each such neighborhood an open-and-closed
cycle insertion over a short parameter interval. This step uses
closedness of the boundary-cycle locus and properness of the cycle
projection. It does not require regularity or genericity of curves
outside $V$. The parametrized implicit function theorem keeps $C$
as the unique degree-one support in the neighborhood along that short
interval. Choose finitely many such intervals covering $[0,1]$.

At an overlap, choose a slice allowed by condition 4. Refine the two
cluster neighborhoods to a common smaller one, as in DIW Proposition
3.19. The difference between a larger and smaller neighborhood is a
compact open-and-closed set of cycles. It contains no degree-one
support. Every remaining connected cycle is therefore a simple
degree-two support. The remaining support space is both compact and
discrete, hence finite. Its cores are unobstructed by the choice of
slice. This is the entire finite cluster decomposition required at
degree two; no higher-degree induction is involved.

The interval contributions are constant, and the refinement
differences add with signs. Summing along the finite cover gives the
endpoint difference as a finite signed sum of degree-two cluster
contributions. All their cores are actual curves at the selected
parameters, inside $V$, so (2) bounds their genera.

Finally, a degree-two core is primitive relative to this energy
cutoff: its proper multiples have degree at least four. Its cluster
contribution is its signed $G_{h_i}(q^2,t)$ modulo $q^3$. One may use
the degree-one case of the elementary comparison, or the base case
of DIW Proposition 3.54. The latter deformation may leave the
triangular class: every additional correction lies above the cutoff,
and the genus $h_i$ of the existing core is fixed. This proves (3).
$\square$

This is a localization of the published proof, not a claim that DIW
state the localized proposition verbatim. Its precise source chain is
Theorem 2.11; Propositions 2.21 and 2.24; Lemma 3.8; Propositions
3.10, 3.17, 3.19 and 3.30; and the degree-one base case of Proposition
3.54. The proof above replaces the global cluster decomposition by
the compact-discrete degree-two argument. It shows why no uncontrolled
perturbation is used when selecting the genera $h_i$.

## 4. Ghosts and the universal exponents

The proposition is an identity of **complete connected generating
series**, not a comparison of smooth-domain Euler numbers. The cycle
insertion pulls back along the stable-map-to-cycle map. Therefore it
includes every contracted component and every stable-map boundary
stratum over the selected support. Contracted components do not create
new supports and cannot move a stable map between two open-and-closed
support insertions.

For example, the coefficient of $q$ in the elementary connected
series is

$$
                  [q]G_h(q,t)=(2\sin(t/2))^{2h-2}.         \tag{4}
$$

It generally has infinitely many source-genus coefficients. Thus
bounding the genera $h_i$ in (3) does not discard higher-genus maps
to the degree-two correction cores. Their complete contributions
are already present in (4). No assertion of a genus bound on all
stable maps is made.

For the endpoint germs, use the existing identification in
Article III, Proposition 3.1, between the universal
elementary exponents and the coefficients of their connected
numerical elementary expansion. Equation (3), uniqueness of that
expansion, and constancy of the degree-one core sign imply

$$
                 n_{2,h}(D_1)=n_{2,h}(D_0)
                    \quad (h>3g-2).                       \tag{5}
$$

Equivalently, their balanced lifts satisfy

$$
 \widetilde P_1\widetilde P_0^{-1}
       =\prod_{h\leq3g-2}F_h(q^2)^{c_h}\pmod{q^3}.        \tag{6}
$$

The passage to (6) uses the already established endpoint factorization
and its identification of the same integer exponents. DIW alone is a
numerical theorem and does not supply a universal evaluation on
arbitrary parameter families. No new such evaluation is being assumed
here.

If the initial endpoint is the equal-spin elementary germ, its
$n_{2,h}$ vanish identically. For a path satisfying the proposition,
(5) then gives

$$
                      n_{2,h}=0\qquad(h\geq3g-1).         \tag{7}
$$

In genus two this is the desired $n_{2,h}=0$ for $h\geq5$.
It remains conditional on constructing the path with conditions 1--4.
The conclusion is stronger than the smooth-cover linear wall
comparison and does not identify individual wall contributions.

## 5. An unobstructed core path for the genus-two even-spin endpoints

Let $a,b$ be distinct even theta characteristics on a genus-two curve.
They have $H^0(a)=H^0(b)=0$: a line of degree one has at most one
independent section, whereas even parity forces the dimension to be
even. Fix $D_b=\bar\partial_b+R_b$ in the second coordinate. On the
first underlying smooth line, join $D_a$ to an identified copy of
$D_b$ as follows.

First lower the mass in $\bar\partial_a+sR_a$ from one to zero.
The spin norm identity gives injectivity for $s>0$, and $H^0(a)=0$
gives it at $s=0$. The real index is zero, so these operators are
invertible. Next join the holomorphic line $a$ to $b$ inside
$\operatorname{Pic}^{g-1}(C)\setminus\Theta$, where $\Theta$ is the
effective theta divisor. This complement is path connected: the
divisor has real codimension two in the connected complex torus and
can be avoided by a perturbation of a path. Choose smooth
representatives of the resulting Dolbeault operators on the fixed
smooth line bundle. Every operator on this segment has zero kernel
and index zero. Finally raise the spin mass on $b$ from zero to one.

The resulting base path is invertible. Arbitrary lower off-diagonal
terms preserve base invertibility because the two diagonal base
operators are invertible. Small constrained perturbations of the
first diagonal term also preserve it. Therefore the ordinary core
wall of DIW Proposition 3.37 can be avoided in this even-spin case.
There is no need to establish a two-parameter triangular unfolding
of that wall before using Proposition 3.30.

This argument uses nonspecial spin endpoints. In higher genus an
even theta characteristic can have nonzero even-dimensional section
space, so the mass-to-zero segment needs separate attention. For
opposite base signs an everywhere invertible base path is impossible.

## 6. Why a global isolating tube and a projective compactification are different

Article III, Proposition 16.5,
requires one compact open-and-closed insertion for the entire
selected family. That is a sufficient hypothesis for equality of
the endpoint products with no correction factors. It is stronger
than what is needed for (3). Here the neighborhoods change along
the path, and their differences produce exactly the degree-two
correction factors. Linear scaling at a wall obstructs a fixed
isolating disk; it does not obstruct this finite refinement argument
once nonlinear perturbations give fixed-slice discreteness.

Compactifying the fibers, by itself, would not isolate an affine
local contribution. In the ruled surface compactification of a
positive line $L$, let $C_0$ and $C_\infty$ be its zero and infinity
sections and let $f$ be a fiber class. Then

$$
 C_0=C_\infty+(\deg L)f,\qquad C_\infty^2=-\deg L.
$$

An affine class $m[C_0]$ can have a boundary configuration containing
$mC_\infty$ and $m\deg(L)$ fiber components. Its zero intersection
with the infinity divisor does not exclude this configuration:
the negative self-intersection of the infinity component cancels
the positive fiber intersections. Thus positivity of intersections
alone cannot prove confinement. A compactification is useful as an
ambient space for compactness, but its boundary configurations must
not be silently identified with the original local invariant.

The small common tube in the compact realization, together with
changing Hausdorff cluster neighborhoods, avoids this extra global
problem. Nonlinear confinement may still be useful for a different
construction, but it is not a missing hypothesis of the localized
Proposition 3.30 argument.

## 7. Relation to the Bai--Swaminathan comparison

[Bai--Swaminathan](https://arxiv.org/html/2106.01206v4) define their
linear double-wall invariant separately from elementary-cluster
coefficients; their general identification is Conjecture 1.10.
Definition 6.5 and Theorem 6.7 impose a finite smooth-domain
degeneracy locus. Attaching contracted components to a low-genus
kernel violates that hypothesis in higher source genera. The
proposition here does not apply their Euler wall theorem to such
nodal loci. It instead uses compact support insertions and the
existing elementary comparison at a primitive correction core.

There is therefore no contradiction with the unresolved general
comparison in that source. The proposed result concerns a specially
controlled family with an a priori support bound and an everywhere
unobstructed degree-one core.
