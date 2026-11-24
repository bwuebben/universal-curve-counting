# A5. Split double-cover bounds and the all-genus wall comparison

Let $C$ have genus $g\ge2$, and let $N=L_1\oplus L_2$ with
$\deg L_j=g-1$. Write $D=D_1\oplus D_2$ for a real normal
Cauchy--Riemann operator. At a super-rigid endpoint, distinguish:

1. $L_{2,k}(D)$, the absolute Euler number on the **compactified**
   connected stable-map space of degree two and source genus $k$.
2. $n_{2,h}(D)$, the elementary-factor exponent, recovered from all these
   Euler numbers by the triangular rule in
   Article III, Proposition 3.1.
3. $w_{2,h}(D)$, Bai--Swaminathan's linear double-wall invariant.

These have different definitions. In particular, “absence of smooth
genus-$h$ kernels” is not a statement about the compactified Euler number
in source genus $h$.

The primary-source versions checked were
[Bai--Swaminathan, arXiv:2106.01206v4](https://arxiv.org/html/2106.01206v4)
and [Doan--Ionel--Walpuski, arXiv:2103.08221v2](https://arxiv.org/html/2103.08221v2),
the latest versions displayed by arXiv on the investigation date.
The former defines $w$ using generic paths (Definitions 5.14 and 5.17),
and constructs a global degree-two invariant (Theorem 1.6).
Its identification with elementary-cluster coefficients is explicitly
**Conjecture 1.10**, proved there only in source genus zero (Theorem 6.12).
Theorem 6.7 treats a finite smooth-domain degeneracy locus under Definition
6.5. The discussion following Theorem 1.11 identifies removal of these
restrictions as the obstacle to the general comparison.

Doan--Ionel--Walpuski Theorem 3.25 gives cluster transport with corrections
in proper multiples. Propositions 3.48 and 3.54 give integral elementary
coefficients with finite genus support; they do not give the bound four
here. Their proof uses generic almost-complex paths, without preserving a
normal splitting or a holomorphic projection. These are the numerical
inputs used in the proof of Article III, Proposition 3.1.

## 2. A split branched wall cannot be a simple wall

Put

$$
                     h=2g-1+r,\qquad r>0.
$$

A smooth double cover $f:Y\to C$ of genus $h$ has $2r$ branch points.
Fix $C$ and consider a path of split operators $D_t=D_{1,t}\oplus D_{2,t}$.
Suppose $D_0$ is invertible on $C$ and

$$
                         \dim_{\mathbb R}\ker f^*D_0=1.
$$

**Proposition.** Such a crossing cannot satisfy the isomorphism condition
in Bai--Swaminathan Definition 6.5. The failure of surjectivity has real
dimension at least $2r$.

**Proof.** Since the kernel splits, after relabeling it lies entirely in
$f^*L_1$, and $f^*D_{2,0}$ is injective. Riemann--Roch gives, for both lines,

$$
 \operatorname{ind}_{\mathbb R}(f^*D_{j,0})
 =2\bigl(2(g-1)+1-h\bigr)=-2r.
$$

Consequently

$$
 \dim\operatorname{coker}f^*D_{1,0}=1+2r,
 \qquad
 \dim\operatorname{coker}f^*D_{2,0}=2r.
 \tag{1}
$$

The base operator is invertible, so the reduced kernel and cokernel agree
with the unreduced spaces in this calculation. A deformation of $f$ and
of the split path gives a family of operators still preserving the two
pulled-back lines. Extend a kernel generator in the first line. Applying
the family of operators, differentiating, and projecting to cokernels
produces no component in the second cokernel. The same conclusion is
independent of choices of extension and trivialization, since their
changes lie in the image of the central operator.

The map required by the simple-wall hypothesis is

$$
 T_f\mathcal M_h(C,2)\oplus\mathbb R
 \longrightarrow
 \operatorname{Hom}_{\mathbb R}
       (\ker f^*D_0,\operatorname{coker}f^*D_0).
 \tag{2}
$$

Both sides have real dimension $4r+1$: the smooth double-cover locus
has $2r$ complex branch parameters. But the image of (2) lies in the
first cokernel block, of dimension $2r+1$. The other block has the
strictly positive dimension $2r$, proving the assertion. The argument
also applies if the complex structure of the base varies, because that
variation still preserves the direct sum. $\square$

For $g=2$ and $h=4$, the source and target of (2) each have dimension five,
while its image has dimension at most three. Thus the first branched
wall already fails the published simple-wall hypothesis. If both lines
have nonzero kernels, the total kernel has dimension at least two and
the simple-wall hypothesis fails earlier.

This does not prove that the split degeneracy locus has its naively
expected dimension: additional constraints can make it more singular.
It proves the failure of the required linear surjectivity without a
dimension heuristic. Off-diagonal operator variations are exactly the
directions that can supply the missing block. They leave the setting in
which the scalar zero count was proved.

## 3. Étale crossings and the reference operator

### A split étale wall is not an embedded double graph

Suppose now $g=2$, $f:Y\to C$ is connected étale of degree two,
$D=D_1\oplus D_2$ is base-invertible, and $\ker f^*D$ has real
dimension one. Again the kernel lies in one line. Its generator $u$ is
odd: the invariant part descends to a base kernel, which is zero.
If $\alpha$ is the flat sign line of $f$, the section $u$ descends to
a solution on $L_j\otimes\alpha$, a complex line of degree one.
The scalar similarity principle gives exactly one simple zero downstairs.

At the two points of $Y$ over that zero, the graph $y\mapsto(f(y),u(y),0)$
has the same image. It is therefore not injective. The simple zero gives
distinct tangent lines, so this is the one-double-point case of the
line-operator note. The failure is not fixed by demanding that the path
be transverse within the split locus.

The exclusion of nonembedded kernel graphs in Bai--Swaminathan's
Definition 5.14, via Definition 5.10, is thus a genuine restriction in
this problem. The all-genus kernel-free statement for smooth covers
of genus at least five does not put the entire split path into their
admissible path space.

### No balanced split holomorphic reference is double-cover rigid

**Proposition.** Let $g\ge2$ and let $L$ be a holomorphic line bundle of
degree $g-1$. Either $H^0(C,L)\ne0$, or some connected smooth double
cover $f:Y\to C$ satisfies $H^0(Y,f^*L)\ne0$.

**Proof.** Riemann--Roch gives

$$
 h^0(C,L^2)-h^0(C,K_C\otimes L^{-2})=g-1>0.
$$

Choose $0\ne s\in H^0(C,L^2)$ and form the finite spectral curve
$w^2=s$ in $L$. Normalize it. If the normalization is connected,
its projection is a connected smooth degree-two cover, and the
tautological section $w$ pulls back to a nonzero holomorphic section
of $f^*L$. This remains true when the original spectral curve has
multiple zeros or nodes. If the normalization is disconnected, its
two components have projection degree one. The tautological coordinate
on either component gives a global square root of $s$, hence a nonzero
section of $L$. $\square$

For a complex-linear operator on $L_1\oplus L_2$ preserving the splitting,
the proposition supplies a base or double-cover kernel. Hence the
complex-linear double-cover-rigid reference used to normalize $w$ cannot
be chosen in this balanced split locus. Replacing it by the equal-spin
**real** elementary operator is possible as a reference for the local
factorization, but requires a comparison of the resulting normalization;
it is not Bai--Swaminathan's defining reference.

## 4. Contracted components and exact linear transport

There are two different compactness issues.

First, let $f:Y\to C$ have a nonzero kernel section $u$ at a wall.
Attach a smooth contracted curve $T$ of genus $a\ge1$ at a point
$p\in Y$, mapping it constantly to $f(p)$. Extend $u$ over $T$ by
the constant value $u(p)$. On a contracted component the pullback
operator is the usual Dolbeault operator with constant vector-space
coefficients, so the extension is a kernel section of the nodal map.
The source genus is now $g(Y)+a$. Varying the attachment and the marked
curve $T$ gives a positive-dimensional degeneracy locus.

Thus a single low-genus wall is simultaneously a non-simple wall in
arbitrarily high source genera. This is an actual construction of nodal
kernel strata, not just the possibility that ghosts might occur. It
does not assert that these strata contribute nonzero Euler numbers or
are limits of smooth kernel maps. Their contribution is precisely what
an all-genus comparison must compute or cancel.

Second, at a kernel of an **exact linear** graph equation, every real
multiple $\lambda u$ is a solution. For an odd nonzero section these
give nonconstant degree-two graph curves tending to the zero section.
For any disk tube of radius $\rho>0$, scaling to
$|\lambda|\|u\|_\infty=\rho$ gives a curve contained in the closed tube
and meeting its boundary. Therefore one cannot transport the isolated
point insertion $2[C]$ through this wall using the boundary-free tubes
employed at super-rigid endpoints in Article III, Proposition 15.4.

A nonlinear perturbation or another compact local problem is needed
at the wall. Calculation A1 proves the zero estimate for smooth nonlinear
equations with independent line coordinates. It does not prove spatial
isolation or transversality for such a perturbation, nor does it cover
coupling between the coordinates. The required virtual boundary
contribution therefore remains to be determined.

The published exclusion of ghosts in limits of embedded curves along a
generic path (Bai--Swaminathan Theorem 1.8) does not remove either issue:
the split path above violates the genericity conditions, and the
compactified Euler integral includes ghost strata independently of
whether they are limits of embedded curves in that generic family.
