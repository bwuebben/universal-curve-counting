# C1. Triple covers with two holomorphic spin sections

Let $C$ be a smooth curve of genus two, with hyperelliptic involution
$\iota$, and let $a,b$ be distinct even theta characteristics. They
are non-effective. Fix their square isomorphisms with $K_C$ and let
$S_\theta(x,p)$ be the Szegő kernel of $\theta$, normalized to have
diagonal residue one. The difference of their squares is the symmetric
holomorphic bidifferential

$$
 \mathscr D(x,p)=S_a(x,p)^2-S_b(x,p)^2
                   \in H^0(C\times C,K_C\boxtimes K_C).
$$

Write $Q(p)=\mathscr D(p,p)\in H^0(C,K_C^2)$.

**Theorem.** Smooth connected degree-three covers $f:Y\to C$ with
$g(Y)=5$ and

$$
                 H^0(Y,f^*a)\ne0,\qquad H^0(Y,f^*b)\ne0             \tag{1}
$$

are classified, up to isomorphism over $C$, by the distinct zeros
of $Q$. The cover indexed by $p$ has exactly one ramification point,
of index three, over $\iota p$. Both spaces in (1) have dimension one.
The monodromy is $S_3$ and the cover has no nontrivial automorphism
over $C$.

The section $Q$ is the pullback of a nondegenerate binary quadratic
under the canonical double cover $C\to\mathbb P^1$. Consequently
there are four such covers generically and between two and four
distinct covers for every fixed pair. A zero at a Weierstrass point
has order two in $Q$, but still gives one smooth cover; it is not
discarded. The statement classifies geometric points and their
ordinary stack stabilizers, not intersection multiplicities of the
kernel loci.

For every resulting cover, attaching a pointed elliptic curve at its
triple ramification point gives a complete curve of simultaneous
holomorphic kernels inside the reduced main genus-six locus $M_3$ of
[Calculation B3](../excess/second-excess-reduction-in-all-degrees.md). Thus the
common kernel locus there is nonempty even after the
full étale-lift summands have been removed. This is a geometric
boundary obstruction to a disjoint-kernel argument. Calculation C2 evaluates the common-kernel intersection, and Article V, Corollary 1.2 computes the elementary cubic exponents.

## 1. The trace-free bundle determines an elementary modification

Put $f_*\mathcal O_Y=\mathcal O_C\oplus T$. Its trace-free bundle
has rank two and degree $4-g(Y)=-1$. Since $a,b$ have no base sections,
any sections in (1) give a map

$$
                     a^{-1}\oplus b^{-1}\longrightarrow T.         \tag{2}
$$

It has generic rank two. Otherwise the joint graph would lie in a
line subbundle $M\subset a\oplus b$ with both projections nonzero.
Then $\deg M\le0$: degree one would identify $a$ and $b$. The graph
has degree three over $C$, since the field extension has prime degree
and a degree-one graph would give base sections. Its monic spectral
equation has arithmetic genus $4+3\deg M\le4$. Its normalization
cannot have genus five, a contradiction.

Thus (2) has a torsion cokernel of length one, at a point $p$.
For $E=T^*$ it gives

$$
 0\longrightarrow E\longrightarrow a\oplus b
                       \longrightarrow\mathbb C_p\longrightarrow0.
                                                                    \tag{3}
$$

The quotient uses both summands. If, for example,
$E=a(-p)\oplus b$, the four summands of
$\operatorname{Sym}^3E\otimes(\det E)^{-1}$ are
$b(-2p),a(-p),b,a(p)$. Only the last has a section, so every nonzero
cubic datum is a pure cube with identically zero discriminant.
It cannot define the generic algebra of a smooth connected cover.
The other axis is the same with $a,b$ interchanged.

We use the cubic-algebra convention in which $E$ is the dual
trace-free bundle, so cubic data lie in
$H^0(\operatorname{Sym}^3E\otimes(\det E)^{-1})$. This is the dual
of the trace-free-bundle convention in
[Miranda, Theorem 3.6](https://www.math.colostate.edu/~miranda/preprints/TripleCoversInAG.pdf).

## 2. The Szegő extension condition

Choose a coordinate $z$ centered at $p$ and local frames $X,Y$ of
$a,b$ whose squares are $dz$. Normalize the meromorphic sections
$A=S_a(\cdot,p)$ and $B=S_b(\cdot,p)$ using the corresponding fiber
frames. Their scalar expansions are

$$
 A=z^{-1}+c_a z+e_a z^2+O(z^3),\qquad
 B=z^{-1}+c_b z+e_b z^2+O(z^3).                                   \tag{4}
$$

The constant term vanishes by antisymmetry of the even-spin Szegő
kernel in compatible coordinate frames. Each space
$H^0(\theta(p))$ has dimension one, and the residue determines its
meromorphic section. These are standard consequences of the prime-form
expression for the Szegő kernel; a primary account is
[Matone--Volpato, §1, p. 741](https://research.unipd.it/retrieve/e14fb267-790a-3de1-e053-1705fe0ac030/NUPHB11328.pdf).
Their prime-form convention gives residue $-1$; multiplying by the
constant sign $-1$ gives the normalization used here and does not
change any squared kernel.

Write the nonaxial kernel in (3) at $p$ as $X-rY$, $r\ne0$.
A local frame of $E$ is $k=X-rY$, $l=zY$. A global cubic on $E$
maps to meromorphic cubic data on $a\oplus b$ with at most a simple
pole at $p$. Its coefficient spaces are, in order,
$H^0(b(p)),H^0(a(p)),H^0(b(p)),H^0(a(p))$. The residue must be a
multiple of $(X-rY)^3$. If the residue vanished, all coefficients
would be holomorphic and hence zero. Thus every nonzero possibility
is a scalar multiple of

$$
             \Phi_r=BX^3-3rAX^2Y+3r^2BXY^2-r^3AY^3.               \tag{5}
$$

Since $k\wedge l=zX\wedge Y$, the determinant correction multiplies
the cubic by $z$. In the modified frame it is

$$
 zB k^3+3r(B-A)k^2l
       +6r^2\frac{B-A}{z}kl^2+4r^3\frac{B-A}{z^2}l^3.              \tag{6}
$$

Its extension condition is exactly $c_a(p)=c_b(p)$.
The difference $d_p=A^2-B^2$ is a holomorphic one-form: its double
pole cancels and (4) has no simple pole in either square. Its value
at $p$ is $2(c_a-c_b)$ in these frames. Thus (6) extends if and
only if $Q(p)=0$.

The one-form $d_p$ is never identically zero. Equality $A^2=B^2$
would make their nonzero meromorphic sections have the same divisor,
and their quotient would trivialize $a\otimes b^{-1}$. This would
contradict $a\ne b$.

## 3. The actual algebra is smooth

The discriminant of (5), normalized to $r=1$, is
$-108(A^2-B^2)^2$. After modification it is
$-108z^{-2}d_p^2$. Its zero divisor is consequently
$2\operatorname{div}(d_p)-2p$. To verify smoothness one must also
inspect the algebra, rather than just this discriminant divisor.

In Miranda's convention write a binary cubic as
$-\beta X^3+3\alpha X^2Y-3\delta XY^2+\gamma Y^3$.
Where $\beta$ is invertible, his multiplication formulas eliminate
one trace-free generator and give

$$
 U^3+3(\beta\delta-\alpha^2)U
                 +3\alpha\beta\delta-2\alpha^3-\beta^2\gamma=0.
$$

Substituting (5), or (6) near $p$, reduces this to

$$
                         U^3-3d_pU+2Ad_p=0.                       \tag{7}
$$

Near $p$ the modified coefficient $\beta=-zB$ is a unit.
Although $A$ has a pole, $Ad_p$ is holomorphic when $Q(p)=0$.

If $p$ is not Weierstrass, the divisor of $d_p$ is
$p+q$ with $q=\iota p\ne p$, both simple. At $p$, write
$d_p=d_1z+O(z^2)$, $d_1\ne0$. The special fiber of (7) is
$U^3+2d_1$, with three distinct roots. The cover is étale there.
At $q$, neither $A$ nor $B$ vanishes: otherwise $d_p$ would have
a double zero at $q$ as well as a zero at $p$. The unmodified
$\beta=-B$ is therefore a unit near $q$. For a coordinate $t$ at
$q$, write $d_p=\mu t+O(t^2)$. The derivative in $t$ of (7) at
its triple root is $2A(q)\mu\ne0$. This proves that the source
is smooth and the map has ramification index three.

If $p$ is Weierstrass, then $\operatorname{div}(d_p)=2p$.
For $d_p=d_2z^2+O(z^3)$, $d_2\ne0$, equation (7) has the form

$$
                      U^3+2d_2z+O(z^2)+O(z^2U)=0.                 \tag{8}
$$

It again defines a smooth source with index-three ramification,
now over $p=\iota p$. Away from the indicated point the
discriminant is nonzero. No zero-cubic fiber occurs: its leading
modified coefficient at $p$ is one, and at $q\ne p$ it is a
nonzero cube. A component containing a point of local degree three
already has the full degree of the cover, so the source is connected.
Finally $\deg T=-1$ gives genus five by Riemann--Roch.

## 4. Uniqueness, number of covers and monodromy

The parameters $r$ and the nonzero scalar in (5) are absorbed by
diagonal automorphisms of $a\oplus b$ in the cubic-algebra
correspondence. They do not give additional covers.

To show that $p$ is intrinsic, use
$T\subset a^{-1}(p)\oplus b^{-1}(p)$, whose permitted polar parts
form a line with both coordinates nonzero. After tensoring by $a$,
the first coordinate of a global section lies in
$H^0(\mathcal O_C(p))=\mathbb C$ and has no pole. The residue
condition then kills the other polar part; its remaining holomorphic
section lies in $H^0(a\otimes b^{-1})=0$. Hence
$h^0(T\otimes a)=1$, and similarly $h^0(T\otimes b)=1$.
The two sections therefore recover (2) up to their individual scalars
and recover its determinant zero $p$. This proves the bijection.

Write $\mathscr D=\sum M_{ij}\omega_i\boxtimes\omega_j$ in a basis
of $H^0(K_C)$, with $M$ symmetric. The matrix is nonsingular.
Otherwise a kernel direction is the canonical evaluation direction
at some $p$, since $C\to\mathbb P^1$ is surjective. Then
$d_p$ would be identically zero, contradicting Section 2.
The nondegenerate quadratic $v^{\mathsf T}Mv$ has two distinct
roots on $\mathbb P^1$. Their inverse images give the stated
two, three or four distinct zeros of $Q$. This also proves
existence for every fixed pair, without a genericity assumption.

The discriminant divisor is $2q$, while
$\det E=ab(-p)$. Its sign local system is therefore

$$
               \det E(-q)=ab\otimes K_C^{-1}=a\otimes b^{-1}.
$$

This nontrivial two-torsion character makes the transitive monodromy
$S_3$, rather than its cyclic subgroup. Its centralizer in the
permutation group on three sheets is trivial. Each cover consequently
has ordinary stack weight one. No local intersection or virtual
multiplicity follows from this stabilizer calculation.

## 5. Simultaneous kernels on the main genus-six boundary

Let $f:Y\to C$ be one of these covers, with triple ramification point
$r$ over $q$. Attach any pointed smooth elliptic curve $(E,e)$ at
$r$. Both sections in (1) extend over the tail by their matching
constant value. The resulting genus-six stable map has both kernels.

It lies in $M_3$, the closure of smooth four-simple-branch triple
covers. Expand the target at $q$ and use a degree-three function on
$E$ with sole pole $3e$ and four distinct simple finite branch
points. Such a function exists: the embedding by $|3e|$ realizes
$E$ as a smooth plane cubic with flex tangent $3e$; projection
from a general point of that tangent line has the required pencil.
Avoiding the other finitely many flex tangents gives simple finite
branching. A line cannot be tangent to a smooth plane cubic at two
distinct points, so their branch values are distinct. Attach its
pole to $r$, with matching index three. Fixed-target admissible-cover
smoothing as in Calculation B3 gives the desired smooth covers. The
stabilized limit is exactly the given tail map.

Properness of the gluing map includes the nodal elliptic tail.
For each cover this gives a closed copy of
$\overline{\mathcal M}_{1,1}$ in the simultaneous-kernel locus
inside $M_3$. The core cover and its attachment are intrinsic,
and its deck group is trivial, so there is no additional covering
gerbe. Different core covers give disjoint such substacks.

These maps are in the fundamental-group-surjective summand: at
degree three a nontrivial étale factorization would have degree
one upstairs, incompatible with the ramified core. Thus the
common-kernel boundary persists after the complete étale reduction
in Calculation B3. Its positive-dimensional intersections are exactly
why smooth finite-cover kernel exclusion cannot by itself evaluate
the main mixed Euler integral. The theorem makes no claim that
these elliptic-tail curves exhaust every contracted-boundary
common-kernel locus.

## 6. An exact example and scope

For $C:y^2=\prod_{j=1}^6(x-j)$, take the even theta partitions
$123\mid456$ and $124\mid356$. The hyperelliptic Szegő-square
formula gives

$$
 \mathscr D(x,t)=\frac{2xt-7x-7t+17}{y(x)y(t)}\,dx\,dt.
$$

Its diagonal zeros have coordinates

$$
                            t=\frac{7\pm\sqrt{15}}2.
$$

Neither is a branch value, and each has two lifts to $C$.
Thus this example has four covers and four corresponding
elliptic-tail boundary curves. The formula follows by subtracting
the partition numerators in
[D'Hoker--Schlotterer, Appendix C, (C.3)--(C.4)](https://scoap3-prod-backend.s3.cern.ch/media/files/69384/10.1007/JHEP12%282021%29063_a.pdf).
It is an exact algebraic calculation, not a numerical theta-function
approximation.

The exact identities are reproduced by
`code/cubic_szego_construction.py`,
using `uv run --no-project --with sympy==1.14.0 python
code/cubic_szego_construction.py`. The check verifies the frame change,
discriminant, actual algebra coefficients and explicit bidifferential;
it does not assign any kernel-locus intersection multiplicity.
