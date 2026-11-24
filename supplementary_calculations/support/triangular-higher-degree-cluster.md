# A7. Prime-degree cluster comparisons and the first composite-degree obstruction

## 1. Hypotheses and the invertible core path

Let $C$ have genus $g\ge2$, set $c=g-1$, and fix two non-effective theta
characteristics $a,b$. Thus $H^0(C,a)=H^0(C,b)=0$, and both have even
parity. The case $a=b$ is elementary; the interesting case is $a\ne b$.
Use fixed smooth lines $L_1,L_2$ of degree $c$, and the triangular equations

$$
 \bar\partial_{f^*L_1}u_1+f^*A_t(u_1)=0,
 \qquad f^*D_bu_2+f^*F_t(u_1)=0.                 \tag{1}
$$

The pullbacks include $df$ on the form factor. Require $A_t(0)=F_t(0)=0$,
fix the prescribed normal first jets along the zero section, and keep
the positive-spin endpoints fixed in parameter collars. The fixed
second operator is injective on every finite cover, including connected
stable sources, by the spin/Stokes argument and matching at nodes.

The construction of an invertible first-line path in Calculation A2 works in
every genus under the present non-effectivity hypothesis. Lower the spin
mass on $b$ to zero, join $b$ to $a$ through
$\operatorname{Pic}^{g-1}(C)\setminus\Theta$, and restore the mass.
The complement is path connected, because a path in the connected complex
torus can avoid the real-codimension-two theta divisor, including its
lower-dimensional singular strata. All the operators have index zero
and no kernel. The core sign consequently stays positive. This argument
does not justify replacing non-effectivity by even parity on an arbitrary
fixed curve.

Write

$$
 H_m=1+mc,\qquad B_m=1+\binom{m+1}{2}c.        \tag{2}
$$

The scalar-difference proof in
[the degree-two cluster investigation](triangular-cluster-transport.md)
applies in every degree: if the first-coordinate difference vanished on
an off-diagonal component, injectivity of the second operator would force
the entire difference to vanish. For a reduced simple graph of degree $m$
this gives

$$
 R+2b\le m(m-1)c,\qquad H_m\le h\le B_m       \tag{3}
$$

for a smooth embedded support of genus $h$. Here $R$ is the ramification
degree and $b$ counts nodes when a reduced nodal support is considered.
The upper bound concerns support genus, not source genera of maps with
contracted components.

For a chosen cutoff $M$, assume a common compact tame realization and
the following properties for **all actual simple curves** contained in
the chosen tube and of degrees at most $M$:

1. Every such curve is embedded, and each fixed-parameter moduli space
   is discrete in the Gromov topology.
2. Distinct such supports are disjoint. It suffices to impose this on
   pairs whose degrees sum to at most $M$.
3. There is a dense set of parameter values on which all these curves
   are unobstructed.
4. The zero-section operator is invertible throughout, and the fixed
   endpoint germs contain only the zero-section support.

These are localized versions of the cluster isolation hypotheses, not
merely surjectivity at smooth-cover linear walls. Proving them in the
restricted class (1), with the fixed core jets and endpoint collars,
is the separate projected-jet problem. A common tube isolating the whole
path is not assumed.

## 2. The degree-three comparison

Let $\mathsf L_{a,b}$ and $\mathsf L_{b,b}$ be the complete connected
local Gromov--Witten series at the endpoints. Let $G_h(q,t)$ be the
connected elementary series, whose degree-one coefficient is
$(2\sin(t/2))^{2h-2}$. Assume the hypotheses of Section 1 with $M=3$.
Then there are integers $c_{d,h}$ such that

$$
 \mathsf L_{a,b}-\mathsf L_{b,b}
 =\sum_{h=H_2}^{B_2}c_{2,h}G_h(q^2,t)
  +\sum_{h=H_3}^{B_3}c_{3,h}G_h(q^3,t)
       \pmod{q^4}.                                  \tag{4}
$$

Here is the localized refinement argument. Choose small Hausdorff
neighborhoods around the core with closures inside the tube. Bounded-energy
cycle compactness, discreteness, and countability of the embedded support
spaces give arbitrarily small radii with no cycle on the boundary.
The core is the unique nearby degree-one support. Compact open-and-closed
cycle insertions transport over short parameter intervals, and finitely
many such intervals cover the path. On each overlap choose an unobstructed
slice and refine the two neighborhoods to a common smaller one.

The difference contains no degree-one support or any multiple of one.
Disjointness excludes connected supports with several distinct irreducible
components. A remaining connected cycle of total degree at most three
is therefore a simple degree-two or degree-three curve: a multiple of a
degree-two curve already has degree at least four. For each of these
two degrees the difference is compact and discrete, hence finite. All
selected cores are actual curves for the triangular equations, so (3)
bounds their genera.

Each correction core is primitive relative to this cutoff. Its complete
cluster series is its sign times $G_h(q^d,t)$ modulo $q^4$. For this last
step the usual elementary comparison may leave the triangular class,
because every further proper-multiple correction is above the cutoff.
This proves (4). The compact insertion contains all contracted components
and all stable-map boundaries over each support.

The primary cluster inputs are Doan--Ionel--Walpuski, Theorem 2.11,
Propositions 3.10, 3.17, 3.19, 3.30, and the primitive case of Proposition
3.54. Equation (4) is our localization and degree filtering of those
arguments, not a theorem stated verbatim there.
[Primary source](https://arxiv.org/html/2103.08221v2#S3)

The endpoint identification of numerical elementary coefficients with
the universal balanced exponents is already supplied by
Article III, Proposition 3.1 and Corollary 3.2.
It follows that

$$
 n_{2,h}=0\quad(h>3g-2),\qquad
 n_{3,h}=0\quad(h>6g-5).                         \tag{5}
$$

No new evaluation on arbitrary analytic families is used. The existing
calculations also give, when $a\ne b$,

$$
\begin{aligned}
 n_{2,H_2}&=-2^{2g-1}, &
 n_{2,H_2+1}&=c\,2^{2g-3},\\
 n_{3,H_3}&=-2\left(6^{2c}-2^{2c}\right), &
 n_{3,H_3+1}&=\frac c6\left(11\cdot6^{2c}-3\cdot2^{2c}\right).
\end{aligned}                                                     \tag{6}
$$

The cubic entries use the unramified and first-excess results of
Article V, §§4–6, independently of the
conditional cutoff argument. The positive sign in these formulas is
the consequence of the non-effective endpoints.

## 3. Divisor filtering gives every prime degree

There is a stronger consequence than applying (4) only at cutoff three.
Let $p$ be prime and assume Section 1 through degree $p$. Perform the
same refinement argument, but retain only the coefficient of $q^p$.
On an unobstructed overlap slice, a connected degree-$p$ cycle has the
form $rD$ for one embedded support $D$, since distinct supports cannot
meet. If $d$ is the degree of $D$, then $rd=p$. The refinement difference
has no degree-one support, so $d=p$ and $r=1$.

Consequently this degree-$p$ difference is a compact discrete set of
simple degree-$p$ curves, all with genera between $H_p$ and $B_p$.
Their contributions at this coefficient are primitive. Thus

$$
 q^p
    =\sum_{h=H_p}^{B_p}c_{p,h}(2\sin(t/2))^{2h-2},
 \qquad n_{p,h}=0\quad(h>B_p).                   \tag{7}
$$

In the elementary logarithmic expansion, only divisors of $p$ contribute
to degree $p$. The degree-one endpoint factors agree, so the coefficients
in (7) are exactly $n_{p,h}$. No unresolved degree-four comparison enters
the argument for $p=5$, for example. Conversely, (7) is not an all-degree
statement: the prime-divisor argument is precisely what fails for four.

## 5. What changes at cutoff four

Assume now the actual-curve hypotheses through degree four. A refinement
difference still has no degree-one supports, but it can contain $2D$
where $D$ is an embedded degree-two curve. Its cluster contribution
cannot be replaced, through degree four, by its primitive term alone.

For a regular degree-two correction core $D_i$ of genus $h_i$, let
$K_i(q^2,t)$ denote its full cluster contribution. Existing elementary
decomposition gives a finite expansion through this cutoff of the form

$$
 K_i(q^2,t)=\sigma_iG_{h_i}(q^2,t)
                  +\sum_r b_{i,r}G_r(q^4,t)\pmod{q^5}.          \tag{15}
$$

The $b_{i,r}$ describe its own double-cover corrections. They are not
the double-cover coefficients of the original split-spin core.
The direct degree-four correction cores satisfy $h\le B_4$, but the
usual elementary comparison used to obtain (15) may create further
degree-four cores outside the triangular class. The original support
estimate then supplies no bound for their genera. Qualitative finiteness
is already available; the explicit bound $B_4$ is what is missing.

This is also why a prime-degree conclusion does not automatically
extend to all degrees by ordinary cluster induction. The general
comparison does not retain the local projection data used by (3).
The distinct all-genus double-wall invariant of Bai--Swaminathan is
related to elementary cluster coefficients by their Conjecture 1.10;
it cannot be substituted as an established comparison for these new
normal operators.
[Primary source](https://arxiv.org/html/2106.01206v4#S1.SS3.SSS1)

## 6. The inherited normal equations and an exact recursive bound

Suppose a degree-$m$ embedded correction $u:Y\to L_1\oplus L_2$ has
genus $h$, and its scalar projection $\bar u:Y\to L_1$ is immersed.
Then the normal operators have an exact sequence

$$
 0\longrightarrow(f^*L_2,f^*D_b)
 \longrightarrow(N_u,D_u^N)
 \longrightarrow(N_{\bar u},Q)\longrightarrow0.                 \tag{16}
$$

Let $R=2(h-1)-2mc$ be the ramification degree of $f:Y\to C$. The
degrees and real indices are

$$
\begin{array}{c|cc}
 &\deg&\operatorname{ind}_{\mathbb R}\\ \hline
 f^*L_2&mc&-R\\
 N_{\bar u}&2(h-1)-mc& R.
\end{array}                                                     \tag{17}
$$

The second line degree follows from $c_1(TL_1)[C]=-c$ and the normal
sequence of the immersed scalar curve. The first operator is injective
on every cover of $Y$, since such a cover composes with $f$ to a cover
of $C$. For $R>0$ it is not surjective on $Y$: its cokernel has real
dimension $R$. The lower off-diagonal coupling can pair this cokernel
with the positive-index kernel of $Q$ to make the full operator
invertible. Removing that coupling is not an innocuous simplification.

After a smooth complex splitting of (16), the homogeneous normal
equation is triangular. More generally, consider a deformation of this
normal equation that preserves the lower line and its injective operator,
with the first equation depending only on the quotient coordinate.
The same scalar-difference argument gives, for a reduced relative
degree-$k$ graph of genus $h'$,

$$
\begin{aligned}
 h'-1
 &\le k(h-1)+\binom{k}{2}\bigl(2(h-1)-mc\bigr)\\
 &=k^2(h-1)-\frac{mk(k-1)}2c.                                \tag{18}
\end{aligned}
$$

If $h\le B_m$, this becomes exactly

$$
 h'-1\le\frac{mk(mk+1)}2c=B_{mk}-1.                         \tag{19}
$$

Thus the quadratic genus estimate is stable under repeated passage to
these triangular normal models. No loss in the constant is forced by
the arithmetic of the recursion. In the first composite case $m=k=2$,

$$
 h'\le4h-3-2c\le1+10c=B_4.                                 \tag{20}
$$

For a genus-two base, the possible degree-two cores have $h=3$ or $4$.
Their induced relative-double bounds are respectively $7$ and $11$.
Calculation A2 proves scalar immersion for both types, so (16)--(20) have
no extra projected-immersion hypothesis in this particular first
composite case.

In larger degree or base genus, full embeddedness does not alone give
scalar immersion. At a critical point of the scalar projection the
vertical line can meet the tangent of the full curve; its map to the
normal bundle can vanish, and saturation introduces a divisor.
The inherited injective operator then cannot be identified with the
same line operator without further analysis. Consequently projected
immersion, or a substitute controlling this torsion, is an additional
hypothesis for the recursive use of (18).

## 7. Why this recursive estimate is not yet a recursive comparison

For a ramified correction, $mc<h-1$. Its preserved line in (16) is a
pullback spin line satisfying

$$
 (f^*L_2)^2=f^*K_C=K_Y(-R),                                  \tag{21}
$$

not a theta characteristic of $Y$. In particular it is topologically
different from either line of the standard equal-spin elementary model
in genus $h$.

There is a direct obstruction to ending a deformation at that standard
model while keeping a complex invariant line of degree $mc$. Write the
equal-spin operator as $D_\theta\oplus D_\theta$ on
$\theta\otimes\mathbb C^2$, with the same positive mass and metric on
the two copies. If a complex line subbundle $E$ is preserved by this
operator, apply it to arbitrary complex multiples of a local generator
of $E$. Separating the coefficients of the multiplier and its conjugate
shows that the projective line $E\otimes\theta^{-1}$ is holomorphic
and invariant under standard conjugation. It gives a holomorphic map
$Y\to\mathbb P^1$ with image in $\mathbb R\mathbb P^1$, hence is
constant. Therefore

$$
                       \deg E=\deg\theta=h-1.                \tag{22}
$$

This contradicts $\deg E=mc$ when $R>0$. The contradiction is invariant
under complex bundle identifications. It excludes this particular
endpoint while preserving the distinguished line; it does **not** prove
that no other operator in the triangular class can have elementary
contribution, nor that cancellations cannot yield the desired final
bound by another method.

A useful sufficient theorem would be an elementary comparison for the
normal extensions (16), with relative corrections selected inside
triangular models retaining the injective lower operator and the degrees
in (17), together with a suitable terminal model whose full contribution
is known. Subject to projected immersion and the requisite isolation and
genericity, (18)--(19) would then propagate $n_{m,h}=0$ for $h>B_m$
through the cluster recursion. Constructing that comparison or an
alternative replacing it is genuine additional work. The standard
primitive comparison and the existence of the normal sequence alone do
not provide it.
