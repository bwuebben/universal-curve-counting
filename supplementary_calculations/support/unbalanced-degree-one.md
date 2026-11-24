# A10. The degree-one class for two line bundles of unequal degrees

Let $C$ be a smooth compact complex curve of genus $g\ge2$, let $A B=K_C$ be a pair of holomorphic lines with $\deg A=g-1+c$ and $c\ge0$, and give the two summands weights $r,s\in\mathbb Z\setminus\{0\}$, respectively. Write $X_1(A,B;r,s)$ for the degree-one insertion in the rational complex cycle algebra $H$, normalized by the signed weights of the finite Borel tests. Then
$$
             X_1(A,B;r,s)=(s/r)^c\,x_{g,1}.
$$
In particular the full degree-one antidiagonal class is $(-1)^c x_{g,1}$. The calculation uses compactly supported relations on spaces of complex cycles. For an oriented Cayley bundle in the positive maximal flat normal model, its image under $\Phi$ gives the degree-one coefficient in the smooth universal series by Article III, Appendix E. For the equivariant insertion, see [the antidiagonal local-curve construction, §2](https://doi.org/10.5281/zenodo.24432700).

## Reduction and a two-circle proof

Deform $A$ in its connected Picard variety, keeping $B=K_CA^{-1}$. Degree-one cycles are sections: their sole component maps isomorphically to $C$ on its normalization and affine fibers carry no compact curves. For the deformation choose sufficiently many fixed distinct base points on each summand to detect every section throughout the path (more than its nonnegative degree suffices). Their common zero is the compact family of core sections over the finite Borel base and path. The compact prism relation and stabilization by appending tests give deformation invariance, including at a jump of $h^0$. Paths can be taken piecewise real analytic in the Picard variety with a Poincaré bundle and locally trivialized evaluation lines. A general $B$ of degree $g-1-c$ has $h^0(B)=0$, since the effective locus has dimension at most $g-1-c<g$ when the degree is nonnegative. Serre duality and Riemann–Roch now give
$$
 H^1(A)=0,\qquad \dim H^0(A)=\dim H^1(B)=c.
$$
Choose $c$ independent fiber evaluations giving coordinates on $H^0(A)$. At this generic splitting these evaluations by themselves cut out the zero section, so the normalized insertion may be computed with this shorter tuple. For $c=0$ the empty tuple already gives $x_{g,1}$.

Put $z=r/s$. On the two-circle base $\mathbb P^i\times\mathbb P^{c-i}$ let $Y_i$ be the insertion with unit weights on the respective axes and the same $c$ tests on $A$. The proof over the product of projective spaces in Article II, §2 uses only the number of separating tests and their characters, hence gives $X_1(z)=\sum_{i=0}^c z^{i-c}Y_i$ with these degrees as well. For $i>0$ the second projective factor has dimension less than $c$. Choose $c$ sections $v_j$ of its hyperplane line with no common zero and choose a basis $e_j$ of $H^1(B)$, represented on a two-set Stein cover of $C$. Deform the second line to the additive torsor whose transition adds $t\sum_jv_je_j$, independent of the first coordinate. This is an ordinary holomorphic family, with complex $t$ in a neighborhood of $[0,1]$. At $t\ne0$ its restriction to the graph of every first-coordinate section has nonzero cohomology class, so there is no degree-one cycle. At zero the test zero set is the family of core sections over the compact base. Keeping the first projection and its tests on the prism proves $Y_i=0$. At $z=1$, scalar normalization for the bundle of degree $2g-2$ gives $Y_0=x_{g,1}$. The formula follows for all signed weights by the two-circle identity of Article II, §2. For the scalar normalization, see [Pardon, Remark 4.1 and Definition 4.7](https://arxiv.org/pdf/2308.02948v3).

## Direct affine deformation and signs

The single-circle model has base $P=\mathbb P^c$, first line $A\boxtimes\mathcal O_P(r)$ and second line $B\boxtimes\mathcal O_P(s)$. Write $V=H^0(A)$ and $W=H^1(B)=V^\vee$. Its tests are the $c$ coordinates of the first section $u$, taking values in $V\otimes\mathcal O_P(r)$, with normalization $1/r^c$. Choose $\delta\in H^1(BA^{-1})=H^0(A^2)^\vee$ such that
$$
 D_\delta:V\longrightarrow W,\qquad a\longmapsto a\delta
$$
is invertible. Indeed sum the point evaluations of sections of $A^2$ in squared line frames at $c$ ordinary points where evaluation on $V$ is an isomorphism. This produces the symmetric bilinear form of matrix $E^{\mathsf T}E$, for $E$ the invertible evaluation matrix.

Write $n=|s|$. If $s>0$, take $\sigma_j=Z_j^n-Z_0^n$, $1\le j\le c$, as sections of $\mathcal O_P(s)$. If $s<0$, take their conjugate-linear metric duals with the Fubini–Study metric. Both systems have precisely $n^c$ simple common zeros, all in $Z_0\ne0$. At each zero the tuple has local real degree $\operatorname{sgn}(s)^c$. Indeed its untwisted holomorphic derivative is invertible, and metric duality on $c$ complex lines has real sign $(-1)^c$. Derivatives of the metric multiply section values that vanish at the point. Choose a real-analytic section $\tau$ of $\mathcal O_P(s-r)$ nonzero in this chart, using a power of $Z_0$, its metric dual or the constant section according to the sign of $s-r$.

Use a disc and the complement of its center as a two-set Stein cover of $C$, with Stein intersection, and fix overlap representatives of the $e_j$ and $\delta$. On the total space of the first line, with its tautological coordinate also denoted by $u$, deform the second line to an affine line bundle by adding to its transition
$$
         t\left(\sum_{j=1}^c\sigma_j e_j+\tau u\delta\right).
$$
Both terms take values in the pulled-back $B\boxtimes\mathcal O_P(s)$. They are holomorphic in the curve and fiber coordinates and real analytic on $P\times[0,1]$. This gives a family in the complex analytic diagram over real-analytic base simplices. Explicitly the complexification of the base contains the open neighborhood of the real diagonal in $P\times\overline P$ where $\sum Z_i W_i\ne0$, with independent homogeneous coordinates $Z,W$. There the metric dual of $F\in H^0(\mathcal O(n))$ extends as $\overline F(W)Z^{\otimes n}/(\sum Z_iW_i)^n$ in the tautological power of the first factor; $\overline F$ here conjugates the polynomial coefficients. The same formula applies to the negative-power choice of $\tau$. The gluing is holomorphic on this neighborhood and in complex $t$. Thus signs from conjugate sections are compatible with the analytic simplicial prism, including when $r$ and $s$ have opposite signs.

The degree-one cycles in each threefold are sections over $C$. They project to $u\in V\otimes\mathcal O_P(r)$ and lift exactly when the torsor class
$$
         t\left(\sum_j\sigma_j e_j+\tau D_\delta u\right)
          \in W\otimes\mathcal O_P(s)
$$
vanishes. A lift is then unique by $H^0(B)=0$. The test imposes $u=0$. At $t=0$ its zero set is all of $P$ and at $t>0$ its support is precisely the $n^c$ constant sections of the parameter interval corresponding to the zeros of $\sigma$. At each such parameter the lift of $u=0$ is visibly the zero section in both local charts. This describes a compact supported prism. Near the endpoint zeros its actual cycle space is the graph
$$
           u=-\tau^{-1}D_\delta^{-1}\Big(\sum_j\sigma_j e_j\Big)
$$
over $P$. For continuity of the graph parametrization, restrict the overlap translation to the graph and use its expression in the finite span of the $e_j$ and of $a\delta$ for $a\in V$. Choose, once, local cochains bounding a basis of the subspace of coboundaries in that span. In local parameter frames the resulting lift varies real analytically. The continuous projection to first-coordinate sections gives the inverse in the topology of cycles.

The local degree of the $u$ tests on this graph is $\operatorname{sgn}(s)^c$ at each zero: the remaining factors of their derivative are complex isomorphisms. In the oriented base graph their Thom class is this multiple of the base Gysin of the unit in the marked fiber, as in [the geometric product formula, §6](https://doi.org/10.5281/zenodo.24455518). At a zero the target itself is the total space of the extension of $A$ by $B$ with connecting homomorphism $\pm\tau D_\delta$, in the fixed parameter frames. This follows directly from the linear overlap term, as the constant translations vanish there. The cohomology sequence gives $H^0(N)=H^1(N)=0$ for the extension $N$. It has degree $2g-2$, and its only degree-one cycle is its zero section. Thus its unit in $H$ is $x_{g,1}$ by the empty-tuple scalar construction. The total normalized contribution is $n^c\operatorname{sgn}(s)^c/r^c=(s/r)^c$ times this class. This proves the full class identity in degree one. Higher multiplicities require their own geometric calculation.
