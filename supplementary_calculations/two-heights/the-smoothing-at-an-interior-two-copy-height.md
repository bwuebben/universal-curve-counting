# F1. Smoothing two copies of a ramified double cover at unequal heights

Let $D_t$ be a smooth fixed-symbol path of real Cauchy–Riemann operators on a complex rank-two bundle of degree $2g-2$ over a closed curve $C$ of genus $g\ge2$. Suppose $D_0$ is invertible and $f:Y\to C$ is the transverse double-cover crossing of [Calculation E2](../polar/first-quartic-contribution-at-a-ramified-double-wall.md) with precisely two simple branch values. Its genus is $h=2g$. Its kernel generator $k$ is odd and its graph embedded. At a chosen ramification point $r$ take $w=c+z^2$, $k=vz+O(|z|^3)$ and assume $B(c)v$ is nonzero modulo $\mathbb Cv$. Here $D_0=\bar\partial+B(w)d\bar w$ in a frame holomorphic for the complex-linear part, with $B$ antilinear. Write $W_f$ for the isomorphism from the two complex branch velocities and real time to $P_f=\operatorname{coker}f^*D_0$, evaluated on $k$. Our reduced equation convention is $L K=j(e)$.

In source genus $2h=4g$ a pair of copies of $f$ has just one direct joining node, by genus additivity and stability. Write their heights as $ak,bk$ with $a,b$ real. For $|a|\ne|b|$ with both nonzero the sheets meet only at the two ramification points: at a regular value the two values of $k$ are nonzero opposites. At either ramification join the full kernel is $\mathbb R^2$. Independent lower deck involutions and interchange of the components act on its projectivization, with ratio transformations $b/a\mapsto-b/a,a/b$. The coarse quotient is $[0,1]$ by the ratio of the smaller to the larger absolute height. The projective isotropy has order two on the interior and four at either endpoint of this interval at the fixed ramification join. At ratio zero the point belongs to the marked zero-double family with upper cover $f$. At ratio one choose the two heights equal. The two flags can then be any common point of $Y$ because the graph is embedded. This is the split-nodal boundary of the two-branch double-cover family over $Y$, with the inherited section. The lower involution identifies $y$ and its conjugate, giving coarse curve $C$. Both ends of this description concern actual points of the projective incidence.

**Theorem.** On any compact subinterval of the open height interval, the two-copy family is the entire nearby zero set at its joining point, also after omitting the equation comparing the two reconstructed time coordinates. With the ordering and coordinates specified below, the projected thirteen real equations have transverse local Brouwer degree $-3$ at each fixed interior height. The nontrivial plane degree comes from a conjugate cubic in the lower flag displacement. This is a degree in directions normal to the height interval. The full absolute Euler change includes the remaining time equation and the two endpoint families.

## The fixed-target and Fredholm coordinates

Work before the finite quotients, ordering the copies and taking $a>b>0$; one may set $a=1$. Center $w$ at the old lower branch value. In the plumbing of the node take
$$
    zy=q,\qquad w=x+z^2+2pz+y^2+2sy.
$$
The two nodal-side maps have respective branch displacements $x_1=x-p^2$ and $x_2=x-s^2$. Denote the other two branch displacements by $u_1,u_2$ and put $\lambda_i=(x_i,u_i,t)$. The node-preserving constraint is exactly $x_1+p^2=x_2+s^2$. These are the two-smooth-component charts (3) of [Calculation E6](../polar/the-quartic-contribution-at-a-marked-double-cover.md), here with both normalized covers near $f$. The normalized branch deformations are unobstructed, their marked evaluations are transverse using either selected branch value, and the new coordinate $q$ smooths the node. Thus $x,u_1,u_2,p,s,q$ are six complex coordinates and $t$ is real. Use fixed outer annuli and the unique nearby lifts of the target-coordinate changes to the unramified source collars. On side one the map difference from the reference has size $O(|q|^2+|sq|)$ there, and on side two $O(|q|^2+|pq|)$, with the same bounds on first derivatives. Outside those collars the normalized maps are unaltered.

Use independent odd thick supplements $j_i:P_f\to\Omega^{0,1}$ to the two component operators $L_i$, and odd real tests $\ell_i$ on regular lower sheets, with $\ell_i(k)=1$. The map $(K_i,e_i)\mapsto (L_iK_i-j_i(e_i),\ell_i K_i)$ is invertible through this ordinary smooth family. Write $k_i$ for its normalized solution, $E_i^0$ for its obstruction. In the source coordinate centered at its moving branch point each $k_i$ is odd. One has $E_i^0=W_f\lambda_i+O(|\lambda_i|^2)$ in common frames. The matching supplement is $j_M(d)=L_1(\chi H(w)d)$ on the first side, with $H$ a real frame of $D_t$-homogeneous local sections, $H(x)=1$, and $\chi=1$ on the surgery disk. Its annular derivative is on an unchanged thick part and $\ell_1(\chi H d)=0$. All forcings and tests have support outside the inner neck and the collar identifications. The two tests prescribe $\ell_1 K=a$, $\ell_2 K=b$. Set $M=e_M$. On the nodal source the normalized solution is
$$
     V_1=ak_1+\chi H M_0,\quad V_2=bk_2,\qquad M_0=bk_2(s)-ak_1(p).
$$
Here evaluations at $p,s$ are relative to the respective branch centers. Near the moving attaching points we also write these sections as $V_1(z),V_2(y)$, with $m=V_1(0)=V_2(0)$. The full index is $-12$, and at the two-copy point the two-dimensional kernel and fourteen-dimensional cokernel are given by normalization, the latter an extension of $P_f\oplus P_f$ by the complex two-plane $N_c$. These supplements and tests therefore give a continuously invertible stabilized operator near the wall. The weighted first Sobolev norm controls the section sup norm, uniformly across the neck; smooth interior estimates apply on fixed thick disks. For the zeroth-order term and its continuity on the half-collars, use the bound on the derivatives of the holomorphic total-space map in the larger of $|z|,|y|$, as in Calculation E6. [Bai–Swaminathan, Appendix B](https://arxiv.org/html/2106.01206v4#A2) · [Siebert, Proposition 4.15](https://arxiv.org/pdf/dg-ga/9608005v2#page=52).

## The component equations and motion of the selected branch

Write $P_r(d)$ for the class of $L_f(\chi d/z)$ taken off its pole. Its representative extends at the simple ramification point as in [Calculation E5](../polar/polar-kernels-on-two-branched-double-covers.md). The selected branch derivative has the precise normalization
$$
                       W_r(\delta)=P_r(-\delta v/2).                           \tag{1}
$$
Indeed move that branch value by a real parameter times $\delta$, holding the other branch and time fixed and choosing the supplement on unchanged thick sheets. Differentiate the normalized stabilized solution at fixed target value off the branch. In the inner chart its argument is the nearby root of $z^2-\epsilon\delta$. Differentiating at zero gives $-\delta v/(2z)$ plus a regular Sobolev remainder: $k-vz$ and its first two derivatives are respectively $O(|z|^3),O(|z|^2),O(|z|)$, and the derivative of the smooth family at fixed source coordinate is smooth. The remainder is $O(1)$ with bounded first weak derivatives. Away from the mark this exterior derivative solves the original operator with forcing $j(W_r\delta)$. Subtracting the cutoff pole gives a global Sobolev section, so passage to the cokernel proves (1). This is a real-linear identity; the polar class is entirely in the selected complex branch-velocity plane in $W_f$ coordinates.

The two actual normalized component rows satisfy, as all parameters tend to the wall,
$$
 e_1=aE_1^0-2b W_r(q)+o(|q|),\qquad
 e_2=bE_2^0-2a W_r(q)+o(|q|).                                                  \tag{2}
$$
Use the upper Green calculation of Calculation E6 once on each half. The rank of $L_i$ can jump here, so use the dual of the normalized component inverse above. A dual test for $e_i$ has smooth adjoint equation with a multiple of $\ell_i$ as its only interior source, and its pairing with $j_i$ is the desired coordinate. It annihilates $j_M$ on the first side since $\ell_1(\chi H d)=0$ and the test is smooth across the mark. The normalization term is the same as on the original nodal solution and cancels in the difference. The tests vary smoothly with $\lambda_i$ and are homogeneous on the joining disk. Truncating the other side at a fixed small source radius $\rho$, the boundary term pairs the test at the mark with $q$ times the opposite section's integral $\oint K(\zeta)d\zeta/(2\pi i\zeta^2)$ on $|\zeta|=\rho$. This integral tends to the opposite height times $v+O(\rho)$. The map differences on the receiving side are $q^2/\zeta^2+2p_j q/\zeta$, where $p_j=p$ or $s$. The derivative of the first term has integral $O(\rho|q|)$. The possible $|p_j q|\log(1/|q|)$ term in the second is $O(|p_jq|)$: freeze the smooth coefficient and use the two nonzero angular modes of order two of $K$. As in Calculation E6 their integrated absolute values in $d\log|\zeta|$ are bounded by the section sup norm and the integrated forcing in the radial mode ODE, bounded by $C\int(r+|q|^2/r^2+|p_jq|/r)\,d\log r$ on $|q|/\rho<r<c_1$. The thick collar difference is $O(|q|^2+|p_jq|)$. Letting $\rho$ decrease gives the two incoming pole classes with positive sign in the boundary convention $P_r=[L_f(\chi d/z)]$: the exterior derivative of a section with prescribed pole contributes that off-pole class to $e_i$. Equation (1) proves (2). All estimates are uniform along sequences of bounded heights tending to a fixed pair.

Use $W_f^{-1}$ on $e_1/a,e_2/b$. Retain both complex branch coordinates from each copy and the first real time coordinate. These nine real output coordinates have leading expressions
$$
  x-p^2-2(b/a)q,\quad x-s^2-2(a/b)q,\quad u_1,\quad u_2,\quad t,
                                                                                \tag{3}
$$
with errors $O(\max_i|\lambda_i|^2)+o(|q|)$. The unused coordinate can be taken as the second time minus the first. In particular its linear smoothing term vanishes. If the nine retained rows vanish, the original and any straight interpolation to (3) satisfy $\max_i|\lambda_i|=O(|q|)$: each individual branch displacement and the shared time are bounded by $C|q|+C\max_i|\lambda_i|^2$. Subtraction now gives
$$
             s^2-p^2=2(b/a-a/b)q+o(|q|).                                     \tag{4}
$$
On an interior compact height interval this implies $|q|=O(\varepsilon^2)$ with $\varepsilon=\max(|p|,|s|)$.

## A comparison section on the neck

In the common smooth target frame put
$$
                         V=V_1(z)+V_2(y)-m.
$$
On each outer collar cut off the term from the opposite side minus $m$ and identify with the original section outside. The resulting global section has prescribed scalar normalizations and the nodal supplement coefficients. On the first half $|z|\ge\sqrt{|q|}$, write $y=q/z$. In the homogeneous collar equations the cross terms from the zeroth-order operator are bounded by
$$
 C\left[|y|(|z|+|p|)+|z|(|y|+|s|)|q|/|z|^2\right].
$$
Indeed $|V_j(\zeta)-m|\le C|\zeta|$, and $|dw_1/dz|\le C(|z|+|p|)$, $|dw_2/dz|\le C(|y|+|s|)|q|/|z|^2$. Replacing the target coefficients by their values at the two reference maps has the same or smaller bound: their respective displacements are $O(|y|^2+|s||y|)$ and $O(|z|^2+|p||z|)$, and the sections are bounded. The corresponding estimate holds on the second half. The outer cutoff derivatives are $O(|q|)$ and the change of map and collar identification costs $O(|q|^2+\varepsilon|q|)$. Hence the residual is bounded by $C(|q|+\varepsilon\sqrt{|q|})$ also in the weighted target norm (the half-collar area has bounded mass). The uniform inverse gives
$$
          \|K-V\|_\infty+|M-M_0|\le C(|q|+\varepsilon\sqrt{|q|}).             \tag{5}
$$
In particular on the zeros of (3) up to the indicated errors the matching row to order one is $v(bs-ap)+O(\varepsilon^2)$. Choose $N_c=\mathbb Cv\oplus\mathbb C B(c)v$. If its first plane row vanishes, or is linearly interpolated to $v(bs-ap)$ and vanishes, then $s=(a/b)p+O(\varepsilon^2)$. Consequently either all normal variables vanish or, along any sequence of candidate zeros,
$$
             s/p\longrightarrow a/b,\qquad
             q/p^2\longrightarrow -a/(2b) .                                  \tag{6}
$$
When heights vary these statements mean convergence of the differences to zero. Thus the calculation still required is on a cone $R_0=\sqrt{|q|}\asymp\varepsilon$, $\lambda_i=O(\varepsilon^2)$; flatter approach to $q=0$ has already been excluded.

## The matching logarithm

We compute $M-M_0$ to error $O(\varepsilon^3)$ on the comparable cone. Let $\Lambda$ be a complex-linear covector on the local rank-two bundle and use the real functional $\operatorname{Re}\Lambda$ at the common marked value. Define the evaluation test on each nodal normalization using the dual of the normalized inverse: evaluate the solution of $(L_i V-j_i(e),\ell_iV)=(g,h)$ at the mark. Its forcing functional annihilates $j_i$ (for such input take $V=0$), and its scalar functional has value $\operatorname{Re}\Lambda(k_i(p_i))$ on $1$, where $p_1=p,p_2=s$. In our trace convention the forcing functional is $-\operatorname{Re}\int h_i g\,dA/\pi$ in local coordinates, with $h_i(z)dz$ a differential with leading pole $\Lambda\,dz/z$. The normalization follows by applying the functional to $L_i(\chi H d)$, for which it takes the value $\operatorname{Re}\Lambda(d)$. Equivalently the small clockwise boundary in the punctured component contributes $-\pi\operatorname{Re}\Lambda(d)$ to the area form version of Green's identity. Its only additional adjoint source is proportional to the thick test $\ell_i$. These assertions follow first by duality on $W^{1,m}$, $m>2$, where point evaluation is continuous; the distributional adjoint equation and the Cauchy kernel give the displayed pole.

Choose a frame holomorphic for the complex-linear part of $D_t$ on the disk, smoothly in $t$. If $H_\Lambda(d)=\overline{\Lambda(B_t(x)d)}$, the local adjoint expansion, with $\zeta$ centered at the flag, is
$$
 h_i(\zeta)=\Lambda/\zeta+
       2(\zeta+p_i)H_\Lambda\log|\zeta|^2+U_i(\zeta),\qquad
 |U_i(\zeta)|\le C(|\zeta|+\varepsilon)
       \quad (R_0\le|\zeta|\le\rho).                                         \tag{7}
$$
Here $\rho$ is any fixed sufficiently small radius. To verify it, the formal *real* adjoint of $L_iV=\partial_{\bar\zeta}V+B_t(w_i)V\,\overline{w_i'}$ gives
$\partial_{\bar\zeta} h_i(d)=w_i'\overline{h_i(B_t(w_i)d)}$ off the mark in this frame. Thus the singular source beyond the pole is $2(\zeta+p_i)H_\Lambda/\bar\zeta$. Subtract the displayed logarithmic primitive, with local cutoffs for it and the pole. The inhomogeneous remainder is bounded by $C(|\zeta|+\varepsilon)^2(1+|\log|\zeta||)$ on the disk, in addition to the bounded adjoint multiplier times $U_i$. This uses $|w_i-x|\le C(|\zeta|^2+|p_i||\zeta|)$. The regular remainder has locally uniform $W^{1,m}$ norm for finite $m>2$: the cut-off singular parts now have uniformly $L^m$ ordinary forcing, and one applies the normalized smooth-component adjoint inverse, including its finite-dimensional constraints. The data also vary smoothly in these norms when the branch and mark move in coordinates centered at the mark. In particular $U_i(0)=O(\varepsilon)$ on our cone. At $p_i=0$ the evaluation is invariant under the component deck, the normalized supplements and tests are odd, and the coefficient of its invariant differential is odd; both singular terms have odd coefficients as well. This gives $U_i(0)=0$ at that mark. Smooth dependence of the subtracted forcing in $L^m$ controls the displacement: after differentiation the ratios $(w_i-x)/\bar\zeta$ which occur in the regular forcing are uniformly bounded, as are their parameter derivatives on the fixed disk.

For a quantitative modulus use the Cauchy integral for $U_i(\zeta)-U_i(0)$, separating the disk of radius $2|\zeta|$ from its complement. With $r'=|\zeta|$ the complement contributes at most $C r'[1+\varepsilon |\log r'|+\varepsilon^2|\log r'|^2]$: the kernel difference is at most $C r'/|y|^2$ there, and the forcing is at most $C[(|y|+\varepsilon)+(|y|+\varepsilon)^2(1+|\log|y||)]$. On the inner disk the two single Cauchy kernels give the same bound by integration (or rescaling to radius two). The smooth boundary term is $O(r')$. For $r'\ge R_0\asymp\varepsilon$ this proves (7). A Hölder bound $O(r'^\alpha)$ with $\alpha<1$ alone would give a weaker order in the boundary pairing.

For $r'=|\zeta|$ write
$$
       T_i(r')=\operatorname{Re}\operatorname{avg}_{|\zeta|=r'}
                                      [\zeta h_i(\zeta)K(\zeta)] .
$$
The half-circles at $r'=R_0$ describe the same seam on the smoothed surface, with opposite angular parametrizations. On each exterior use the reference component operator $L_i$. The thick forcing for $e_i$ pairs to zero by construction; the scalar normalization and, on the first component, the matching supplement have exactly their nodal values as functionals. Green's identity then yields
$$
 \operatorname{Re}\Lambda(M-M_0)
  =T_1(R_0)-T_2(R_0)
   -\frac1\pi\operatorname{Re}\left(\int_1 h_1(L-L_1)K\,dA
                                   -\int_2 h_2(L-L_2)K\,dA\right).             \tag{8}
$$
Here $L$ is the actual pullback and each integral is on the corresponding half with its outer collars identified. The notation includes the smooth coordinate pullbacks there. In fact truncating a component yields $T_i$ equal to its scalar-normalization functional plus its forcing functional for $LK$ plus $\operatorname{Re}\int_i h_i(L-L_i)K/\pi$. On subtraction the scalar terms give $-\operatorname{Re}\Lambda(M_0)$ and the $LK$ terms give $\operatorname{Re}\Lambda(M)$. This proves (8) and fixes its signs. In particular a possible logarithmic integration constant in one component has not been specified by an assumed holomorphic neck expansion; (8) compares both constants in the actual normalized solution.

Put $d_0=q(bp-as)$ and $\beta=\Lambda(B(c)v)$. On the comparable cone (5) gives, at the seam,
$$
          K=v(az+by)+m+O(\varepsilon^2),\qquad |m|=O(\varepsilon).
$$
Indeed the branch-centered $k_i$ have no quadratic term by oddness, their linear coefficients are $v+O(|\lambda_i|)$, and $H(w_1)-H(x)=O(|z|^2+|p||z|)$ near the mark. The leading terms from $\Lambda/\zeta$ in $T_1-T_2$ cancel identically for any single-valued $K$. In the next term $\zeta(\zeta+p_i)H_\Lambda(K)$ the constant $m$ has no zeroth angular mode after multiplication. The only constant products are $p z\,bq/z$ on the first side and $s y\,aq/y$ on the second. At the seam the $U_i$ terms are $O(\varepsilon^3)$, and replacing $K$ in the logarithmic term by its displayed approximation costs $O(\varepsilon^4|\log\varepsilon|)$. Replacing $B_t(x)$ by $B(c)$ costs less, since $t,x=O(\varepsilon^2)$. Hence
$$
          T_1(R_0)-T_2(R_0)=
                 -4\operatorname{Re}(\bar\beta d_0)\log(1/R_0)
                 +O(\varepsilon^3).                                           \tag{9}
$$

On the first half the target-map difference is $G=q^2/z^2+2sq/z$. The logarithmic part of the first area term with its minus sign in (8) comes from $-2sq/z^2$ in $G'$, $\Lambda/z$ in $h_1$ and $avz$ in $K$. In polar coordinates $-(1/\pi)\int$ is $-2\int \operatorname{avg}[r'^2\cdot]\,d\log r'$; this term therefore gives $+4\operatorname{Re}(\beta a\overline{sq})\log(1/R_0)$. On the second side the corresponding term before subtraction is $+4\operatorname{Re}(\beta b\overline{pq})\log(1/R_0)$. For the remainder use, on $R_0\le r'\le\rho$,
$$
  |K|\le C r',\quad
  |K-avz|\le C(\varepsilon+|q|/r'+r'^3+\varepsilon^2),\quad
  |h_1-\Lambda/z|\le C(r'+\varepsilon)(1+|\log r'|).
$$
The first two estimates follow from (5) and the same Taylor expansion on each half; the second half interchanges the indices. The quadratic pole derivative contributes at most $C|q|^2/R_0=O(\varepsilon^3)$. The replacement of $K$ in the simple-pole derivative costs at most
$C|sq|\int_{R_0}^\rho(\varepsilon+|q|/r'+r'^3+\varepsilon^2)dr'/r'^2=O(\varepsilon^3)$.
The regular adjoint terms have an integrable radial bound after extracting $|sq|$. Variation of $B_t$ between $w$ and $w_1$ is $O(|q|^2/r'^2+|sq|/r')$; multiplied by $|w_1'|\le C(r'+\varepsilon)$ it contributes $O(\varepsilon^3)$. Freezing the smooth factor in the simple-pole term at the central value has the same bound, using $|w-c|\le C(r'^2+\varepsilon r'+\varepsilon^2)$. The fixed collars cost $O(|q|^2+\varepsilon|q|)$ even in the differential-part comparison by uniform thick elliptic estimates. This proves that the signed area terms together give $-4\operatorname{Re}(\bar\beta d_0)\log(1/R_0)+O(\varepsilon^3)$. They have the same sign as the seam term (9).

As $R_0^2=|q|$, the full vector matching equation on this cone is
$$
                M=M_0-4B(c)v\,\overline{q(bp-as)}\log(1/|q|)
                                   +O(\varepsilon^3).                          \tag{10}
$$
This coefficient agrees with $-2B(c)v$ times the conjugate residue of $(az+by+m')dw$ times the logarithm, for any constant $m'$. The residue is $2q(bp-as)$. The simple calculation of that residue alone would not control the two normalized evaluation constants in (8).

## Isolation and the projected degree

With $\lambda_i=O(\varepsilon^2)$, the component of $M_0$ modulo $\mathbb Cv$ is $O(\varepsilon^3)$ by oddness and smooth dependence of the branch derivative. On a possible zero of the retained nine rows and of the $\mathbb Cv$ row, (6) gives
$$
              \frac{q(bp-as)}{p^3}+\frac{a(b^2-a^2)}{2b^2}\ \longrightarrow\ 0.
$$
The limit scalar for the residue is nonzero on every compact interior interval. Thus (10) excludes the vanishing of the remaining complex matching quotient on such a sequence. If $p=s=0$, (4) implies $q=0$ locally and the retained branch/time equations give $\lambda_i=0$. This proves isolation transversely in the projected system and in the full system.

For the degree use ordered representatives $a>b>0$, complex domain orientation on $(u_1,u_2,x,q,s,p)$ and increasing $t$. In the target order the other two branch rows, the first time row, the two selected branch rows, and the matching planes in the frames $v$ and $B(c)v$. Each normalized component's coordinates are supplied by $W_f^{-1}$ after division by its height. The omitted real row is the second time minus the first. This specifies the relative-time coorientation for the projected statement.

First interpolate the retained nine coordinates to (3), and the first matching coordinate to $bs-ap$. Keep the matching quotient unchanged. Any sequence of noncentral zeros of this homotopy tending to the height interval would obey (4)–(6), as proved above uniformly in the interpolation parameters. Equation (10) rules it out. Next interpolate the matching quotient to
$$
      \frac{4a(b^2-a^2)}{b^2}\,\bar p^{\,3}\log(1/|p|)\,[B(c)v],
$$
extended by zero at $p=0$. On the zeros of the other interpolated rows, now exact, its ratio to the actual quotient tends to one by (6), (10) and $\log(1/|q|)/\log(1/|p|)\to2$. Equivalently their difference is $o(|p|^3\log(1/|p|))$ along that set. This interpolation too has a zero-free normal link for sufficiently small radius. Finally replace the positive logarithm by one. The first block is triangular: the two selected branch equations solve complex linearly for $(x,q)$ with determinant $2(b/a-a/b)\ne0$ after the terms in $p,s$ have been specified. The first matching row solves complex linearly for $s$ at fixed $p$. The independent other branch rows and increasing time have positive sign. All these blocks have local real degree $+1$, and the nonzero scalar times $\bar p^3$ has degree $-3$. This proves the stated projected degree and explains why its computation uses no estimate for the discarded relative-time row beyond (2). Reversing an output coorientation would change a projected signed degree.

At the equal-height end the determinant just used vanishes and the split tower boundary moves; at the zero-height end division by the second amplitude fails and the marked zero-double family enters. The compact support is accordingly not the disjoint union of the two endpoint families and separately countable interior heights. Neither a division of this projected degree by a point stabilizer nor an integration of it over the open interval yields the full quartic Euler change.
