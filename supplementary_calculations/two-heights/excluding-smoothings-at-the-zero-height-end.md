# F5. The zero-height end of the two-branch quartic wall

Let $f:Y\to C$ be the transverse odd double with two branch values of [Calculation F1](the-smoothing-at-an-interior-two-copy-height.md). Assume the base operator is invertible, the graph of its kernel generator $k$ is embedded, and $[B(c)v]\ne0$ in $N_c/\mathbb Cv$ at the ramification join under consideration, where $v=k_z(0)$ in $w=c+z^2$. Use the height perturbation of [Calculation F4](separating-the-two-heights-on-a-common-parameter.md). Near its zero-height endpoint the two effective times are constant shifts on the respective pairs of sheets over the support of the perturbing tensor. They may be treated as independent small inputs $t_1,t_2$ for the local calculation.

**Proposition.** A smoothed joined map with a projective kernel cannot tend to the zero-height endpoint at which the upper cover is $f$, even when the two effective times vary independently and tend to zero. At every sufficiently small nonzero perturbation parameter the entire nearby zero set consists of the zero-height marked family, with its joining node retained.

## The two length scales

Use the ordered chart $zy=q$, $w=x+z^2+2pz+y^2+2sy$, centered at the old joining value, and the two real thick normalization tests $a,b$ from Calculation F1. In a projective chart take $a$ bounded away from zero and bounded above, and $b\to0$. Write $\lambda_1=(x-p^2,u_1,t_1)$ and $\lambda_2=(x-s^2,u_2,t_2)$, where $u_i$ is the other branch displacement. All these parameters approach zero. The supports of the two time shifts are disjoint from the surgery disk. The smooth-component estimates, stabilized inverse and Green formulas in [Calculation F1](the-smoothing-at-an-interior-two-copy-height.md) therefore apply with the respective times. In particular the *unscaled* component rows on a zero satisfy
$$
 a\big(W_f\lambda_1+O(|\lambda_1|^2)\big)=2bW_r(q)+o(|q|),\qquad
 b\big(W_f\lambda_2+O(|\lambda_2|^2)\big)=2aW_r(q)+o(|q|).
$$
Their remainders are uniform along bounded heights tending to a fixed pair, including $(a,0)$ after passage to a subsequence. This is the estimate on smooth tests before either amplitude is divided out.

Suppose $q\ne0$. The second row excludes $b=0$ since $W_r$ is injective. Inverting the component derivative yields $\lambda_1=o(|q|)$ and $|b||\lambda_2|=O(|q|)$. The Taylor remainder $b\,O(|\lambda_2|^2)$ is therefore $o(|q|)$. The selected branch rows give
$$
                   2a q=b(p^2-s^2)+o(|q|).                                    \tag{1}
$$
Put $\varepsilon=\max(|p|,|s|)$ and $\rho_0=\sqrt{|q|}$. In particular $|q|\le C|b|\varepsilon^2$. Let $M$ be the actual matching row and $M_0=b k_2(s)-a k_1(p)$ its value on the nodal problem with the same component parameters. The glued comparison gives $|M-M_0|\le C(|q|+\varepsilon\rho_0)$. By oddness $k_i(j)=vj+O(|\lambda_i||j|+|j|^3)$ for the respective branch-centered argument $j$. Using $M=0$ we get
$$
              |p|\le C\big(|b||s|+|q|+\varepsilon\rho_0\big).
$$
Since $\rho_0\le C\sqrt{|b|}\varepsilon$, this forces $p=o(\varepsilon)$, and hence $\varepsilon=|s|$ on the sequence. Returning to (1) gives
$$
  q/(b s^2)+1/(2a)\longrightarrow0,\qquad
  \rho_0\asymp\sqrt{|b|}|s|,\qquad p=o(\rho_0),\qquad
  |\lambda_2|=O(|s|^2).                                                       \tag{2}
$$
For the third relation divide the preceding inequality for $|p|$ by $\rho_0$ only after the second relation is known. In particular $|b s|=o(\rho_0)$. No comparison of the logarithms of $|b|$ and $|s|$ has been made.

## The matching term

We retain the seam at radius $\rho_0$ on both normalizations. Let $K$ be the true normalized section, and $m=b k_2(s)$ the common central value of the nodal comparison section. On this seam
$$
                 K=m+v(az+by)+o(\rho_0),\qquad |K|\le C\rho_0.
$$
Indeed (5) of Calculation F1 bounds the discrepancy from the glued nodal section by $O(|q|+|s|\rho_0)=o(\rho_0)$. In that glued section the smooth odd component derivatives converge at the joining point to the indicated linear terms. The matching cutoff on the first side equals a smooth local homogeneous frame times $M_0$ and changes from its marked value by $O(|M_0|(|z|^2+|pz|))$. On either half annulus of fixed outer radius, writing $r'$ for the larger source radius, the same comparison gives $|K|\le C r'$ since $|m|=O(|bs|)\le C\rho_0$. On the first half its more precise form is
$$
        |K-avz|\le C(\rho_0+r'^3)+o(1)r'\qquad(\rho_0\le r'\le\rho).
$$
Here the opposite section varies by $O(|bq|/r')$, and the Taylor error on the first component after its marked value is removed is $O((|\lambda_1|+|p|^2)r'+r'^3)$. These estimates control the section on the full logarithmic annulus.

For an evaluation test with polar coefficient $\Lambda/z$ on the first component, or $\Lambda/y$ on the second, use expansion (7) of Calculation F1. Its regular remainder satisfies
$$
              |U_i(\zeta)|\le C(|\zeta|+|p_i|),
              \qquad p_1=p,\quad p_2=s.
$$
This bound extends to radii much smaller than $|s|$. In fact at $p_i=0$ the regular constant is zero by deck symmetry for every nearby component parameter. The $L^m$ differentiated cutoff equation gives $U_i(0)=O(|p_i|)$. The Cauchy integral bound there gives the modulus
$C r'(1+|p_i||\log r'|+|p_i|^2|\log r'|^2)$ with uniform constants. Since $r'|\log r'|^j$ is bounded for $j=1,2$, this implies the displayed estimate throughout the annulus. In particular the full adjoint $h_i$ has $|h_i|\le C/r'$.

Set $d_0=q(bp-as)$ and $\beta=\Lambda(B(c)v)$. In the seam term of the signed Green identity (8) of Calculation F1 the simple evaluation poles cancel exactly for the actual section. The logarithmic terms have circular average $-4\operatorname{Re}(\bar\beta d_0)\log(1/\rho_0)$, up to $o(|qs|\log(1/\rho_0))$. The constant $m$ has zero circular average after multiplication by $\zeta(\zeta+p_i)$. The seam error $o(\rho_0)$ is multiplied by $O(\rho_0(\rho_0+|s|)|\log\rho_0|)$ and the regular terms contribute $O(\rho_0^2(\rho_0+|s|))=O(|qs|)$.

In the first signed area integral the derivative of $2sq/z$ and the terms $\Lambda/z,avz$ give $+4\operatorname{Re}(\beta a\overline{sq})\log(1/\rho_0)$. Replacing $avz$ by the true section in that leading term costs
$$
 C|sq|\int_{\rho_0}^{\rho}
      \big(\rho_0+r'^3+o(1)r'\big)\frac{dr'}{r'^2}
   =O(|sq|)+o(|sq|\log(1/\rho_0)).
$$
The regular and logarithmic adjoint terms left after the simple pole have absolute value at most $C(r'+|p_i|)(1+|\log r'|)$, integrable in $dr'$ against the bounded factor $|K|/r'$. The second simple-pole derivative $-2pq/y^2$ costs at most $C|pq|\log(1/\rho_0)=o(|sq|\log(1/\rho_0))$ by (2). The two quadratic-pole derivatives cost $O(|q|^2/\rho_0)=O(|q|\rho_0)$. For a changed coefficient the displacement on a half annulus is $O(|q|^2/r'^2+|p_j q|/r')$ and the reference map derivative is $O(r'+|p_i|)$. Integrating with $|h_iK|\le C$ gives $O(|qs|+|q|^2(1+|s||\log\rho_0|))$. Freezing the smooth coefficient in the first principal integral gives $o(|qs|\log(1/\rho_0))+O(|qs|)$ by $w\to c$ on shrinking inner disks and a bounded integral on their complements. The outer collar errors are $O(|q|^2+|s q|)$. Thus the seam and first area terms have the same leading sign, as $bp=o(s)$.

In vector form the resulting equation is
$$
       M=M_0-4 B(c)v\,\overline{q(bp-as)}\log(1/|q|)
                       +o(|qs|\log(1/|q|)).                                  \tag{3}
$$
The image of $M_0$ in $N_c/\mathbb Cv$ is
$O(|p|^3+|p||\lambda_1|+|b||s|^3+|b s||\lambda_2|)=O(|qs|)$ by (2). The logarithmic coefficient in this quotient is bounded below in norm by a positive multiple of $|qs|$, since $a$ stays away from zero and $[B(c)v]\ne0$. This contradicts $M=0$, proving the exclusion.

At $q=0$ the two component rows are independent. If $b\ne0$ each forces its own effective time to vanish by the lower transversality. For $\delta\ne0$ this is impossible on the distinct plateaux of Calculation F4. If $b=0$, the lower row forces its own time and branch motions to vanish and vertical matching forces its flag to be the ramification point. The upper cover can be any marked degree-two cover there. This proves the zero-set assertion. When $\delta=0$ the same nodal calculation with $b\ne0$ also forces both covers to be $f$ and the flags to meet at ramification for unequal absolute heights.

## The contribution after displacement

For a generic arbitrarily small real perturbation of $D_0$ preserving $f$ and its given section, the only ordinary upper incidence on the unmarked double surface and the only augmented one on either prescribed-flag curve are at the distinguished cover $f$. To check genericity, remove $f$ from each smooth incidence chart. Any other upper cover has a different character on the complement of the branch divisors: equal characters would identify the two covers there and hence by normalization everywhere. On their common cover a nonzero ordinary or polar solution and the lower section cannot be proportional over $\mathbb R$ on a whole testing disk, by the real Leibniz rule, unique continuation, and the different characters. The projective transversality argument of [Calculation E5](../polar/polar-kernels-on-two-branched-double-covers.md) now gives projection index $-1$ on every such chart. Its étale and trivalent tests apply unchanged. At $f$ the full ordinary branch/time isomorphism and its quotient by the selected branch plane for each pole give isolated crossings even allowing time, as in [Calculation F3](the-two-endpoint-degrees-at-a-ramified-double-wall.md). Compactness and stabilization then give regular ordinary and polar upper families at every sufficiently small nonzero upper time. The finite-chart Sard–Smale conclusion uses the countable cover of the smooth locus, and the preceding local isolation controls the deleted point. The same transversality and local openness allow a smooth perturbation. The support of the tensor detecting the time crossing can be chosen after this perturbation.

Assume the projected jet condition holds at both lower ramification points. At fixed sufficiently small $\delta\ne0$ the marked surfaces for the two flags then contribute in increasing original time
$$
      -\eta_f\left(4 F_2(D_\theta)+
          \sum_{r}E(D_\theta;f(r),\mathbb C k_z(r))\right),\qquad
           \theta=\delta\big(\chi(0)-\chi(2)\big).                            \tag{4}
$$
In a neighborhood of either compact marked surface, neighborhoods of its two sheet heights over the test support are disjoint. The plateaux therefore give the fixed lower homogeneous family and the displaced upper family there, also at collisions or contracted bridges of the upper cover. The equation off the support is their common equation. At fixed $\delta$ the lower normalized derivative is $W_f$, and the upper cokernel and principal-part quotient are honest bundles by the regularity just proved. The standard marked comparison of [Calculation E6](../polar/the-quartic-contribution-at-a-marked-double-cover.md) applies: its two transverse lower rows give the original sign, its last two plane rows give the marked Euler integral, and the distinguished nonzero component orders the two roles even when the genera agree. The scale slice can be changed to the chart's one lower thick test because the equation here is homogeneous with its coefficients on the plateaux. The lower two-sheet framing and the upper marked stabilizers are consequently exactly those of Calculation E6. Uniform exclusion above treats the only point at which the upper injectivity of the limiting operator fails. On the complement of a small chart there the upper bounds in Calculation E6 hold uniformly by compactness. Thus no smoothed solutions can accumulate on the zero-height family as $\delta$ tends to zero. The number in (4) concerns the displaced, excisable marked families. Recovering the full weight at the original wall still requires the equal-height family.
