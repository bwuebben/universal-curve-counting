# F7. Excluding unequal heights near the split tower

Use the exceptional two-branch crossing and the local quartic chart of [Calculation F1](the-smoothing-at-an-interior-two-copy-height.md). Thus the graph of the odd lower generator $k$ is embedded, the lower Hurwitz-and-time derivative $W_f$ is an isomorphism, and $[B(c)v]\ne0$ at the ramification point under consideration, where $v=\partial_z k(r)$. We study the linear pullback equation with two independently assigned times $t_1,t_2$ on the two thick components. As in [Calculation F4](separating-the-two-heights-on-a-common-parameter.md), take $D_t=D_0+tK$ with $K$ supported over a disk of regular values away from the neck. Independent times therefore prescribe a single operator on the smoothed domain, with the original equation on the neck. These are the frozen equations that occur in comparing a perturbation depending on height with its values on two reference sheets.

Use positive measured amplitudes $a,b$ near the same fixed positive value, and put
$$
 h_0=a-b,\quad m=(p+s)/2,\quad n=(p-s)/2,\quad
 R=\sqrt{|q|},\quad L=|m|+R,\quad \varepsilon=\max(|p|,|s|).
$$
Here $p,s$ are the two attaching coordinates, $zy=q$, and
$$
 w=x+z^2+2pz+y^2+2sy.
$$
The thick tests for the amplitudes are exchanged by interchanging the two components. Choose the two copies of the same odd lower supplement, and use the symmetric matching supplement $(L_1\chi H/2,-L_2\chi H/2)$. The cutoff sections on the respective components have the same prescribed central value and vanish on the amplitude-test disks. As in Calculation F1 their operator images are supported on unchanged thick collars. The frame $H$ consists of local homogeneous sections on the target disk with $H(x)=1$. The matched nodal obstruction is still $M_0=bk_2(s)-ak_1(p)$, with normalized nodal sections $ak_1+\chi H M_0/2$ and $bk_2-\chi H M_0/2$. Their common marked value is the average of the two uncorrected values. Write $e_1,e_2,M$ for the two full component cokernel rows and the complex rank-two matching row of the smoothed, normalized equation. Interchange sends $(e_1,e_2,M)$ to $(e_2,e_1,-M)$ exactly, with the central cokernels identified on the two copies.

**Theorem.** Fix a constant $C_0$. Near the split join, if all these frozen rows have norm at most $C_0|h_0|L^4$, then either the two amplitudes are equal or all the geometric and time parameters vanish. In the equal-amplitude case the resulting exact solution is a double cover of the original embedded lower graph at time zero. The same assertion holds with $a,b$ in a compact neighborhood of their common positive limiting amplitude. The theorem concerns the neighborhood of this split join; it makes no claim about regularity of upper covers elsewhere.

## The symmetric tower

The involution interchanging $z,y$ has local quotient coordinate $Z=z+y+m$. With $d=z-y$ one has the exact identities
$$
            d^2=(Z-m)^2-4q,\qquad
            w=(x-m^2-2q)+Z^2+2nd.                                              \tag{1}
$$
When $n=0$ and the other two branch positions coincide, the quotient extends over the complements, which are two copies of the same exterior. Indeed on the outer annuli the identification is the unique nearby lift of the target coordinate to the old exterior covering. Interchange preserves this identification, and the quotient disk in the coordinate $Z$ is attached by the usual moving simple-branch construction. The quotient is the lower double with selected branch displacement $x-m^2-2q$ and the common other branch displacement. When these displacements and the common time are zero, pullback of $k$ gives a solution for each common amplitude, on both the smooth and nodal upper covers in (1). This supplies an exact family along which the matching and all antisymmetric component rows vanish. The matching row also vanishes for equal inputs before imposing the symmetric component equations by uniqueness of normalized stabilization.

## Parameter derivatives with the node coordinate fixed

Set $Q=|q|$. Write $x_1=x-p^2$, $x_2=x-s^2$, let $u_i$ denote the other branch displacements, and put $\lambda_i=(x_i,u_i,t_i)$, $\Lambda=\max_i|\lambda_i|$, $\lambda_-=\lambda_1-\lambda_2=(-4mn,u_1-u_2,t_1-t_2)$. Use the notation $W_r$ for the selected branch summand of $W_f$ and $E(\lambda_i)$ for the obstruction of a single normalized copy. We use the following uniform estimates as all geometric and time parameters tend to zero:
$$
 \frac{e_i}{a_i}
     =E(\lambda_i)-2\frac{a_j}{a_i}W_r(q)+F_i,qquad F_i=o(Q),\qquad
 |F_1-F_2|\le o(Q)|h_0|+C Q(|n|+|\lambda_-|).                               \tag{2}
$$
Here $a_1=a,a_2=b$, and the functions are zero on $q=0$ after the nodal term is removed. For the matching correction,
$$
 |M-M_0|\le C\big[(|h_0|+|\lambda_-|)(Q+\varepsilon R)
                                      +|n|(R+Q)\big].                         \tag{3}
$$
The factors multiplying the remainders in these comparisons are essential when $h_0$ tends to zero faster than any power of $R$.

We give the parameter estimates in detail. Keep the complex plumbing parameter $q$ fixed while differentiating. Source identifications use the same $zy=q$ cylinder, the fixed simple-branch source coordinates on the other disks, and moving identifications on collars of fixed size where the projections are unramified. These identifications and any fixed number of their derivatives in $p,s,x,u_i$ are bounded on the smaller collars, independently of $q$. They are obtained by local lifts on a slightly larger collar, as in [Calculation E6](../polar/the-quartic-contribution-at-a-marked-double-cover.md). Smooth cutoffs transport the comparison to the unchanged thick regions. The changing principal symbols in these identifications are supported on those fixed thick collars. The inner equation has the ordinary Dolbeault symbol. On the half with $|z|\ge R$, its multiplier is given by the smooth $B(w)$ and $dw/dz=2z+2p-2q^2/z^3-2sq/z^2$. Its derivatives in the indicated parameters are uniformly bounded there, since $|q|/|z|^2\le1$; the same holds on the other half. This also bounds the derivative in multiplication norm from the first weighted Sobolev space to its equation space, by the uniform section sup bound and the bounded half-collar measure. The thick tests and supplements have uniformly bounded parameter derivatives. The normalized stabilized inverse is already uniformly bounded by continuity through the nodal fiber. For $q$ fixed, differentiating the inverse uses $\partial T^{-1}=-T^{-1}(\partial T)T^{-1}$ in these identifications. Consequently both the stabilized sections and their first derivatives in these parameters and the amplitudes are bounded, uniformly in sup norm and in first Sobolev norm. Higher interior estimates hold on fixed smooth collars. No derivative with respect to $q$ at the nodal fiber is used. [Siebert, Proposition 4.15](https://arxiv.org/pdf/dg-ga/9608005v2#page=52) · [Bai–Swaminathan, Appendix B](https://arxiv.org/html/2106.01206v4#A2).

The undifferentiated $o(Q)$ in (2) is the two-sided Green projection of Calculation F1, (2). The symmetric change of matching supplement preserves both component tests: they pair to zero with the component-exact matching forces, whose primitives miss the amplitude normalizations. The two times occur on the thick pieces only. With fixed operators the stabilized solution and both rows are linear in the two prescribed amplitudes. Applying the undivided component expansion to the two basis amplitudes gives the $o(Q)$ bound also for the derivatives of $F_i$ in the amplitudes on the compact range $a,b>0$.

For each of the other parameter derivatives the bound $O(Q)$ suffices. In the component Green formula subtract the nodal component obstruction at the same parameters before differentiating. The normalized smooth-component adjoints and their parameter derivatives are uniformly smooth at the marked disk, including their inhomogeneous amplitude-normalization source off it. The latter's contribution cancels with the nodal solution. The boundary term on a fixed opposite source circle is $q$ times a bounded functional of $K$ and of the local test, also after differentiation, by the preceding inverse estimate. Here one fixes its small positive radius once for the $O(Q)$ bound. The comparison of coefficients on the upper truncation uses $G=q^2/y^2+2pq/y$, or its counterpart on the other side. The quadratic-pole derivative integrates to $O(Q)$ and a derivative in $p$ can replace $2pq/y$ by $2q/y$. For the derivative of the simple pole, freeze the coefficient and the test at the current marked value. The resulting constant angular test pairs only with the modes of order $2$ or $-2$ of the section or its parameter derivative. Both have uniformly bounded integrals in $d\log|y|$: the radial equation $(\partial_{\log|y|}-j)a_j=f_j$ has $L^1$ right side bounded by integrating a constant times $|y|+Q^2/|y|^2+Q/|y|$ over $Q/\rho<|y|<c_1$. Differentiating the equation adds bounded sections times the derivative of its multiplier, with the same integrable radial bound. For $j\ne0$ integration from the stable end of the cylinder bounds the $L^1$ norm of $a_j$ by its endpoint sup bound and $\|f_j\|_1$, divided by $|j|$. Subtracting the frozen factors, including their parameter derivatives, costs $O(Q)$: their differences from their values at the marked point are bounded by positive powers of $|y|$ or $|q/y|$ on this annulus, with uniformly bounded $d\log|y|$ integrals. The remaining coefficient differences times $dw_i$ and the collar variations cost $O(Q)$ exactly as in Calculation F1. The matrix pairing the thick representatives with the component tests may be fixed to the identity for that reference component; they are on unchanged sheets. This proves the claimed $O(Q)$ derivative bounds of the nodal-subtracted component expressions.

For (3) use the trial section of Calculation F1, $V_1(z)+V_2(y)-m_0$ on the neck, where $V_i$ are now the symmetrically matched nodal solutions and $m_0$ is their common value. It is joined to the $V_i$ in the outer collars. Its residual and its amplitude derivatives are $O(Q+\varepsilon R)$ in the weighted equation norm. This estimate remains true under a unit derivative in the component branch displacements or the times with the attaching coordinates held fixed. In the neck a branch-center translation only varies the smooth target coefficient and smooth nodal sections; their differences $V_i(\zeta)-m_0$, also after such differentiation, are $O(|\zeta|)$. At an attaching-coordinate derivative there are extra terms of sizes $O(|y|)$ and $O(|z|Q/|z|^2)$ on the $z$-half. Both are $O(R)$; they arise from differentiating $|dw_1|$ or $|dw_2|$ in the cross terms before taking absolute values. The outer cutoff and collar errors remain $O(Q)$ under these derivatives. Applying the differentiated inverse therefore gives $O(R+Q)$ for an attaching derivative of the difference from the trial, and $O(Q+\varepsilon R)$ for an amplitude, time or other branch derivative. These bounds apply to the section difference in sup norm and to the obstruction coefficients.

To subtract swapped data, keep $q,m$, the average amplitude, and the mean component parameters fixed, and interpolate linearly in the antisymmetric amplitudes, $n$, and the differences of the other branches and times. For the selected mean displacement use $x=\tfrac12(x_1+x_2)+m^2+n^2$ along the segment; the extra common derivative $2n\dot n$ is bounded. Then $\lambda_-$ is linear along this interpolation and $\varepsilon$ remains bounded by a constant times its endpoint size. At the symmetric midpoint $F_1=F_2$ and $M=M_0=0$. Integrating the derivative bounds proves (2) and (3).

## The scale forced by the lower rows

Suppose now that the output rows are $O(|h_0|L^4)$ and the data are approaching the split join. Equation (2) first gives $\Lambda=O(Q+|h_0|L^4)$, since $E(\lambda)=W_f\lambda+O(|\lambda|^2)$. On subtracting the normalized component rows the smooth Taylor error is $O(\Lambda|\lambda_-|)$. Thus
$$
 \begin{aligned}
 |\lambda_-|&\le C\{Q(|h_0|+|n|)+|h_0|L^4\},\\
 -4mn&=2(b/a-a/b)q+
           o(Q)|h_0|+O\big(Q|n|+(Q+\Lambda)|\lambda_-|+|h_0|L^4\big).
 \end{aligned}                                                                  \tag{4}
$$
Note that $\lambda_-$ includes both real time coordinates through their difference; no inference from a single frozen time is being made here.

At each branch, oddness in the branch-centered coordinate gives $k_i(p_i)=v p_i+O(|p_i|^3+|\lambda_i||p_i|)$ with the same smooth Taylor estimates after differentiating. Subtract the two expansions in $M_0$. Its leading component along $\mathbb Cv$ is $v(-h_0m-(a+b)n)$, and the remainder has size bounded by a constant times
$$
 |h_0|(\varepsilon^3+\varepsilon\Lambda)
    +|n|(\varepsilon^2+\Lambda)+\varepsilon |\lambda_-| .
$$
Using (3), (4) and the assumed row tolerance, one obtains
$$
                 |n|\le C|h_0|(|m|+Q).                                         \tag{5}
$$
Indeed the terms on the right before absorption are $|h_0|[|m|+Q+\varepsilon R+\varepsilon^3+\varepsilon\Lambda+L^4]$, together with $o(|n|)$ and $O((\varepsilon+Q+\varepsilon R)|\lambda_-|)$. Insert (4), $\varepsilon\le |m|+|n|$, and $L^4\le C(|m|^4+Q^2)$ and absorb the small multiples of $|n|$.

If $h_0=0$ then (4)–(5) give $n=0$ and $\lambda_-=0$. The rows are exactly zero. The unique normalized section is invariant under interchange and descends to the normalized lower cover in (1): evaluate on regular sheets and use boundedness and removability at its branch values and nodal limits. The lower crossing forces its selected displacement $x-m^2-2q$, other branch displacement and time to vanish. This proves the tower assertion.

For $h_0\ne0$, divide the second equation of (4) by $h_0$. Since $(a/b-b/a)/h_0=(a+b)/(ab)$ stays bounded below, (5) gives $Q\le C(|m|^2+L^4)+o(Q)$, hence $Q=O(|m|^2)$. If $m=0$, then $q=n=0$, and the zero tolerance and the component equations imply all displacements and times vanish. Otherwise (4)–(5) and the matching leading component refine to
$$
 \begin{aligned}
 n&=-\frac{h_0}{a+b}m+O(|h_0||m|^2),&
 |\lambda_-|&=O(|h_0||m|^2),&\Lambda&=O(|m|^2),\\
 \frac{q}{m^2}&=-\frac{2ab}{(a+b)^2}+o(1).
 \end{aligned}                                                                  \tag{6}
$$
In the last step $L^4=O(|m|^4)$ and the coefficient of $q$ after division remains nonzero. Consequently $R\asymp |m|$ and $\varepsilon=O(|m|)$ uniformly even when $|h_0|$ is exponentially small. These conclusions require only the indicated component rows and the first matching plane.

## The logarithm after subtraction

On the cone (6) put $\mu=|m|$. The matching expansion with the factor of height difference is
$$
  M=M_0-4B(c)v\,\overline{q(bp-as)}\log(1/Q)+O(|h_0|\mu^3).                   \tag{7}
$$
We verify the error estimate, rather than subtract two unconditional $O(\mu^3)$ error bounds. Join the symmetric inputs to the given ones by the segment used above, parametrized now by its signed height difference. Along it $q,m$ are fixed, $R\asymp\mu$ and all data obey
$$
 |p_i|=O(\mu),\quad |\lambda_i|=O(\mu^2),\qquad
 |\dot a_i|=O(1),\quad |\dot p_i|=O(\mu),\quad |\dot\lambda_i|+|\dot x|=O(\mu^2).
                                                                                \tag{8}
$$
A dot denotes differentiation along this segment. The component times need not agree along it. It is unnecessary to solve any of the obstruction equations at its intermediate points.

The exact angular matching identity of Calculation F1, (8), holds with the symmetric supplement. Use opposite normalized point evaluations at the two flags; the two thick matching forces contribute half the same evaluation each after taking their difference. This gives precisely the coefficient $M$. Component supplements drop out and the reference normalizations give $M_0$. In its notation let $H_i(\zeta)d\zeta$ be the evaluation adjoint on one component with leading coefficient $\ell/\zeta$ for a fixed complex covector $\ell$, and set $H_{B,\ell}(d)=\overline{\ell(B(x)d)}$. In the fixed complex frame near the mark,
$$
 H_i(\zeta)=\ell/\zeta+2(\zeta+p_i)H_{B,\ell}\log|\zeta|^2+U_i(\zeta).
$$
On $R\le r'=|\zeta|\le\rho$ both $U_i$ and $\dot U_i$ are bounded by $C(r'+\mu)$. To see the differentiated bound, subtract the displayed pole and logarithmic primitive at the current parameters before taking the derivative. The pole has zero derivative; the derivative of the coefficient of the logarithm is $O(\mu+r'\mu^2)$. The ordinary adjoint equation for the remainder has multiplier bounded by $C(r'+\mu)$, whose derivative under (8) is $O(\mu+r'\mu^2)$. Its local inhomogeneous term and the differentiated one are bounded by $C(r'+\mu)^2(1+|\log r'|)$, since $|w_i-x|\le C(r'^2+\mu r')$ and $|\dot w_i-\dot x|\le C\mu r'$. Derivative forcing from the multiplier on $U_i$ has an additional bounded factor. On the circle $r'\ge R\asymp\mu$ the Cauchy difference estimate in Calculation F1, (7), applies also to this differentiated equation. At the mark the residual's regular value vanishes at $p_i=0$ for every normalized component by deck parity. Its parameter dependence is smooth after the explicit pole and log subtraction in a coordinate centered at the flag, by the smooth normalized component inverse; the subtracted forcing is smooth in every finite $L^k$. It follows from (8) that both $|U_i(0)|$ and $|\dot U_i(0)|$ are $O(\mu)$. This gives the claimed bound. One may first use the uniform sup estimate for $U_i,\dot U_i$ from the component inverse to bound the right side in the Cauchy estimate, then use their values at zero. The time change is supported off the entire disk.

The matched trial and its parameter derivative both differ from the respective actual sections by $O(\mu^2)$ in sup norm: use the residual derivative bounds proving (3), with the attaching derivatives multiplied by $O(\mu)$. In particular on the middle circle the identity
$$
              K=v(az+by)+m_0+O(\mu^2)
$$
has the same error order after differentiating, with $|m_0|+|\dot m_0|=O(\mu)$. Here $K$ now denotes the normalized section in the matching calculation. On the first half-annulus, $|K|+|\dot K|\le C r'$ and
$$
 |K-avz|+|\dot K-\dot a vz|\le
                   C\big(\mu+Q/r'+r'^3+\mu^2\big);
$$
the analogous estimate holds on the other half. These are the expansions of the two odd nodal sections about the branch centers, with their matched constants. Their moving first derivatives in the other-branch, selected-branch and time coordinates carry the $O(\mu^2)$ weights in (8).

For specificity, write $d_0=q(bp-as)$. In the difference of the two middle-circle integrals the pole terms cancel for both $K$ and $\dot K$, since the annulus and its middle circle are fixed in this differentiation. The logarithmic primitive yields $-4\operatorname{Re}(\overline{\ell(B(c)v)}d_0)\log(1/R)$ as in Calculation F1. Its differentiated error is $O(\mu^3)$: the regular adjoint and its derivative have size $O(\mu)$, the section and its derivative have size $O(\mu)$, and the circle factor has size $R$; the $O(\mu^2)$ error in the section times its logarithmic coefficient costs $O(\mu^4|\log\mu|)$. The common constant and its derivative have no zeroth angular mode in the logarithmic product. Freezing $B(x)$ at $c$ and differentiating this replacement incurs the same bound by $|x|+|\dot x|=O(\mu^2)$.

For the area terms on the first half the displacement is $G=q^2/z^2+2sq/z$; thus $\dot G=2\dot s q/z$. The term with $-2sq/z^2$ in $G'$ against the leading $\ell/z$ and $avz$ is integrated exactly before differentiation. On the opposite half interchange the indices. Their signed sum contributes another $-4\operatorname{Re}(\overline{\ell(B(c)v)}d_0)\log(1/R)$. The remainder from the quadratic pole and its derivative is $O(Q^2/R)=O(\mu^3)$ by the half-annulus bounds. For the simple pole, replacing the section or its derivative by its leading term costs at most
$$
  C Q\mu\int_R^\rho
            (\mu+Q/r'+r'^3+\mu^2)\frac{dr'}{r'^2}=O(\mu^3),
$$
also when $s$ is replaced by $\dot s$. The adjoint after subtraction of its pole and the derivative of this remainder are bounded by $C(r'+\mu)(1+|\log r'|)$, giving an integrable bound after extracting $Q\mu$. Variations of $B(w)-B(w_i)$ have the original $O(Q^2/r'^2+Q\mu/r')$ bound, since the derivative of the displaced part is $O(Q\mu/r')$ and the derivative of the smooth coefficient is bounded; variations of $w_i'$ have size $O(\mu)$. The thick-collar comparison, including its derivative under (8), costs $O(Q^2+Q\mu)$, by smoothness there. These estimates give $O(\mu^3)$ for the derivative of the entire error in the signed area terms. The domains of the radial integrals, in particular both limits of every logarithm, are fixed along this segment.

Thus after taking out the summed logarithm the derivative of $M-M_0$ has an error of size $O(\mu^3)$. At the symmetric midpoint the actual difference, its nodal counterpart and the explicit logarithmic term all vanish. Integration over the parameter segment of length $|h_0|$ proves (7). The use of moving-component evaluation tests is important here; continuity of an unframed stabilization alone would not justify division by the height difference.

Under (6) the projection of $M_0$ to $N_c/\mathbb Cv$ is $O(|h_0|\mu^3)$ by the nodal Taylor estimate preceding (5). On the other hand
$$
           bp-as=-h_0m+(a+b)n=-2h_0m+O(|h_0|\mu^2).
$$
Since $Q\asymp\mu^2$ and $[B(c)v]\ne0$, (7) gives a lower bound by a positive multiple of $|h_0|\mu^3\log(1/\mu)$ for this matching quotient, for all sufficiently small nonzero $\mu,h_0$. This contradicts the prescribed $O(|h_0|L^4)$ tolerance. The preceding cases $h_0=0$ and $m=0$ complete the proof.

The conclusion identifies exact equal-height zeros and rules out approach by unequal frozen heights with nonzero geometry. Applying it to a nonlinear normalized displacement requires a comparison with the frozen rows at the tolerance used in the theorem. The normalization constraint, target perturbation, and the stabilizers in a compact Euler-class comparison are additional data.
