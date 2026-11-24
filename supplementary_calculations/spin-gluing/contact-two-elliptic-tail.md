# D2. A contracted elliptic tail at a contact of order two

In the étale smoothing of [Calculation D1](even-contact-spin-matching.md), move the attaching point of a contracted elliptic tail into the nonseparating neck. Write $g\ge2$ for the genus of the smooth target. Its normalization $B$ at the node has genus $g-1$. The connected double cover $Y\to B$ branches just at $P,Q$ over the two ends of the neck. Fix the target, its spin data, the identification on the complement of the neck and a choice of smoothing root as in Calculation D1.

There are three conclusions. A separating contracted component contributes a hyperbolic Hodge summand to the spin complex, also in nodal families. On the balanced expansion of this particular contact the pullback of the exceptional target spin has one complex bridge kernel per spin. The corresponding defect on the untailed family has length one in each of two cohomological degrees. Its first obstruction to smoothing singles out one of the two null lines of the residue plane. Changing the target spin exchanges them. Finally, the signed Euler integral of the smooth elliptic-tail component specializes on the compact parameter space of attaching points. In the exceptional target-spin extension its bridge part has magnitude $1/12$ per smoothing root. That part depends on the extension. Below we distinguish this specialization of an explicit bundle from a relative virtual contribution.

## 1. Moving the attachment through the node

The smooth cover has plumbing $uv=s$, $xy=t=s^2$, $(x,y)=(u^2,v^2)$. Write $(a,b)$ for the point of attachment, so $ab=s$. Near $a=b=0$ introduce a source rational bridge with coordinate $Z$ and a target rational bridge with coordinate $W$. The map and end plumbings are
$$
 W=Z^2,\qquad u=aZ,\quad v=b/Z,\qquad x=a^2 W,\quad y=b^2/W.
$$
These equations specify charts at both ends, for example $uZ^{-1}=a$, $xW^{-1}=a^2$ on the first. The attaching point is $Z=W=1$. At $a=b=0$ the normalization $Y$ joins the bridge at its two ends and the pointed tail $F$ attaches at $1$. Contracting the target bridge collapses the entire source bridge to the node. Relative expansions and their automorphisms provide the ambient local model for these equations. [Jun Li, expanded degenerations](https://arxiv.org/abs/math/0009097)

Over the fixed target with labeled branches, a bridge scaling preserving the attaching point is the identity, which also fixes the target scaling by $W=Z^2$. The involution on $Y$ fixes $P,Q$ and acts on $(a,b)$ by $(-a,-b)$ while fixing the marked bridge in the central expanded chart. The pointed tail has its ordinary stack stabilizer. At smooth parameters the deck group acts freely on the moving point of an étale cover. The quotient by the deck group is taken once on the whole parameter space. The chart is over the specified smoothing root $s$. The equation $t=(ab)^2$ accounts separately for the two roots.

The first-order obstruction to smoothing the tail attachment, in the tangent direction of the expanded target, is the map
$$
 T_1\mathbb P^1_Z\otimes T_eF\longrightarrow
 T_1\mathbb P^1_W\otimes H^1(F,\mathcal O_F),\qquad
 d(Z^2)_1\otimes\partial_e.
$$
Here $e$ is the attaching point on $F$ and $\partial_e$ is the principal-part connecting map for $\mathcal O_F(e)$. In the chosen coordinates the derivative is two. For a stable one-pointed curve of arithmetic genus one, $\partial_e$ is an isomorphism, including on the irreducible nodal tail. This follows from $H^1(\mathcal O_F(e))=0$ and $H^0(\mathcal O_F(e))=\mathbb C$. Equivalently, the coefficient linear in the plumbing parameter of the tail is the principal part with this derivative. This identifies the attachment smoothing with its tangent obstruction on the expanded chart. On contraction the bridge map is constant, so the same tangent calculation cannot be used there. The other two nodes of the expanded cover have contact order two.

## 2. Separating a contracted component in the spin complex

Here is a family statement independent of the covering degree. Let $Z=S\cup F$ over a base $T$ be obtained by gluing proper flat families of connected nodal curves along smooth sections, with a single separating node between $S$ and $F$. Suppose $F$ has constant arithmetic genus $k$. Write $\pi$ for the curve projections, $L$ for the line on $S$, $\ell$ for its value at the section and $L_Z$ for its extension by $\pi_F^*\ell$ on $F$ with the given attaching identification. Write $\mathbb E_F=\pi_{F*}\omega_{F/T}$. Given a map $\sigma_S:L\to\omega_{S/T}L^{-1}$, extend it by restriction to the core followed by extension of regular core differentials by zero. This defines $\sigma_Z:L_Z\to\omega_{Z/T}L_Z^{-1}$. Let $K_Z=[L_Z\xrightarrow{\sigma_Z}\omega_{Z/T}L_Z^{-1}]$ in degrees 0 and 1 and define $K_S$ similarly. Then
$$
 R\pi_{Z*}K_Z\simeq R\pi_{S*}K_S\ \oplus\
 \big((\ell\otimes\mathbb E_F^\vee)\oplus(\ell^{-1}\otimes\mathbb E_F)\big)[-1].
 \tag{1}
$$
The last two summands are dual under the symmetric middle-cohomology pairing and each is isotropic. This decomposition uses the specified separating component and section and commutes with base change.

Indeed the normalization sequence expresses $R\pi_{Z*}L_Z$ as the homotopy fiber of the difference of the two evaluations onto $\ell$. Evaluation from the contracted component has a canonical right inverse given by constants. Its homotopy fiber is $R\pi_{F*}(\pi_F^*\ell(-e))$, canonically $(\ell\otimes\mathbb E_F^\vee)[-1]$ by evaluation of constants and duality. Thus the first term splits as $R\pi_{S*}L$ and this summand, and restriction to $S$ is the first projection. Relative Serre duality gives the dual splitting, shifted by $[-1]$, for the second term of $K_Z$. Extension of regular core differentials by zero is dual to restriction of sections to $S$, as follows from the residue description of the dualizing sheaf at a separating node. Consequently the pushforward differential is supported in the core block. This proves (1) and the hyperbolic pairing on the extra two summands. The calculation uses only normalization, constants and relative duality on proper nodal curves, including nodal $F$.

Over the smooth-target elliptic-tail component the actual normal cokernel for one spin is the first complex line $\ell\otimes\mathbb E_F^\vee$, viewed as real, with the orientation factor of the core. This follows by cancellation of the constant tail sections with the node-matching fiber in the real operator normalization sequence. The second term in (1) supplies its dual in the algebraic complex. For two spins the smooth tail sign is the product of the two core signs. Calculation D1 gives $\sigma=-1$ on each of the two roots for the mixed pair and $\sigma=+1$ for an equal pair. Article V, §§4–6

## 3. The two spin extensions and the bridge kernel

First use the balanced base change $s=d^2$ and resolve $uv=d^2$. The resulting regular semistable core surface $p:\mathcal S\to\Delta_d$ has special fiber $Y\cup T_0$, where the rational bridge $T_0$ meets $Y$ twice and $T_0^2=-2$. On its interior $u=dZ$, $v=d/Z$. On the resolved spin target $xy=s^2$ use the exceptional coordinate $W=x/s$ so $W=Z^2$. The dualizing sheaf pulls back under stabilization of this two-pointed rational bridge: locally its frame is $du/u=dZ/Z$ in the relative logarithmic coordinates. For the étale smoothing the source spin extension $H$ from Calculation D1 therefore satisfies $H^2=\omega_{\mathcal S/\Delta}$, with $H|_{T_0}$ trivial.

For a locally free target spin both its target and source extensions pull back to $H$. For an exceptional target spin let $L$ be the pullback of the Cornalba spin line of degree one on the target exceptional component. In the complement identification there is an inclusion
$$
                        L=H(-T_0)\ \subset H.
 \tag{2}
$$
For example in the first branch frame $\ell_x=u n_x$ from Calculation D1; on the opposite side $\ell_y=v n_y$. At an interior point of $T_0$, $u=dZ$ and $v=d/Z$, so the inclusion has order one along $T_0$. At its left endpoint use the chart $u\xi=d$, where $T_0$ has equation $u=0$ and the frame ratio is $u$; the right endpoint uses $v$ in the same manner. There are no other zeros. On this surface the core differential is the inclusion $H(-T_0)\to H(T_0)$, with the unit from $H^2\simeq\omega$ understood. It agrees with the spin square map and the cover differential on the dense smooth fibers and extends over the surface. The spin data on the target exceptional component are those of [Cornalba](https://mate.unipv.it/cornalba/papers/theta_moduli.pdf).

The family complex for $L$ is quasi-isomorphic to $Q[-1]$ for $Q=H(T_0)/H(-T_0)$. It has the filtration
$$
 0\longrightarrow H|_{T_0}\longrightarrow Q
   \longrightarrow H(T_0)|_{T_0}\longrightarrow0.
$$
The first term is a trivial line on the bridge and the last has degree $-2$. It follows that $p_*Q$ and $R^1p_*Q$ each have length one at $d=0$, and all other direct images vanish. In particular $d$ kills each. Locally over the disk the pushed-forward core complex is noncanonically the sum of two complexes $[\mathcal O\xrightarrow{d}\mathcal O]$ in degree pairs $(0,1)$ and $(1,2)$. Indeed there are no other cohomology sheaves and a discrete valuation ring has homological dimension one. Derived base change gives fiber dimensions $(1,2,1)$ at $d=0$ and zero elsewhere. For $H$ the core complex is acyclic on the whole disk. Equation (1) gives the corresponding correction for a family with a persistent contracted tail. At genus one it adds two dimensions in degree one per spin. The two torsion lengths refer to the parameter $d$ of this semistable family.

The singular-fiber kernel can also be seen directly. The target-pulled $L$ restricts to $f_0^*A_0$ on $Y$ and to $\mathcal O(2)$ on $T_0$. The square map of the exceptional spin is zero on the target bridge. On the source bridge the ordinary component operator is $\bar\partial_{\mathcal O(2)}$. On $Y$ it is the smooth spin operator for the branched double cover $Y\to B$ with ordinary differential. It is injective and its real cokernel has rank two by the spin-residue calculation in Article V, §§4–6. Evaluation of $H^0(\mathcal O(2))$ at both ends is onto; its kernel $H^0(\mathcal O(2)(-0-\infty))$ is a complex line. Normalization thus gives this line as the matched core kernel and the real plane on $Y$ as its cokernel. At the third point on the bridge a constant tail section takes precisely the required value. Adding the tail preserves the kernel and adds $H^1(F,\mathcal O_F)\otimes L|_1$ to the cokernel, as a quotient or direct summand after a splitting. On each smooth étale core the kernel vanished. This is an explicit failure of constant rank for the componentwise target-pulled operator; transport of an ordinary cokernel bundle of that operator across the seam is impossible. The derived calculation specifies the jump without making an identification of real small-eigenvalue spaces during the smoothing.

### The obstruction to extending the bridge section to first order

Write $A^\bullet$ for the surface core complex and $A_0^\bullet$ for its ordinary termwise restriction to $Y\cup T_0$. The first Bockstein
$$
 \beta:\mathbb H^0(A_0^\bullet)\longrightarrow\mathbb H^1(A_0^\bullet)
$$
is the connecting map for the extension modulo $d^2$ by the special fiber, with the factor $d$ divided out. We identify its image with a line in the ordinary orthogonal residue plane of the ramified normalization. The coefficients below use the branch coordinates of Calculation D1.

The sheaf kernel of the central differential consists of bridge sections of $L|_{T_0}$ vanishing at both ends. It has degree zero on the bridge and $H^1=0$. The sheaf cokernel maps onto $M|_{T_0}$ of degree $-2$, where $M=\omega L^{-1}$. At the left node its kernel has coefficients $u$ modulo $u^2$ in the frame $(du/u)\ell_x^{-1}$: the first factor imposes zero value on the bridge and the second is the image of the square map. The right node gives the corresponding coefficient of $v$. Together these identify the kernel with the ordinary ramification quotient on $Y$, since multiplying $du/u$ by $u$ gives the regular branch differential and $df=2u\,du$ there. The exact sequences give
$$
 \mathbb H^1(A_0^\bullet)\simeq H^0(\operatorname{coker}\sigma_0)
     \simeq E_Y.
$$
Here $E_Y$ is the sum of the two residue lines at $P,Q$. Use polar coefficients $c_P/u$ and $c_Q/v$ in the pulled-back branch frames. The residue form has quadratic expression $2(c_P^2+c_Q^2)$.

Choose the branch $\delta=s$ for the target spin parameter. On the bridge the transition is $\ell_y=\varepsilon i Z^{-2}\ell_x$. Its endpoint-vanishing section is
$$
       Z^{-1}\ell_x=(\varepsilon i)^{-1}Z\ell_y.
$$
Near the left node lift the coefficient to $\xi$ in the chart $u\xi=d$. The differential of that lift is $2u^2\xi(du/u)\ell_x^{-1}$. Dividing by $d$ gives $2\,du\,\ell_x^{-1}$, with polar coefficient $1/u$. Near the right node the coefficient is $(\varepsilon i)^{-1}/v$ by the same calculation. These local images determine the Bockstein: differences between lifts give images of the central differential in the cokernel, and $H^1(\ker\sigma_0)=0$ removes any remaining ambiguity. Consequently
$$
       \operatorname{im}\beta_\varepsilon
          =\mathbb C\big(1,(\varepsilon i)^{-1}\big)\subset E_Y.
$$
The two spin choices give the two distinct null lines. Reversing the smoothing-root labeling exchanges them. There is a real orientation test on this finite-dimensional model. For the standard positive real form $\mathbb R^2\subset E_Y$, project to it by taking real parts. On $c=x+iy$ the two projected maps $c\mapsto\operatorname{Re}\big(c,c/(\varepsilon i)\big)$ have matrices $\operatorname{diag}(1,1)$ and $\operatorname{diag}(1,-1)$. Every positive real form gives the same orientation ratio, since the space of such forms is connected and the real projection from either null line stays invertible (its kernel would lie in the imaginary multiple of a positive real form). Complex rescalings of the bridge generator preserve its real orientation. The resulting ratio $-1$ agrees with the independently computed smooth sign. Identifying this projected map with an actual real Fredholm reduction through the degeneration requires a comparison beyond the algebraic first variation.

## 4. The compact integral over the moving attachment

For one prescribed smooth cover and one smoothing root the parameter stack is $[(X_s\times\overline{\mathcal M}_{1,1})/\mu_2]$, with $X_s$ the compact source and $\mu_2$ its deck group. For a pair of evaluation lines $L_1,L_2$ on the attachment curve, the sign-transported Hodge-bundle Euler class is the sign $\sigma$ times the second Chern class of $\mathbb E^\vee\otimes(L_1\oplus L_2)$. On the product of the two curves only $-\lambda(c_1(L_1)+c_1(L_2))$ contributes. Using $\int_{\overline{\mathcal M}_{1,1}}\lambda=1/24$ (including the tail stabilizers), the value over this quotient is
$$
       -\frac{\sigma}{48}(\deg L_1+\deg L_2)
       =-\frac{\sigma(g-1)}{12}.
 \tag{3}
$$
Each line has degree $2g-2$. Thus the absolute contribution per root is $+(g-1)/12$ for the mixed pair. This agrees with Article V, §§4–6. Specialization below uses the smooth orientation sign as a coefficient on the test cycle.

There is an equivalent calculation on the original surface $uv=s$ of moving points. For the exceptional target model the direction $[x:s]=[u:v]$ away from the origin is indeterminate there. Blow up the origin to obtain the proper modified family of attachment curves, with exceptional curve $D$ of self-intersection $-1$. The special fiber cycle is $Y+2D$, since $s=uv$. Here $D$ maps isomorphically to the target exceptional line by $[u:v]$; this is a map on attachment parameters. On the balanced family of §3 the map between the two rational components of an individual cover has degree two. On the blowup the geometric target evaluation line is $\rho^*H_X(-D)$, where $H_X$ is the source extension on $uv=s$, as follows from the two frame ratios $u,v$. For the local deck involution both the exceptional source frame $n_x$ and $u$ change sign, so their product is the invariant target frame. Both extensions descend with their natural linearizations to the quotient stack.

For either exceptional target spin, $\rho^*H_X$ has degrees $(2g-2,0)$ on $(Y,D)$. Its twist by $-D$ has degrees $(2g-4,1)$, because $Y$ meets $D$ twice. The reduced balanced model of §3 gives degrees $(2g-4,2)$ on $(Y,T_0)$ for the target pullback. The core and bridge parts of (3), on the attachment-space quotient and including special-fiber multiplicities, are therefore

| Choice for both evaluation lines | Main normalization | Exceptional curve | Total per root |
|---|---:|---:|---:|
| Source spin extended by the dualizing square root | $-\sigma(g-1)/12$ | $0$ | $-\sigma(g-1)/12$ |
| Pullbacks of exceptional target spins | $-\sigma(g-2)/12$ | $-\sigma/12$ | $-\sigma(g-1)/12$ |

For the two distinct even gluings on a genus-two target the full $+1/12$ occurs on the exceptional curve for the target-pulled evaluation. It occurs on the main normalization for the source extension. In the locally free target model both choices already agree and the bridge degree is zero. The divisor calculation takes place on complete attachment curves times the compact pointed-tail stack. Quotient degrees use the inherited deck action and (3) gives the smooth-component Euler integral. These degrees do not determine an additional automorphism factor in a relative stable-map virtual class.

## Consequences for gluing

The attachment equation, hyperbolic splitting and rank jump distinguish the tangential and normal parts of the problem. The derivative at the tail mark is invertible in the smoothing obstruction even while the target-pulled normal operator has a kernel on the bridge. Twisting by $+T_0$ restores the dualizing spin complex on the untailed core; (1) tracks the accompanying change of the evaluation line on every separating tail. Equation (3) fixes the compact total while the two extensions assign different integrals to the bridge. A relative degeneration theorem must account for the real orientation of the disappearing kernel and cokernel and identify the resulting relative obstruction class with the absolute one under gluing. The algebraic length-one complexes and the attachment-space integral give explicit data for that comparison. Its real determinant gluing and the relative virtual pairing remain open. Lee and Parker's localization of small modes concerns spin Hurwitz parity for different smoothings; their analysis supplies a precedent for investigating this comparison. [Lee–Parker, §§7–10](https://arxiv.org/abs/1212.1825)
