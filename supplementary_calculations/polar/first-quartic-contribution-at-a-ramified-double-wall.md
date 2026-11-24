# E2. The first quartic contribution at a ramified double-cover wall

Suppose $D_t$ is an ordinary double-cover crossing over a fixed curve $C$ of genus $g\ge2$, with base operator invertible and a smooth double cover $f:Y\to C$ simply branched over $2l>0$ distinct values. Write $h=g(Y)=2g-1+l$, let $\xi$ generate the odd kernel line, and suppose its graph $\gamma=(f,\xi)$ is embedded. Assume the Hurwitz-and-time derivative is invertible with absolute sign $\eta_f$, the nonzero projected antilinear jet condition of Article IV, §§2–3 and Appendices A–B holds at every ramification point, and all connected étale doubles $u:T\to C$ have invertible pullback $u^*D_0$. These are generic lower-wall conditions for a fixed-symbol path. For $l=1$ take also the simultaneous twist regularity of [Calculation E1](the-normal-twist-at-a-ramified-quartic-tower.md), which holds for generic paths.

Consider the connected degree-four stable maps at this crossing with a projective normal kernel whose nonzero positive images are real multiples of $\gamma$. Zero-height positive components are allowed. A “contribution” below is the sum of the local, supported Euler changes on this specified incidence. A separate simple quartic image is outside it. Let
$$
 S_C=\sum_{u:T\to C\ {\rm connected\ \acute etale\ double}}\epsilon(u^*D_0).
$$
This is a sum over isomorphism classes; its integers have absolute value one before summation. In particular $S_C$ is odd, since there are $2^{2g}-1$ such doubles. Put $A=D^N_\gamma$ for the normal operator in the homogeneous graph structure as in Calculation E1.

**Theorem.** The least source genus in this incidence is $k_0=h+2g-1=4g-2+l$. Its full supported connected Euler change at that genus is
$$
 J_f(k_0)=
 \begin{cases}
 \displaystyle \eta_f\left(\frac12\sum_{0\ne\lambda\in H^1(Y;\mathbb F_2)}
                  \epsilon(A^\lambda)-2S_C\right),&l=1,\\[5pt]
 -2l\,\eta_f S_C,&l>1.
 \end{cases}                                                        \tag{1}
$$
In particular $J_f(k_0)$ is a half-integer of odd numerator when $l=1$. When $l>1$ it is a nonzero odd multiple of $2l$. No genus below $k_0$ has any support of this type.

**Genus and classification.** For each positive normalization of a graph of nonzero height in our incidence the degree over $C$ is an even positive integer: that branch maps to the embedded graph of $f$. If one component uses all four degrees, its normalization factors through a connected degree-two cover of $Y$, hence has genus at least $2h-1$. The difference $(2h-1)-k_0=l-1$ puts it above $k_0$ except when $l=1$. Equality requires an étale upper double, a smooth source and no further components. Two separate positive components of nonzero height would both have normalization $Y$ and total genus at least $2h>k_0$.

There remains one copy of $Y$ at nonzero height and degree two at height zero. One degree-two zero-height normalization has genus at least $2g-1$ by Riemann–Hurwitz. Equality requires a connected étale double $u:T\to C$. Two degree-one normalizations would have total genus $2g$ and already exceed the budget. The arithmetic genus is the sum of the genera of the normalizations plus the first Betti number of the dual graph. Equality at $k_0$ permits no positive-genus contracted component, loop or self-node. With two positive components a tree of unmarked contracted rational components would require a contracted vertex of valence at most two: a stable tree with only two possible leaves must be a path. Thus the two components join directly at a single node. The embedded odd graph meets height zero precisely at its ramification points. Choose any such $r$ and any $s\in T$ above $f(r)$. The two choices of $s$ give isomorphic maps by the deck symmetry of $u$. Distinct choices of $r$ are inequivalent over the fixed target since their images differ. These maps exist for every $u,r$.

**The local equation for the direct join.** Work at this finite map. At its selected value, the rank-four finite flat algebra over a target disk is the direct product of the cubic nodal algebra and the algebra of the other smooth étale sheet. The idempotents lift under small deformations since the two source neighborhoods are disjoint. In translated coordinates the cubic factor has $zy=q$, $w=(z+p)^2+y$. The parameters are the $2l$ ordinary lower branch positions and increasing time, the complex motion $p$ of the attaching mark and the complex smoothing parameter $q$. The movement of $s$ is determined locally by the common target value since $u$ is étale. The finite-map tangent complex of $u$ has no kernel or cokernel. More explicitly the joined point contributes three complex parameters, and the other $2l-1$ simple branch disks one each; the remaining unramified sheets have zero relative tangent complex. This gives the smooth fixed-target stack chart of complex dimension $2l+2$, equal to $2k_0-2-4(2g-2)$.

Normalization of the vertical family gives kernel $\mathbb R\xi$ and an exact cokernel sequence
$$
     0\longrightarrow N_{f(r)}\longrightarrow P
                         \longrightarrow\operatorname{coker}f^*D_0\longrightarrow0.
$$
In fact the two normalized components have kernels $\mathbb R\xi,0$, the first evaluates to zero at the join, and $u^*D_0$ is invertible. On the cubic factor the projected matching column from the node, modulo the complex line spanned by $v_r=\partial_z\xi(r)$, has leading term a nonzero real multiple of
$$
     [B(f(r))v_r]\,\bar q\log(1/|q|).
$$
The attaching column has leading term $\pm p v_r$. These are the logarithmic estimate and residue of Article IV, §§2–3 and Appendices A–B. To apply that estimate here choose its local interpolation and residue representative on the cubic factor and extend it by zero on the fourth sheet. The singular primitive has the same value on $Y$ as in the cubic calculation and zero on *both* local branches of $T$. It consequently represents exactly the displayed matching fiber; the two sheets of $T$ reconnecting elsewhere changes no residue. On all other branch disks there are just the lower quadratic factor and two smooth unramified sheets. The stabilized global comparison uses thick test forms on each normalization; adjoint unique continuation and independent endpoint evaluation span the cokernel. The block from the lower branch and time directions is the same $W_f$; outside the cubic factor the invertible component contributes its determinant sign $\epsilon(u^*D_0)$. The remainder estimate is the one-node estimate of the cited theorem since the extra local direct factor is smooth and the global stabilized inverse stays bounded. The complex attaching column has positive degree and the antilinear smoothing quotient negative degree.

The projective kernel is isolated by this full set of transverse blocks after the radial substitution $q\mapsto q\log(1/|q|)$. The linking sphere has even dimension $4l+4$; the Euler integral there is twice the signed degree. The only automorphisms of the joined map are the lower deck group of order two: components have different genera, the lower involution fixes $r$, and the nontrivial deck transformation of $T$ moves $s$. Thus the local change per isomorphism class $(u,r)$ is $-\epsilon(u^*D_0)\eta_f$. Summing yields the term $-2l\eta_f S_C$. For $l=1$ add the smooth tower terms from Calculation E1, with their own stabilizers. Each point has an isolating neighborhood in the full fixed-genus incidence and the union is finite. Additivity by excision of the normal Euler change proves (1).

**The lower-degree factor.** The already known change in label two is the exponent $\eta_f$ of $F_h(q^2)$. Its connected contribution in label four starts at source genus $2h-1$ with coefficient $\eta_f(2^{2h}-1)/2$. When $l>1$ this starts above $k_0$. When $l=1$ subtracting it at $k_0$ leaves
$$
     -\eta_f\#\{\lambda\ne0:\epsilon(A^\lambda)=-1\}-2\eta_f S_C.
$$
This subtraction separates the known proper-divisor contribution. A formula for a whole quartic chamber would also need primitive quartic crossings and higher-genus contributions at ramified double walls.

One genus beyond $k_0$ a direct join with a zero-height double may acquire a contracted elliptic component, or the zero-height cover may acquire two simple branch values. Two zero-height identity components now fit the genus budget. When $l=1$, two copies of the nonzero double graph can also meet at ramification, with independent nonzero real heights. These configurations mark the limit of the finite-map calculation above.
