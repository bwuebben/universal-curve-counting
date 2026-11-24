# D1. Spin matching at one contact of order two

We compute the orientation sign on a finite contact chart for the nonseparating cut of Article V, §2. The orientation character belongs to the orthogonal normal bundle of the target curve. In a split representative the two spin lines are $\theta$ and $\theta\eta$, where $\eta$ is the flat real sign line complexified. We choose the cut dual to $w_1(\eta_{\mathbb R})$ and a parallel unit of $\eta$ on its complement. Changing the identification across this cut by a minus sign changes the spin gluing. In the two isotropic normal coordinates this interchanges the lines. That exchange on a pair of complex lines preserves their joint real orientation. The sign computed below concerns the oriented real cokernel on the source.

## Proposition

Let $B$ be a smooth connected compact complex curve of genus at least one with distinct points $p,q$, and let $f_0:Y\to B$ be a smooth connected double cover with a single preimage $P,Q$ of each respective point. Form the nodal target by identifying $p$ with $q$ and the nodal source by identifying $P$ with $Q$. Plumb the cover in coordinates
$$
xy=t,\qquad uv=s,\qquad x=u^2,\qquad y=v^2,\qquad t=s^2.
$$
Keep the complement of the plumbing discs fixed. Let $R$ be the ramification divisor away from $P,Q$, of even degree $2r$. Fix either root $s$ of a small nonzero $t$. Compare the two target spin structures $\theta_+,\theta_-$ agreeing off the cut and differing in the gluing sign. For either target spin type described below, the canonical orientations of the two real spin cokernels on this **same smoothed cover** are opposite under the common orthogonal identification of their ramification residue spaces. At $r=0$ this compares the signs of the index-zero determinants.

In particular the ordered mixed pair $(\theta_+,\theta_-)$ differs in orientation from the equal pair $(\theta_+,\theta_+)$ by $-1$. If $e$ is the sign of this mixed pair on the identity cover, the degree-two normalization divides by $e^2=1$. The comparison holds on charts with the indicated complement identification and on either smoothing branch. The following proof extends the *algebraic orientation test* across the finite nodal cover.

## 1. The two spin plumbings

Write $D=P+Q$. The two types over an irreducible one-node target occur in Cornalba's description of stable spin curves. [Cornalba, §6, Example (6.2)](https://mate.unipv.it/cornalba/papers/theta_moduli.pdf)

* **Locally free target spin (Ramond).** Choose $A$ on $B$ with $A^2=K_B(p+q)$. Take frames squaring to $dx/x$ and $dy/y$. Since $dy/y=-dx/x$ on the neck, either choice $e_y=\varepsilon i e_x$, $\varepsilon=\pm1$, sews a spin line. On the source normalization put $H=f_0^*A$. The logarithmic differential sends the squares of the two pulled-back frames to $2du/u$ and $2dv/v$. The source gluing has the same two choices for either $s$. The two target gluings have opposite spin parity.
* **Exceptional target spin (Neveu–Schwarz).** Choose $A_0^2=K_B$ and frames $\ell_x^2=dx$, $\ell_y^2=dy$. Choose $\delta$ with $\delta^2=t$. The two smooth target identifications are $\ell_y=\varepsilon i\delta x^{-1}\ell_x$. On the smoothed source take $n_x=f^*\ell_x/u$ and $n_y=f^*\ell_y/v$. Their transition is
$$
                         n_y=\varepsilon i(\delta/s)n_x.
$$
Thus on either branch $s=\pm\delta$ the two source identifications differ by sign. The frames extend as generators over the plumbing, matching $H=f_0^*A_0(D)$ on the normalization. Their squares under the differential are again $2du/u$ and $2dv/v$. More explicitly, use $n_x$ as a generator on the source plumbing surface and the displayed constant transition for $n_y$ on the second overlap. The relative dualizing sheaf there has generator $du/u=-dv/v$. The resulting square map is a unit at the node and agrees with the usual differential off the neck. This constructs the extended source line directly. Both smooth target spins specialize on the exceptional model to $A_0$ on $B$ and $\mathcal O(1)$ on the exceptional component. Evaluation of $\mathcal O(1)$ at the two attaching points is an isomorphism, so their common parity is that of $A_0$. In particular an even $A_0$ gives a pair of even target spins.

In each case $H^2=K_Y(D)(-R)$ and the resulting source line $H_\varepsilon$ on the nodal curve $X=Y/(P\sim Q)$ is locally free. Write $M_\varepsilon=\omega_X H_\varepsilon^{-1}$. The differential defines $\sigma:H_\varepsilon\to M_\varepsilon$, invertible at the node with divisor $R$. On the smooth fibers the pair is the usual pulled-back spin and differential. At fixed choice of $s/\delta$ the two matching conditions at the node differ by a sign. Interchanging the two roots interchanges their labels in the exceptional model.

## 2. The matching calculation on the compact normalization

On $Y$ put $M=K_Y(D)H^{-1}$, $E=H^0(M/\sigma H)$ and $F=M_P\oplus M_Q$. Here $E$ is the same space for both gluings. Its symmetric form is the sum over $R$ of the residues of $vv'/\sigma$. The same expression at $P,Q$ pairs the fiber values in $F$, since there is at most a simple pole on each branch. Each branch gives a nonzero diagonal coefficient. Hence $F$ is a nondegenerate complex quadratic plane. The two lines $T_+,T_-\subset F$ of matching values for $M_\varepsilon$ are distinct. They are null because the logarithmic residues of a dualizing differential on $X$ sum to zero at the node. They are therefore the two null lines of $F$.

Restriction and evaluation give a map
$$
 H^0(Y,M)\longrightarrow E\oplus F.
$$
Its kernel is $\sigma H^0(Y,H(-D))$. The image $W$ is isotropic by the global residue theorem. If $Y$ has genus $h$ then $\deg H=h-r$, $\chi(M)=r+1$ and $h^1(M)=h^0(H(-D))$. Thus $\dim W=r+1$, half the dimension of $E\oplus F$.

Choose any maximal isotropic $U\subset E$. For each central gluing define $d_\varepsilon$ as the dimension of the space of global sections of $M_\varepsilon$ with residue class in $U$. Then
$$
 d_\varepsilon=h^0(Y,H(-D))+\dim\big(W\cap(U\oplus T_\varepsilon)\big).
$$
Let $K=U\oplus0$. Since $W=W^\perp$, $\dim(W\cap K^\perp)=\dim(W\cap K)+1$. Indeed the map $W\to K^*$ has rank $\dim K-\dim(K\cap W)$. The image of $W\cap K^\perp$ in $K^\perp/K=F$ is exactly one of its two null lines. One of the two intersections in the formula for $d_\varepsilon$ therefore has one more dimension than the other. In particular $d_+-d_-$ is odd, including when $E=0$.

## 3. Transport to smooth fibers and the determinant sign

The complement identification identifies the two orthogonal residue spaces throughout the plumbings. Extend $U$ there (it may be held fixed for the chosen family). To compare $d$ through jumps of cohomology, choose a divisor $G$ of sufficiently many smooth, distinct auxiliary sections away from the neck and $R$. Shrinking the plumbing disk if necessary, arrange $H^0(H_s(-G))=H^1(M_s(G))=0$ on the central and nearby fibers, with $M_s=\omega_s H_s^{-1}$. Such a choice exists by duality and sufficiently large degree on the irreducible central curve. For either spin gluing, the space
$$
 J_s=H^0\big(M_s(G)/\sigma H_s(-G)\big)
$$
is a nondegenerate orthogonal bundle by the same local residue pairing. It splits into $E$ and the length-two jets at each point of $G$. The image of $H^0(M_s(G))$ in $J_s$ is an isotropic subbundle by the residue theorem and base change. At the central node the two residues cancel. Its rank is $r+\deg G$ by Riemann–Roch. The regular jets at $G$ form an isotropic subbundle of rank $\deg G$. Together with $U$ they give a second maximal isotropic. Their intersection dimension is $d_s$ since regular jets enforce absence of poles at $G$. In a continuous pair of maximal isotropics the intersection dimension modulo two is locally constant. In a hyperbolic basis replacing an isotropic coordinate by its dual reverses the orthogonal determinant, giving the component criterion $(-1)^{n-\dim(A\cap B)}$ in complex rank $2n$. This proves that the central parity of $d_\varepsilon$ persists on its smoothing.

On a smooth fiber set $S=\operatorname{im}(H^0(M_s)\to E)$ and $k=h^0(H_s)$. The residue theorem for the spin operator in Article V, §§4–6 identifies its real cokernel with a positive real form $V\subset E$. With the orientation continued from the complex Dolbeault operator, the real projection $V\to E/S$ has sign $(-1)^k$. Passing from $S$ to $U$ conjugates $r-\dim(S\cap U)$ complex coordinates in the projection. Since $d_s=k+\dim(S\cap U)$, the sign of $V\to E/U$ is
$$
                         (-1)^{r+d_s}.
$$
Both spins are compared with the same complex orientation on $E/U$. Section 2 gives opposite signs. This calculation requires the Fredholm identification only on smooth fibers. The comparison on the singular fiber used dualizing sheaves, the residue theorem and the two matching conditions.

## 4. Smoothing roots, automorphisms and the scope of the test

For each nonzero $t$ the framed plumbing has two roots. The simultaneous deck involution sends $(u,v)$ to $(-u,-v)$, preserving $s$ and both graphs. It acts by $+1$ on the two pulled-back Ramond frames and by $-1$ on both modified exceptional frames. On the smooth spin determinant its action preserves the continued complex orientation, by equivariance along the mass homotopy using an invariant metric. Changing one coordinate sign can relabel the two roots and need not extend to a deck transformation of the connected normalization. The comparison is made separately at each root with the same automorphism weight for both spin choices. There is no extra orientation sign from enumerating the roots. Reversing the chosen unit of $\eta$ on the connected cut complement changes the common residue identification by $-1$ on every summand, of total even dimension, hence preserves the comparison.

For a concrete positive-rank case take $B$ elliptic and $A_0$ a nontrivial two-torsion line. Form a connected double cover branched simply at $p,q$ and two other distinct points, using a square root of their divisor line bundle. This root exists because multiplication by two on $\operatorname{Pic}^0(B)$ is onto. The smoothed target has genus two and both target spins are even. The residual rank is two. Both projection orientations can be compared even though the rational Euler class on the smooth profile stratum vanishes. Indeed every branch point of a smooth double cover is simple, and the vanishing for an even ramification index follows from the odd-rank orthogonal summand in Article V, §§4–6. When $R=0$ the sign recovers the determinant sign of a transposition in the unramified model of Article III, Appendix F.

Lee and Parker pair different smoothing choices at fixed spin data in their cancellation formula for spin Hurwitz numbers. The comparison here varies the target spin at a fixed smoothing of a finite cover and also treats residual simple ramification. [Lee–Parker, §§5 and 10](https://arxiv.org/pdf/1212.1825)

The calculation yields the predicted factor $-1$ on the finite contact chart. A compact relative gluing formula must also treat degenerations of this chart, identify its real normal determinant with that of the expanded or contracted stable map, and pair the obstruction class with the relative virtual class with its stack multiplicities. In particular, on a contracted component $df$ vanishes identically, so the sheaf quotient in §2 no longer has finite support. No integrated boundary coefficient or partition trace follows from the present sign. Calculations D6 and D7 extend the comparison to the compact boundary and evaluate the resulting trace.
