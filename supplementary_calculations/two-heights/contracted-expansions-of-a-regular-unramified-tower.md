# F2. Contracted expansions of a regular unramified tower

Let $C$ be a closed curve of genus $g\ge2$ and $D_t$ a smooth fixed-symbol path on a complex rank-two bundle $N$ of degree $2g-2$. Let $f:Y\to C$ be a smooth connected cover with simple ramification and let $\xi\in\ker f^*D_0$ have embedded graph $\gamma=(f,\xi)$ in the homogeneous graph almost complex structure. Write $h=g(Y)$, $A=D^N_\gamma$, and fix a connected unramified cover $p:X\to Y$ of degree $e$. Put $u=fp$, $d=\deg u$, $H=g(X)=1+e(h-1)$ and $\kappa=p^*\xi$. Assume $u^*D_0$ has kernel $\mathbb R\kappa$ and the derivative $W_u$ on the smooth fixed-target Hurwitz chart and increasing time, evaluated on this line, is an isomorphism onto its vertical cokernel. Denote its absolutely oriented sign by $\eta_u$ as in [Calculation E1](../polar/the-normal-twist-at-a-ramified-quartic-tower.md).

For $k\ge H$ let $G_{u,k}$ be the image of $\overline{\mathcal M}_k(X,1)$ by composition with $u$, together with the projective line of the pullback of $\kappa$ at $t=0$. Its nonconstant normalization is the specified $u$. Every maximal connected contracted subcurve attaches at one point of this normalization, since a degree-one map to $X$ is an isomorphism on its nonconstant normalization. This includes collisions of attachment points and every stable specialization on a contracted block.

**Proposition.** In each fixed genus $G_{u,k}$ is compact and isolated in the projective kernel incidence. Its contribution to the connected absolute normal Euler change is
$$
 \frac{2\eta_u}{|\operatorname{Aut}_C(u)|}\,[v^{2k-2}]x^{H-1},
 \qquad x=(2\sin(v/2))^2.                                             \tag{1}
$$
The isolating neighborhood and time interval may depend on $k$. A contracted block joining different sheets of $p$ is outside $G_{u,k}$.

## Lifting the graph and removing the contracted components

Take a relatively compact tubular neighborhood $U$ of $\gamma$ in $N$. Pull $p$ back along a smooth retraction to obtain a finite covering $q:\widetilde U\to U$. The inverse image of $\gamma$ is one embedded copy of $X$, and the composite $\pi\circ q$ with the bundle projection is holomorphic on the entire tube for the pulled-back graph structures. The tube retracts onto $X$, so projection to $C$ sends its integral homology generator to degree $d$. Every nonconstant closed holomorphic component in the tube has positive projection degree divisible by $d$. Indeed a component with constant projection maps by $q$ to a single complex vector fiber, where it is constant. These assertions persist for the projection-preserving perturbation below, which keeps the standard vertical complex structure.

Properness of the degree-one stable-map stack gives compactness of the central image. Its vertical kernel is exactly the indicated line by restriction to $X$ and matching of constants on each contracted block. Write $G=\operatorname{Aut}_C(u)$. The map to the projective incidence identifies central objects precisely by the action of $G$ on their degree-one target $X$. An identification of two underlying stable maps induces an automorphism of the fixed positive normalization over $C$ and conversely. The automorphism preserves the projective kernel by its uniqueness. Thus the central stack is $[\overline{\mathcal M}_k(X,1)/G]$, retaining all automorphisms of the contracted curves.

Choose a small disk of regular values of $u$ where $\kappa$ is nonzero on a chosen sheet, and order the $d$ degree-one pieces over this disk by their sheets. A point fiber may also contain a contracted subcurve. Over a smaller disk each labeled component still has degree one near the central family. Contracted curves mapping into the disk belong to the same labeled piece as their attaching point, as follows by domain collapse over a larger disk and degree at a regular value. One can use the whole $d!$-sheeted ordering cover. Its restriction to the central stack has $d!/|G|$ copies of the degree-one family, since the deck action on orderings is free. On each copy the positive normalization has a fixed identification with $X$. Fix a real linear normalization $\ell_a$ by integrating a smooth real covector on the first selected sheet against pullback area, with $\ell_u(\kappa)=1$. Contracted attachments contribute zero to this integral. Use the same sheet under the identification on every copy, transporting choices by the corresponding labels.

For a central normalized graph the fixed identification gives a lift to $X\subset\widetilde U$, with each contracted block lifted constantly at its attachment. On a local universal marked domain $E$ the map $1\times q$ is injective on a neighborhood of the graph of this chosen lift in $E\times\widetilde U$. This follows by compactness of the graph and local injectivity in the target variable. In the relative graph topology the graphs of nearby maps, including maps on smoothed domains, lie in the image and hence have a unique nearby lift. The construction commutes with chart identifications by uniqueness. It uses no count of all lifts of a map into a nonnormal cover. [Stable-map relative topology](https://arxiv.org/abs/1910.12264) · [Calculation E3](../polar/a-ramified-double-joined-to-an-unramified-cover.md).

Consider a sequence of projective solutions approaching a point of the compact central locus. Normalize by $\ell_a=1$ and lift locally. There is precisely one positive normalization in each nearby map, by projection degree and the total degree $d$. Forget its nodal flags and consider the unmarked maps from its smooth normalization. They have bounded genus and energy. On a bounded disk bundle containing the tube, a Hermitian connection gives a closed form restricting to the standard fiber form, and addition of a sufficiently large multiple of a base area form tames the graph structures on its closure. Its class comes from $C$ and the pullback area on every closed map is determined by the projection degree. The same form tames the small perturbations below on a fixed smaller lifted tube.

Pass to an unmarked stable limit. Its only nonconstant component maps with degree one onto the embedded $X$, since the original marked sequence converges to the chosen central decoration and its homology class in the tube is $[X]$. A ghost in this unmarked limit would require a critical external branch, two external branches with the same complex tangent, or at least three external branches at a common target point. The embedding excludes each alternative. A self-node is excluded by injectivity at its two preimages. Thus the unmarked positive normalizations converge smoothly to $X$ and have genus $H$. The ghost criterion also applies to a smooth family of target structures by placing the maps at constant parameter in a product with a complex parameter disk, with a block-diagonal almost complex structure. [Ekholm–Shende, Theorem 1.1](https://arxiv.org/pdf/2212.05835#page=2).

The nearby positive graphs are consequently embedded upstairs. Each connected contracted block has a single positive attachment, as all its attaching points have the same lifted target value. Deleting the blocks leaves a smooth normalized solution on the finite-map chart. Over a smoothed source the horizontal deformation complex $[TX\to u^*TC]$ has torsion cokernel at the ramification points and no higher obstruction. Ramification points with the same value give independent summands. The hypothesis on $W_u$ isolates the normalized bare solution on this chart. This proves isolation of the compact central family and excludes any sequence from a differently attached sector. Upper semicontinuity excludes a two-dimensional kernel throughout a sufficiently small neighborhood.

## A displacement fixing the embedded graph

The fixed-time graph deformation sequence upstairs is
$$
 0\longrightarrow\ker L\longrightarrow\ker p^*A\longrightarrow T_u
 \xrightarrow{\partial_\kappa}\operatorname{coker}L
 \longrightarrow\operatorname{coker}p^*A\longrightarrow0,
 \qquad L=u^*D_0.
$$
Here $T_u$ is the complex tangent to the smooth Hurwitz chart. Projection gives the vertical subcomplex, while cancellation of the tangent and domain variations for the embedded lift leaves its normal operator $p^*A$. The branch-only map is injective with cokernel a real line, since adjoining time gives $W_u$. Thus $p^*A$ has kernel and cokernel lines, with kernel the normal radial variation of the graph.

Shrink the test disk about a regular point where this normal field and a nonzero adjoint field are both nonzero. Such points exist by unique continuation. Here vertical tangent projects isomorphically to graph normal. Choose a smooth real zeroth-order tensor $B$ on this sheet pairing the two fields nontrivially after integration. Let $\kappa_i(w)$ be its sheet height. Modify the graph equation in the chosen lifted target neighborhood by
$$
                   \delta B(w)(z-\kappa_i(w)).
$$
Multiply by a vertical cutoff which is one on a smaller neighborhood of the compact part of the sheet meeting the support of $B$. Extend by zero. The modification defines a smooth almost complex structure preserving the projection to $C$ and the standard vertical structure. It fixes the central graph pointwise as a solution. Its normal linearization has nonzero kernel-to-cokernel derivative, since at the graph the derivative factors via vertical projection along the graph tangent and gives the chosen $B$. The perturbed fixed-time graph normal $A_\delta$ is invertible for small nonzero $\delta$.

On normalized projective charts the same homotopy is defined by the linear expression
$$
         L_{a,t}v+\delta B_a\big(v-\ell_a(v)\kappa_i\big),\qquad L_{a,t}=a^*D_t,
$$
with the tensor and forcing pulled back on the selected labeled sheet and zero on the other sheets. Near the normalized central maps the target cutoff there is identically one wherever $B$ is supported. On the affine slice $\ell_a(v)=1$ the displayed expression is exactly the ordinary perturbed graph equation in the lift. The finite-map normalized zero stays at $(u,0,\kappa)$ and stays transverse by openness of the original derivative. The preceding unmarked-limit argument, applied with parameter $\delta$, and this unique bare zero isolate its contracted expansions throughout a sufficiently short homotopy. In particular no boundary contribution is introduced by the displacement.

Write $L^\delta_{a,t}=L_{a,t}+\delta B_a$ for the homogeneous operator of the displaced affine equation. At the central bare map it is injective. A fixed-horizontal kernel vector induces a kernel vector of $A_\delta$. If its normal class is zero, the vertical vector is tangent to the graph and therefore zero over every regular value of $u$, hence everywhere. Restriction and matching then give injectivity on all contracted expansions. The homogeneous operator and the nonzero forcing in the affine equation play different roles here.

For precision about the virtual comparison, take common stabilized horizontal charts for maps to $C$ and lifted graphs, pulling back the same slicing hypersurfaces by $\pi\circ q$. Their horizontal row and stabilization marks coincide. Lift the horizontal obstruction forms using the horizontal splitting of the affine graph equation. With a vertical supplement $B_0 W$ to the range of $L^\delta_{a,t}$, minimal at the central chart point, solve
$$
                 L^\delta_{a,t}v-B_0 w=\delta B_a\kappa_i.
$$
The derivative in $(v,w)$ is invertible by injectivity and the supplement. It stays so on a smaller thickened chart, giving residual rows $s_C,w,\ell_a(v)-1$. The supplement includes test forms on contracted components whenever required. For continuity under smoothing, near a nodal value on the test disk solve the inhomogeneous equation with the sheetwise operator $D_t+\delta B$ on a base subdisk, subtract the pullback primitive with cutoff, and transform the scalar row by the same substitution. The forcing now vanishes near those nodal values and the homogeneous coefficients are smooth coefficients on each labeled sheet. The standard nodal stabilization applies. Common enlargement of the supplements gives the same oriented product Thom class on overlaps. This is the normalized chart comparison in Article IV, §§2–3 and Appendices A–B and [Bai–Swaminathan, §6.2 and Appendix B](https://arxiv.org/html/2106.01206v4). The target lift is a local diffeomorphism in these charts and keeps both horizontal and vertical rows.

The invertible $A_\delta$ gives a unique smooth path of nearby embedded graphs for the unconstrained ordinary equation. The unmarked-limit argument above also applies before imposing the scalar constraint, uniformly on a short interval in this family. Every nearby positive normalization is the moving embedded curve, and each contracted block attaches to it at one point. Thus the ordinary zero locus consists precisely of its contracted expansions. Transversality of the normalized bare zero says that $\ell_a(v)-1$ cuts this path simply. Contracted expansions have the same scalar value as their bare graph. On the ordinary zero locus the scalar is the smooth function of time obtained by evaluating it on this moving curve. Interpolating the scalar row to that function preserves the zero set on a small isolating neighborhood and its oriented class. Cutting by it therefore leaves the ordinary fixed-time degree-one local invariant of the embedded curve, times the oriented scalar intersection. Indeed its fixed-time normal is injective on all degree-one maps by restriction and matching. The ordinary comparison gives $\epsilon(A_\delta)x^{H-1}$, the determinant sign times the degree-one Hodge series. A smaller symplectic tube may be identified with a tube in the projectivization of the symplectic normal bundle plus a trivial line, and the tame structure extended from its closure. This supplies a compact ambient space for the comparison theorem, with unchanged local equations. Its normal bundle has degree $2H-2$ by adjunction, since the ambient first Chern number along the graph is zero. [Zinger, Theorem 1.2](https://arxiv.org/html/0807.0805v3) · [Pandharipande, Theorem 2](https://arxiv.org/abs/math/9811140) · Article IV, Appendices A and D.

In genus $H$ the scalar sign times $\epsilon(A_\delta)$ is the projective bare degree $\eta_u$, by the transverse homotopy. The relative projective Euler class gives half the ordinary Euler change, since the real cokernel on the transverse link is the tangent of an even-dimensional sphere. The horizontal virtual dimension and the real normal index stay even in every genus, so this orientation calibration persists. Finally the calculation gives one copy of the ordinary series per framed central family. The factor $(d!/|G|)/d!$ is $1/|G|$. These observations prove (1).

## The upper permutation representation

Assume now that the lower wall has a one-dimensional kernel and invertible $W_f$ of sign $\eta_f$. Let $V_p$ be the orthogonal complement of the constants in the real flat permutation bundle $p_*\mathbb R$. The proof of [Calculation E1](../polar/the-normal-twist-at-a-ramified-quartic-tower.md) extends to this bundle. Indeed the pushdown of the graph deformation complex splits by $p_*\mathbb R=\mathbb R\oplus V_p$. On the complement it gives
$$
 0\to\ker (f^*D_0)^{V_p}\to\ker A^{V_p}\to T_{\rm red}
 \longrightarrow\operatorname{coker}(f^*D_0)^{V_p}
 \to\operatorname{coker}A^{V_p}\to0.
$$
The reduced horizontal tangent consists of ramification velocities with coefficients in $V_p$, hence has its complex orientation. Time belongs to the constant summand. Thus the inherited line and transversality hypotheses on $u$ hold exactly when $A^{V_p}$ is invertible, and then $\eta_u=\eta_f\epsilon(A^{V_p})$. Deforming the split real complexes to their complex-linear parts identifies the absolute determinant orientations. No regularity assertion for all upper degrees follows from this criterion alone.

If $f$ is the odd lower double cover, suppose these twists are invertible at a fixed upper degree $e$. For any fixed composite $u$ there are precisely two factorizations through the specified $Y$ as maps over $C$, interchanged by its involution $\tau$. A lift to $Y$ is determined on the connected complement of the inverse images of the branch values by its value over one point. Hence $|\operatorname{Aut}_C(u)|$ equals $|\operatorname{Aut}_Y(p)|$ times the order of the subgroup of $\langle\tau\rangle$ stabilizing the isomorphism class of $p$. The involution preserves the twist sign by fiber negation on the graph. Summing (1) over these orbits gives the connected series
$$
 \eta_f\, x^{e(h-1)}
       \sum_{[p:X\to Y]}\frac{\epsilon(A^{V_p})}{|\operatorname{Aut}_Y(p)|}.       \tag{2}
$$
The sum runs over all connected unramified covers of the fixed $Y$ of degree $e$. A nearby index-zero operator $\widehat A$ on $Y$ with positive determinant sign can be obtained by crossing the simple zero of $A$ in a local real zeroth-order direction. The finitely many invertible twists at this degree keep their signs. Each $p^*\widehat A$ is invertible of sign $\epsilon(A^{V_p})$. Consequently the sum in (2) is the least-genus degree-$e$ fixed-target coefficient of $\widehat A$, computed on the finite stack of the unramified covers. The statement uses only this degree and genus for $\widehat A$.

## The exceptional quartic one genus above its minimum

Use the hypotheses of [Calculation E1–096](../polar/first-quartic-contribution-at-a-ramified-double-wall.md) and [Calculation E3](../polar/a-ramified-double-joined-to-an-unramified-cover.md) at a lower double wall with two ramification points. Then $h=2g$, $H=2h-1$, and the first excess genus is $k=2h$. Put
$$
 T=\tfrac12\sum_{0\ne\lambda\in H^1(Y;\mathbb F_2)}\epsilon(A^\lambda),
 \qquad S_C=\sum_{u:T'\to C\ {\rm connected\ unramified\ double}}\epsilon(u^*D_0).
$$
The contribution of the one-ended upper towers in this genus is $-\eta_f T(H-1)/12$. The isolated sectors of Calculation E3, including both distinct joining points and the two zero-height identity components, add $\eta_f\big(S_C(2(H-1)/12+2)+1\big)$. Thus their combined supported weight is
$$
       -\frac{H-1}{12}J_f(H)+\eta_f(2S_C+1),\qquad
               J_f(H)=\eta_f(T-2S_C).                                      \tag{3}
$$
Here $J_f(H)$ is the complete minimum-genus weight of the specified lower incidence in Calculation E2.

Where a supported total for the full exceptional incidence in genus $2h$ is defined, let $R$ denote that total minus (3). This reserves all intersecting residual strata, including the marked zero-height double-cover families, the height ratios between two nonzero graphs and the finite branched upper covers of $Y$ with their nodal limits. The one-ended towers counted in (1) are disjoint from these limits, including a self-node joining the two upper sheets or a contracted block between them.

For the primitive subtraction the already determined degree-two change is $\eta_f F_h(q^2)$. On the elementary genus-$h$ core let $c=(2^{2h}-1)/2$ be its minimal degree-two coefficient. The primitive minimum-genus change contributed by this incidence is $J_f(H)-\eta_fc$. At the next genus the degree-two elementary coefficient is $F_2^{\rm elem}(h)-(H-1)c/12$, where $F_2^{\rm elem}(h)=-(h-1)2^{2h-3}$. Subtracting this proper-divisor contribution and the sine tail of the minimum-genus primitive gives
$$
 \delta n_{4,2h}^{\rm inc}
       = R+\eta_f\big(2S_C+1+(h-1)2^{2h-3}\big).                           \tag{4}
$$
In genus $g=2$ the known additive term in (4) is $\eta_f(2S_C+97)$. Formula (4) identifies the remaining virtual weight needed for the primitive exponent. It assigns no value to $R$. Article V, §§4–6.
