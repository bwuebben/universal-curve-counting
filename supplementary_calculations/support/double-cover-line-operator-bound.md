# A4. A zero bound for double covers and real line operators

This calculation concerns reduced graphs on split real lines. A1 proves the extension to every covering degree.

## 1. Statement

Let \(C\) be a compact connected Riemann surface of genus \(g\), let \(L\to C\)
be a complex line bundle of degree \(d\), and let \(D\) be a smooth real
Cauchy--Riemann operator on \(L\). Absorb its complex-linear zeroth-order
term into a holomorphic structure on \(L\), so that
\[
                         D=\bar\partial_L+R
\]
with \(R\) conjugate-linear. Let \(f:Y\to C\) be a finite holomorphic
double cover with \(Y\) smooth and connected, deck involution \(\iota\),
and branch divisor \(B\). Thus
\[
 f_*\mathcal O_Y=\mathcal O_C\oplus M^{-1},\qquad
 M^2\simeq\mathcal O_C(B),\qquad
 \deg M=r,\quad \deg B=2r.
\]
The pullback operator \(D_f\) includes pullback of the \((0,1)\)-form
factor of \(R\), as in Article III, §14.

**Proposition.** If \(0\ne u\in\ker D_f\) is odd,
\(\iota^*u=-u\), then
\[
                              r\le d.
\]
Writing \(w\) for the tautological section of \(f^*M\), the quotient
\(u/w\) descends to a continuous section \(v\) of \(LM^{-1}\).
It belongs locally to \(W^{1,p}\) for every finite \(p>1\), solves
a real Cauchy--Riemann equation with bounded measurable zeroth-order
coefficient, and has precisely \(d-r\) zeros counted with their positive
integer orders.

In particular, for \(d=g-1\), any smooth double cover carrying a nonzero
odd kernel section has
\[
                         g(Y)=2g-1+r\le 3g-2.
\]
This concerns the odd kernel. A section pulled back from \(\ker D\)
is invariant and is not subject to this ramification bound.

## 2. A short proof of the inequality upstairs

The scalar similarity principle implies that a nonzero solution of a
real Cauchy--Riemann equation on a complex line has isolated zeros of
positive integer order. Their total order is the degree of the line
bundle: this is the usual winding-number description of its first
Chern class. Applied on \(Y\), the total is
\[
                            \deg f^*L=2d.
\]
Every ramification point is fixed by \(\iota\), so oddness forces \(u=0\)
there. There are \(2r\) such points. Therefore \(2r\le 2d\).

This proof already establishes the genus bound without dividing by
the ramification coordinate. In a local coordinate \(w\) centered at a
ramification point, smoothness and the similarity principle give
\[
                   u(w)=a w^k+o(|w|^k),\qquad a\ne0.
\]
The equality \(u(-w)=-u(w)\) makes \(k\) odd. Zeros away from ramification
occur in deck pairs with equal orders. Subtracting one zero at each
ramification point, then dividing the remaining total by two, gives
\(d-r\). The next section shows that these are exactly the zeros of a
well-defined weak section downstairs.

The analytic input used here, and below, is the author-written source
[Wendl, *Lectures on Symplectic Field Theory*, Section 2.5,
Theorem 2.5.3, printed pp. 52--53](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf#page=64).
Its hypothesis allows a zeroth-order coefficient in \(L^p\), \(p>2\),
and a weak \(W^{1,p}\) solution. It supplies a continuous invertible
complex change of frame taking the solution to a holomorphic function.
The local regularity step below is
[Lemma 2.4.7, printed p. 41](https://www2.mathematik.hu-berlin.de/~wendl/Sommer2026/SFT/lecturenotes.pdf#page=53):
an \(L^1\) distributional solution of \(\bar\partial v=F\), with
\(F\in L^p\), is locally \(W^{1,p}\).

## 3. The downstairs equation and its extension at a branch point

Choose holomorphic local frames of \(L\) and \(M\), and write the
branch section as \(s\), so \(w^2=s\). Away from \(B\), oddness makes
\[
                              u=w\,f^*v
\]
for a unique section \(v\) of \(LM^{-1}\). If locally
\(D=\bar\partial+R(z)\,\mathrm{conj}\), its equation is
\[
 \bar\partial v+
 R(z)\frac{\bar w}{w}\bar v
 =
 \bar\partial v+
 R(z)\frac{|s(z)|}{s(z)}\bar v=0.                 \tag{1}
\]
The displayed equation follows after canceling the nonzero derivative
of \(f\) on the unramified locus. The coefficient has absolute value
\(|R(z)|\).

This is an equation on \(LM^{-1}\), not a scalar expression depending
on a chosen square root. For example, if the frame of \(M\) changes by
\(m'=h\,m\), then \(w'=w/h\) and \(v'=h\,v\). Consequently the
antilinear coefficient in (1) changes by \(h/\bar h\), which is the
correct change for the coefficient of \(\bar v'\). Changes of the
frame of \(L\) and of the base coordinate obey the corresponding
usual tensor rules. Thus (1) defines a bounded measurable
antilinear bundle map globally; its values on the finite set \(B\)
may be assigned arbitrarily.

Near a branch point choose coordinates with \(z=w^2\). The actual
pulled-back equation, including its differential factor, is
\[
 \partial_{\bar w}u+
          2\bar w\,R(w^2)\bar u=0.              \tag{2}
\]
This operator has smooth coefficients. A weak kernel section is
therefore smooth. Oddness and Taylor expansion give
\[
                   u(w)=a w+b\bar w+O(|w|^3).
\]
Evaluating (2) at \(w=0\) gives \(b=0\). Hence
\[
                    \frac{u(w)}{w}
                          =a+O(|w|^2)
                          =a+O(|z|).
\]
The quotient is invariant under \(w\mapsto-w\), so it descends
continuously across \(z=0\). In particular it is bounded.

On a punctured disk it solves (1) classically. Integrate by parts on
the complement of a disk of radius \(\epsilon\) against a smooth
test function. The inner boundary term has size at most a constant
times
\(\epsilon\|v\|_\infty\), and tends to zero. The bounded right-hand
side of (1) is integrable. Thus (1) holds distributionally across
the puncture. Lemma 2.4.7 then gives
\[
                         v\in W^{1,p}_{\rm loc}
                       \quad(1<p<\infty).
\]
There is no delta distribution at the branch point.

Apply the scalar similarity principle to (1). Locally \(v=\Phi h\),
where \(\Phi\) is continuous and nowhere zero and \(h\) is holomorphic.
Unless \(v\) vanishes identically, its zeros are isolated, their
orders are positive, and their total is
\[
                       c_1(LM^{-1})[C]=d-r.
\]
If \(v\) has order \(k\) at a branch point, then \(u\) has order
\(2k+1\) upstairs. At an unramified zero of \(v\), each of the two
lifts has order \(k\). These statements also verify directly the
count in Section 2.

## 4. Genus two and split normal graphs

Now suppose \(g=2\), \(N=L_1\oplus L_2\), \(\deg L_j=1\), and
\[
                              D=D_1\oplus D_2
\]
is any pair of smooth real line operators. The operators need not
be the special spin endpoint operators. The complex-linear parts
may define different holomorphic structures from those initially
chosen on the two smooth lines.

For a connected smooth double cover \(f:Y\to C\), let
\[
                    F=(f,u_1,u_2):Y\longrightarrow N,
             \qquad (D_j)_f u_j=0.
\]
Decompose each section into its invariant and odd parts using
\((1\pm\iota^*)/2\). Invariant kernel sections descend to smooth
kernel sections on \(C\): they are bounded near a branch point,
solve the base equation off it, and the same distributional
removal and elliptic regularity argument applies.

If both odd parts vanish, \(F\) factors through \(f\) and a
degree-one graph. It is then not a simple degree-two graph.
If some odd part is nonzero, the proposition gives \(r=0\) or \(1\).

### One pair of branch points

If \(r=1\), choose a nonzero odd component. Its descended section
has degree zero and no zeros. It separates the two sheets everywhere
off the branch divisor. At each branch point its upstairs expression
is \(a w+O(|w|^3)\), with \(a\ne0\), so the vertical derivative of
the graph is nonzero there. The invariant part has zero derivative
in the \(w\)-coordinate at that point. Away from branch points,
the base derivative is already nonzero.

It follows that \(F\) is an injective immersion, hence an embedding,
and its domain has genus
\[
                                 g(Y)=4.
\]

### An étale double cover

If \(r=0\), \(F\) is an immersion because \(df\) never vanishes.
Choose a nonzero odd component. Its descended section has degree
one, hence exactly one simple zero. A pair of sheets can coincide
under \(F\) only over this point; coincidence also requires the
other odd component to vanish there.

At a zero of a solution of a smooth scalar real Cauchy--Riemann
equation, the equation gives \(\bar\partial v=0\) at that point.
A simple zero therefore has a nonzero complex-linear derivative.
The difference of the two local graph sections has this nonzero
derivative in the chosen component. Thus, if a sheet collision
occurs, the two image branches have distinct tangent lines.

The normalization has genus \(3\). There is at most one ordinary
double point in this parametrized sense. Identifying the two
preimages of such a collision produces a nodal complex source of
arithmetic genus \(4\); if there is no collision, the image is a
smooth embedded genus-three curve. This description does not
invoke an algebraic scheme structure on an almost-complex image.

### Two degree-one components

For completeness, consider two distinct kernel graphs over \(C\).
Their difference is a section in the kernel of \(D_1\oplus D_2\).
Some component is nonzero and, by the ordinary scalar zero theorem,
has exactly one simple zero. The two graphs therefore meet at
most once, with distinct tangent lines if they meet.

A connected reduced union must have that one intersection. Its
nodal source has two genus-two components and one node, hence
arithmetic genus \(4\).

Consequently every reduced connected degree-two nodal curve in this
exact split linear graph model, with no contracted components and
whose parametrization is generically one-to-one onto its image, has
arithmetic genus at most four. For an irreducible curve its smooth
normalization is covered by the two cases above. For a reducible
curve, the only degree partition is \(1+1\), covered by the last
case.
