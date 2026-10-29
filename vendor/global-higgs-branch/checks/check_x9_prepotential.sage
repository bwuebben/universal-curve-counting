#!/usr/bin/env sage
"""Exact intersection theory and gravity-decoupling analysis for X_9.

This check is self-contained and uses only exact rational arithmetic in
SageMath.  It verifies:

1. the regular, fine, star triangulation of the boundary of Delta_9 fixed
   in the paper (diagonals (v6,v7) and (v2,v4), star at p), and that the
   resulting ambient toric fourfold is smooth and complete;
2. the triple intersection numbers kappa_{IJK} = D_I.D_J.D_K.(-K) of the
   resolved hypersurface in the basis (D2,D5,D6,D7,D8,E), the second Chern
   numbers c_2.D_I, and chi = -232;
3. that kappa reproduces the two nodal charge rows and the E_2 root row;
4. the K3 fibration Xhat_9 -> P^1 with fiber class D5 ~ D3 + D7 + D8, and
   the surface types of all toric divisors on Xhat_9;
5. the nef cone of the ambient fourfold, the 31 effective toric curve
   classes on Xhat_9, and the local-scale face
   Nef(Xhat_9) n {w.C1 = w.C2 = 0, w|_E = 0} = cone(D0, D1);
6. the annihilators ann(w) = {v : w.v = 0 in H^4} for w = D1, D0 and an
   interior class, the rank of the charge matrix restricted to each, and the
   limiting kinetic norms along J_0 + Lambda*w;
7. the physical-current elimination giving the formal rigid quotient
   XY = beta^4 (the cover-coordinate cubic was corrected on 2026-09-04).

Run from the project root with

    sage checks/check_x9_prepotential.sage
"""

from itertools import combinations, combinations_with_replacement
from hashlib import sha256

# ---------------------------------------------------------------------------
# 1.  Polytope, rays, triangulation
# ---------------------------------------------------------------------------
V = [(1,0,0,0), (0,1,0,0), (-1,-1,0,0), (0,0,1,0), (0,0,0,1),
     (-4,-2,-1,0), (-4,-2,0,-1), (3,1,1,1), (2,1,1,1)]
P_INT = (-2,-1,0,0)       # interior point of the pentagon: E = D_p
F_INT = (-1,0,0,0)        # facet-interior point: D_f misses the hypersurface
pts = [vector(ZZ, x) for x in V + [P_INT, F_INT]]
NR = len(pts)             # 11 rays, indices 0..8 = v_i, 9 = p, 10 = f
names = ['D0','D1','D2','D3','D4','D5','D6','D7','D8','E','Df']

P = Polyhedron(vertices=V)
assert P.dim() == 4 and len(P.integral_points()) == 12
facets = []
for ineq in P.inequalities():
    n = -vector(ineq.A())
    assert ineq.b() == 1                      # reflexive
    facets.append((n, [i for i, x in enumerate(pts) if n * x == 1]))
assert len(facets) == 12

# Heights defining the regular subdivision.  Generic, and chosen so that
#   h6 + h7 < h2 + h3   (diagonal (v6,v7) in the square {v2,v3,v6,v7}),
#   h2 + h4 < h5 + h7   (diagonal (v2,v4) in the square {v2,v4,v5,v7}),
#   p and f low         (star subdivisions at the interior points).
H = {0: 10, 1: 11, 2: 3, 3: 9, 4: 2, 5: 7, 6: 1, 7: 4, 8: 8, 9: -5, 10: -3}
H = {i: QQ(h) + QQ(i + 1) / QQ(1009) for i, h in H.items()}
assert H[6] + H[7] < H[2] + H[3] and H[2] + H[4] < H[5] + H[7]

simplices = set()
for n, on in facets:
    if len(on) == 4:
        simplices.add(tuple(sorted(on)))
        continue
    j = next(k for k in range(4) if n[k] != 0)
    keep = [k for k in range(4) if k != j]
    lifted = [tuple([pts[i][k] for k in keep] + [H[i]]) for i in on]
    Q = Polyhedron(vertices=lifted)
    assert Q.dim() == 4
    for F in Q.faces(3):
        ieq = [e for e in F.ambient_Hrepresentation() if e.is_inequality()]
        assert len(ieq) == 1
        if -vector(ieq[0].A())[-1] < 0:       # lower facet
            vs = [tuple(v) for v in F.vertices()]
            idx = tuple(sorted(on[lifted.index(v)] for v in vs))
            assert len(idx) == 4, "regular subdivision is not a triangulation"
            simplices.add(idx)
simplices = sorted(simplices)
assert len(simplices) == 29
assert set(i for s in simplices for i in s) == set(range(NR))      # fine
assert all(abs(matrix(ZZ, [pts[i] for i in s]).det()) == 1 for s in simplices)

def has_edge(a, b):
    return any(a in s and b in s for s in simplices)

assert has_edge(6, 7) and not has_edge(2, 3)
assert has_edge(2, 4) and not has_edge(5, 7)
PENT = [3, 4, 5, 6, 8]
assert all(has_edge(9, i) for i in PENT)
pent_edges = sorted((a, b) for a, b in combinations(PENT, 2) if has_edge(a, b))
assert pent_edges == [(3, 6), (3, 8), (4, 5), (4, 8), (5, 6)]   # cycle v3 v6 v5 v4 v8

fan = Fan(cones=[list(s) for s in simplices], rays=pts, check=True)
assert fan.is_complete() and fan.is_simplicial() and fan.is_smooth()
X = ToricVariety(fan)

# ---------------------------------------------------------------------------
# 2.  Triple intersections, Chern numbers, charge rows
# ---------------------------------------------------------------------------
HH = X.cohomology_ring()
D = [HH(X.divisor(i)) for i in range(NR)]
aK = sum(D)
c2 = HH(X.Chern_class(2))
c3 = HH(X.Chern_class(3))
assert X.integrate((c3 - c2 * aK) * aK) == -232

kap = {}
for i, j, k in combinations_with_replacement(range(NR), 3):
    kap[(i, j, k)] = X.integrate(D[i] * D[j] * D[k] * aK)

def kappa(i, j, k):
    return kap[tuple(sorted((i, j, k)))]

assert all(kappa(10, j, k) == 0 for j in range(NR) for k in range(NR))

BASIS = [2, 5, 6, 7, 8, 9]
BN = ['D2', 'D5', 'D6', 'D7', 'D8', 'E']
EXPECTED_KAPPA = {
    ('D2','D2','D2'): 8, ('D2','D2','D5'): -2, ('D2','D2','D6'): -3, ('D2','D2','D7'): -2,
    ('D2','D5','D6'): 1, ('D2','D6','D6'): 1, ('D2','D6','D7'): 1,
    ('D5','D6','E'): 1, ('D5','E','E'): -2,
    ('D6','D6','D6'): -1, ('D6','D6','D7'): -1, ('D6','D7','D7'): -1, ('D6','E','E'): -2,
    ('D7','D7','D7'): -3, ('D7','D7','D8'): 4, ('D7','D8','D8'): -6,
    ('D8','D8','D8'): 8, ('D8','D8','E'): -1, ('D8','E','E'): -1, ('E','E','E'): 7,
}
found = {}
for i, j, k in combinations_with_replacement(range(6), 3):
    val = kappa(BASIS[i], BASIS[j], BASIS[k])
    if val != 0:
        found[(BN[i], BN[j], BN[k])] = val
assert found == EXPECTED_KAPPA, found
c2D = [X.integrate(c2 * D[b] * aK) for b in BASIS]
assert c2D == [-4, 24, 26, 18, -4, -2]

def curve(i, j):
    """Class of D_i.D_j.Xhat as a functional on the basis."""
    return vector(QQ, [kappa(i, j, b) for b in BASIS])

C1 = curve(6, 7)                       # exceptional curve of the node (v6,v7)
C2 = curve(2, 4)                       # exceptional curve of the node (v2,v4)
R = curve(4, 9) - curve(3, 9)          # E_2 root class e_1 - e_2 in E
assert list(C1) == [1, 0, -1, -1, 0, 0]
assert list(C2) == [-1, 1, 0, 1, 0, 0]
assert list(R) == [0, 1, -1, 0, 0, 0]
Q = matrix(QQ, [C1, C2, R])            # charge row of v is Q*v

# toric curves of E = S_7 in the cyclic order v3, v6, v5, v4, v8
e_self = {i: kappa(i, i, 9) for i in PENT}
assert e_self == {3: -1, 6: 0, 5: 0, 4: -1, 8: -1}
assert all(kappa(9, 9, i) == -2 - e_self[i] for i in PENT)   # K_E.C = -2 - C^2

# ---------------------------------------------------------------------------
# 3.  Linear equivalence, K3 fibration, surface types
# ---------------------------------------------------------------------------
Mrel = matrix(QQ, [[pts[i][k] for i in range(NR)] for k in range(4)])
piv = [0, 1, 3, 4]
free = [2, 5, 6, 7, 8, 9, 10]
assert Mrel[:, piv].det() == 1
T = -Mrel[:, piv].inverse() * Mrel[:, free]

def to_basis(a):
    """Class on Xhat_9 of the toric divisor sum a_i D_i, in basis coordinates."""
    a = vector(QQ, a)
    out = vector(QQ, [0] * 6)
    for idx, b in enumerate(BASIS):
        out[idx] += a[b]
    for r, pi in enumerate(piv):
        for c, fr in enumerate(free):
            if fr != 10:
                out[BASIS.index(fr)] += a[pi] * T[r, c]
    return out

def unit(i):
    return to_basis([1 if k == i else 0 for k in range(NR)])

CLS = {names[i]: unit(i) for i in range(NR)}
assert list(CLS['D0']) == [1, 4, 4, -3, -2, 2]
assert list(CLS['D1']) == [1, 2, 2, -1, -1, 1]
assert list(CLS['D3']) == [0, 1, 0, -1, -1, 0]
assert list(CLS['D4']) == [0, 0, 1, -1, -1, 0]
assert CLS['Df'] == 0

def trip(u, v, w):
    return sum(kappa(BASIS[I], BASIS[J], BASIS[K]) * u[I] * v[J] * w[K]
               for I in range(6) for J in range(6) for K in range(6))

def L(w):
    """The map v -> w.v in H^4, as the matrix (kappa_{IJK} w^K)."""
    return matrix(QQ, 6, 6, lambda I, J: sum(
        kappa(BASIS[I], BASIS[J], BASIS[K]) * w[K] for K in range(6)))

# K3 fibration: m = -e_3^* is nonnegative on v5 only and nonpositive elsewhere
m = vector(ZZ, (0, 0, -1, 0))
assert [m * x for x in pts] == [0, 0, 0, -1, 0, 1, 0, -1, -1, 0, 0]
assert not any(any(m * pts[i] > 0 for i in s) and any(m * pts[i] < 0 for i in s)
               for s in simplices)
assert CLS['D5'] == CLS['D3'] + CLS['D7'] + CLS['D8']
assert L(CLS['D5']) * CLS['D5'] == 0           # D5.D5 = 0 in H^4: trivial normal bundle
assert not any(has_edge(5, i) for i in (3, 7, 8))

# Surface invariants: K^2 = D^3, chi = c_2.D + D^3, chi(O) = (K^2 + chi)/12
SURF = {}
for i in range(NR - 1):
    K2 = kappa(i, i, i)
    chi = X.integrate(c2 * D[i] * aK) + K2
    SURF[names[i]] = (K2, chi, (K2 + chi) / 12, L(unit(i)).rank())
EXPECTED_SURF = {
    'D0': (83, 229, 26, 4), 'D1': (14, 94, 9, 4), 'D2': (8, 4, 1, 2),
    'D3': (1, 11, 1, 3), 'D4': (0, 12, 1, 4), 'D5': (0, 24, 2, 3),
    'D6': (-1, 25, 2, 4), 'D7': (-3, 15, 1, 4), 'D8': (8, 4, 1, 2),
    'E': (7, 5, 1, 3),
}
assert SURF == EXPECTED_SURF, SURF
assert all(kappa(5, 5, b) == 0 for b in BASIS)          # K_{D5|X} numerically trivial
assert curve(6, 6) == C1                               # K_{D6|X} = C1, the (-1)-curve
assert kappa(2, 4, 4) == -1 and kappa(2, 2, 4) == -1   # C2 = (-1)-section of D2|X = F_1
assert kappa(2, 5, 5) == 0 and kappa(2, 2, 5) == -2     # D2 n D5 = fiber of the ruling
assert kappa(8, 8, 9) == -1                            # R_{v8} = (-1)-section of D8|X = F_1
assert kappa(6, 7, 7) == -1 and kappa(6, 6, 7) == -1   # C1 has normal bundle O(-1)+O(-1)
assert L(CLS['D2']) * CLS['E'] == 0                    # D2 and E are disjoint on Xhat_9
assert Q * CLS['D2'] == vector(QQ, [1, -1, 0])
assert Q * CLS['D3'] == vector(QQ, [1, 0, 1])
assert Q * CLS['D4'] == vector(QQ, [0, -1, -1])
assert Q * CLS['D8'] == 0 and Q * CLS['E'] == 0
assert Q * CLS['D0'] == 0 and Q * CLS['D1'] == 0

# ---------------------------------------------------------------------------
# 4.  Wall curves, nef cones, the local-scale face
# ---------------------------------------------------------------------------
walls = {}
for s in simplices:
    for t in combinations(s, 3):
        walls.setdefault(t, []).append([i for i in s if i not in t][0])
wall_vecs = []
for t, ab in walls.items():
    assert len(ab) == 2
    a, b = ab
    c = matrix(QQ, [pts[i] for i in t]).T.solve_right(-(pts[a] + pts[b]))
    vec = vector(QQ, [0] * NR)
    vec[a] = 1; vec[b] = 1
    for i, ci in zip(t, c):
        vec[i] = ci
    assert sum(vec[i] * pts[i] for i in range(NR)) == 0
    wall_vecs.append(vec)
assert len(wall_vecs) == 58

nefP = Polyhedron(ieqs=[[0] + list(v) for v in wall_vecs])
assert nefP.dim() == 11 and nefP.n_lines() == 4
S_rays = [to_basis(r.vector()) for r in nefP.rays()]
S = Polyhedron(rays=S_rays)                 # restriction of the ambient nef cone
assert S.dim() == 6
assert sorted(map(list, [r.vector() for r in S.rays()])) == sorted([
    [0, 1, 0, 0, 0, 0], [0, 1, 1, -1, -1, 0], [1, 2, 2, -1, -1, 1],
    [1, 4, 4, -3, -2, 2], [2, 5, 4, -3, -2, 2], [2, 5, 5, -3, -2, 2]])
assert S.contains(CLS['D0']) and S.contains(CLS['D1'])

two_cones = sorted(set(t for s in simplices for t in combinations(s, 2)))
assert len(two_cones) == 40
eff = {(i, j): curve(i, j) for (i, j) in two_cones if curve(i, j) != 0}
assert len(eff) == 31
N = Polyhedron(ieqs=[[0] + list(cv) for cv in eff.values()])
assert all(N.contains(r.vector()) for r in S.rays())

LE = L(CLS['E'])
assert LE.rank() == 3
annE = LE.right_kernel()
assert annE.dimension() == 3
eqs = [[0] + list(C1), [0] + list(C2)] + [[0] + list(LE.row(I)) for I in range(6)]
Hsub = Polyhedron(eqns=eqs)
assert Hsub.dim() == 2
faceS = S.intersection(Hsub)
faceN = N.intersection(Hsub)
face_rays = sorted(map(list, [r.vector() for r in faceS.rays()]))
assert face_rays == sorted(map(list, [r.vector() for r in faceN.rays()]))
assert face_rays == sorted([list(CLS['D0']), list(CLS['D1'])])
# two effective curves already cut the quadrant out of the plane span(D0, D1)
assert eff[(1, 8)] * CLS['D1'] == 1 and eff[(1, 8)] * CLS['D0'] == 0
assert eff[(0, 2)] * CLS['D1'] == 0 and eff[(0, 2)] * CLS['D0'] == 4

# ---------------------------------------------------------------------------
# 5.  Dynamical vectors along each decoupling direction
# ---------------------------------------------------------------------------
def ann(w):
    return L(w).right_kernel()

def qrank(space):
    if space.dimension() == 0:
        return 0
    return (Q * matrix(space.basis()).T).rank()

A1 = ann(CLS['D1']); A0 = ann(CLS['D0']); Ai = ann(CLS['D0'] + CLS['D1'])
assert trip(CLS['D1'], CLS['D1'], CLS['D1']) == 14
assert trip(CLS['D0'], CLS['D0'], CLS['D0']) == 83
assert trip(CLS['D0'], CLS['D0'], CLS['D1']) == 49 and trip(CLS['D0'], CLS['D1'], CLS['D1']) == 27
assert A1 == span([CLS['D2'], CLS['E']], QQ)
assert A0 == span([CLS['D8'], CLS['E']], QQ)
assert Ai == span([CLS['E']], QQ)
assert (qrank(A1), qrank(A0), qrank(Ai)) == (1, 0, 0)
# curves contracted by D1: those in D2|X, in E, and C1; by D0: those in D8|X, in E, C1, C2
contracted_D1 = sorted(k for k, cv in eff.items() if cv * CLS['D1'] == 0)
contracted_D0 = sorted(k for k, cv in eff.items() if cv * CLS['D0'] == 0)
assert contracted_D1 == [(0,2),(2,4),(2,5),(2,6),(2,7),(3,9),(4,9),(5,9),(6,7),(6,9),(8,9)]
assert contracted_D0 == [(1,8),(2,4),(3,8),(3,9),(4,8),(4,9),(5,9),(6,7),(6,9),(7,8),(8,9)]
# the K3 fiber grows along both rays
assert trip(CLS['D1'], CLS['D1'], CLS['D5']) == 4 and trip(CLS['D0'], CLS['D0'], CLS['D5']) == 14

# Kinetic norms N_v = -J.v^2 + (J^2.v)^2/(4 Vol), Vol = J^3/6, along J = J0 + Lambda*w
J0 = sum(r.vector() for r in S.rays())       # interior of the ambient nef cone: Kaehler on Xhat_9
print("J0 =", list(J0))
Lam = var('Lam')
def norm_along(w, v):
    J = J0 + Lam * w
    vol = trip(J, J, J) / 6
    return -trip(J, v, v) + trip(J, J, v) ** 2 / (4 * vol)
LIMITS = {}
for name, w in [('D1', CLS['D1']), ('D0', CLS['D0'])]:
    lim = {bn: limit(norm_along(w, CLS[bn]), Lam=oo) for bn in ['D2', 'D8', 'E', 'D5', 'D6']}
    LIMITS[name] = lim
    if name == 'D1':
        assert lim['D2'] == -trip(J0, CLS['D2'], CLS['D2']) > 0
        assert lim['E'] == -trip(J0, CLS['E'], CLS['E']) > 0
        assert lim['D8'] == +oo and lim['D5'] == +oo and lim['D6'] == +oo
        assert trip(J0, CLS['D2'], CLS['E']) == 0                   # block-diagonal limit
    else:
        assert lim['D8'] == -trip(J0, CLS['D8'], CLS['D8']) > 0
        assert lim['E'] == -trip(J0, CLS['E'], CLS['E']) > 0
        assert lim['D2'] == +oo and lim['D5'] == +oo and lim['D6'] == +oo
        blk = matrix(QQ, [[-trip(J0, CLS[a], CLS[b]) for b in ('D8', 'E')] for a in ('D8', 'E')])
        assert blk.is_positive_definite()
# the K3-fiber vector has no J.v^2 term at all: its norm is vol(fiber)^2 / Vol exactly
assert all(trip(w, CLS['D5'], CLS['D5']) == 0 for w in [J0, CLS['D0'], CLS['D1']])
N5 = norm_along(CLS['D1'], CLS['D5'])
assert (N5 - trip(J0 + Lam*CLS['D1'], J0 + Lam*CLS['D1'], CLS['D5'])**2
        / (4 * trip(J0 + Lam*CLS['D1'], J0 + Lam*CLS['D1'], J0 + Lam*CLS['D1']) / 6)).simplify_full() == 0

# ---------------------------------------------------------------------------
# 6.  The formal rigid quotient in physical E2 current operators
# ---------------------------------------------------------------------------
Rg = PolynomialRing(QQ, 'h1,k1,h2,k2,Jp,Jm,t1,t2,be,Xv,Yv')
h1, k1, h2, k2, Jp, Jm, t1, t2, be, Xv, Yv = Rg.gens()
I = Rg.ideal([t1 - h1*k1, t2 - h2*k2, Jp*Jm - be**2,
              Xv - h1*h2*Jm, Yv - k1*k2*Jp, t1 + be, t2 + be])
Jel = I.elimination_ideal([h1, k1, h2, k2, Jp, Jm, t1, t2])
assert Jel == Rg.ideal([be**4 - Xv*Yv])

# ---------------------------------------------------------------------------
# 7.  Report
# ---------------------------------------------------------------------------
payload = repr((tuple(map(tuple, pts)), tuple(simplices),
                tuple(sorted((k, int(v)) for k, v in kap.items()))))
print("X_9 intersection-theory SHA-256:", sha256(payload.encode('ascii')).hexdigest())
print("fan: 29 unimodular cones, smooth complete ambient fourfold; chi(Xhat_9) = -232")
print("nonzero kappa_{IJK} on (D2,D5,D6,D7,D8,E):")
for key, val in EXPECTED_KAPPA.items():
    print("   %s.%s.%s = %d" % (key + (val,)))
print("c_2.D_I =", c2D)
print("charge rows from kappa: C1 =", list(C1), " C2 =", list(C2), " root =", list(R))
print("K3 fibration: m = (0,0,-1,0), fiber class D5 = D3 + D7 + D8, D5.D5 = 0 in H^4")
print("surfaces (K^2, chi, chi(O), rank of restriction of H^2):")
for k, v in SURF.items():
    print("   %s|X: %s" % (k, v))
print("local-scale face of the nef cone: cone(D0, D1), with D0 =", list(CLS['D0']),
      "D1 =", list(CLS['D1']))
print("ann(D1) = span(D2, E), charge rank 1;  ann(D0) = span(D8, E), charge rank 0;"
      "  ann(D0 + D1) = span(E), charge rank 0")
print("Kaehler class J0 =", list(J0))
print("limits along J0 + Lambda D1:", LIMITS['D1'])
print("limits along J0 + Lambda D0:", LIMITS['D0'], " with J0.D8.E =", trip(J0, CLS['D8'], CLS['E']))
print("formal rigid quotient in physical currents: XY = beta^4")
print("All X_9 rigid-limit checks passed.")
