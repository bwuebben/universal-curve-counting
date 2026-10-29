"""Rational Frobenius periods and genus-zero GV invariants for resolved X9.

Run: sage -c "load('checks/check_x9_mirror_periods.sage')"
Optional degree: MIRROR_DEGREE=6 sage -c "load('checks/check_x9_mirror_periods.sage')"
Implements the Gamma-series formula and master equation (4.30) of
Demirtas et al., arXiv:2303.00757, with consistent total-degree truncation.
It computes a large-volume expansion, not analytic continuation to the face.
"""
import contextlib
import io
import os
from functools import lru_cache

with contextlib.redirect_stdout(io.StringIO()):
    load('checks/check_x9_prepotential.sage')

nef = matrix(ZZ,sorted([list(r.vector()) for r in S.rays()])).T
assert abs(nef.det()) == 1
glsm = nef.inverse()*matrix(ZZ,[CLS[n] for n in names[:-1]]).T
qtotal = sum(glsm.columns())
assert glsm*matrix(ZZ,pts[:-1]) == 0
assert qtotal == vector(ZZ,[0,0,2,1,0,0])
assert nef.T*C2 == vector(ZZ,[1,0,0,0,0,0])
assert nef.T*C1 == vector(ZZ,[0,0,0,0,1,0])
assert nef.inverse()*J0 == vector(ZZ,[1]*6)

# Removing the facet-interior ray must retain the same chamber on X.
# Reconstruct the coarser ambient fan, rather than assuming this from the
# rank of the GLSM matrix. Its quotient singularities miss the hypersurface.
coarse_simplices = set()
for normal,on_full in facets:
    on = [i for i in on_full if i!=10]
    if len(on)==4:
        coarse_simplices.add(tuple(sorted(on)))
        continue
    omit = next(k for k in range(4) if normal[k])
    keep = [k for k in range(4) if k!=omit]
    lifted = [tuple([pts[i][k] for k in keep]+[H[i]]) for i in on]
    poly = Polyhedron(vertices=lifted)
    for face in poly.faces(3):
        inequalities = [e for e in face.ambient_Hrepresentation() if e.is_inequality()]
        assert len(inequalities)==1
        if inequalities[0].A()[-1]>0:
            indices = tuple(sorted(on[lifted.index(tuple(v))] for v in face.vertices()))
            assert len(indices)==4
            coarse_simplices.add(indices)
coarse_fan = Fan(cones=sorted(coarse_simplices),rays=pts[:-1],check=True)
assert coarse_fan.is_complete() and coarse_fan.is_simplicial()
coarse_walls = {}
for simplex in coarse_simplices:
    for wall in combinations(simplex,3):
        coarse_walls.setdefault(wall,[]).append(next(i for i in simplex if i not in wall))
wall_relations = []
for wall,opposite in coarse_walls.items():
    assert len(opposite)==2
    indices = list(wall)+opposite
    relation = matrix(QQ,[pts[i] for i in indices]).T.right_kernel().basis()[0]
    if relation[3]<0:
        relation = -relation
    assert relation[3]>0 and relation[4]>0
    expanded = [QQ(0)]*10
    for i,x in zip(indices,relation):
        expanded[i] = x
    wall_relations.append([0]+expanded)
coarse_nef_raw = Polyhedron(ieqs=wall_relations)
class_map = matrix(QQ,[CLS[n] for n in names[:-1]]).T
coarse_nef = Polyhedron(rays=[class_map*r.vector() for r in coarse_nef_raw.rays()])
assert coarse_nef == S

@lru_cache(None)
def harmonic(n,p=1):
    return sum((QQ(1)/QQ(k)^p for k in range(1,int(n)+1)),QQ(0))

def gamma_coefficients(Q,n):
    """(c, grad c, rational Hessian c), with the universal zeta(2) removed.

    Equivalently differentiate the rational function c(n+rho)/c(rho).
    Negative factorial arguments are zeros of reciprocal Gamma, not
    terms to discard before taking first and second derivatives.
    """
    r = Q.nrows()
    columns = Q.columns()
    qt = sum(columns)
    ks = [vector(ZZ,n)*column for column in columns]
    k0 = vector(ZZ,n)*qt
    assert k0>=0
    negative = [i for i,k in enumerate(ks) if k<0]
    if len(negative)>2:
        return QQ(0),zero_vector(QQ,r),zero_matrix(QQ,r,r)
    pref = QQ(factorial(k0))/prod(factorial(k) for k in ks if k>=0)
    for i in negative:
        pref *= (-1)^(int(-ks[i])-1)*factorial(-ks[i]-1)
    if len(negative)==2:
        u,v = [columns[i] for i in negative]
        return QQ(0),zero_vector(QQ,r),pref*(u.column()*v.row()+v.column()*u.row())
    gradient = qt*harmonic(k0)-sum((column*harmonic(k) for column,k in zip(columns,ks)
                                   if k>=0),zero_vector(QQ,r))
    if len(negative)==1:
        i = negative[0]
        u = columns[i]
        gradient -= u*harmonic(-ks[i]-1)
        return QQ(0),pref*u,pref*(u.column()*gradient.row()+gradient.column()*u.row())
    hessian_log = -harmonic(k0,2)*(qt.column()*qt.row())
    hessian_log += sum((harmonic(k,2)*(column.column()*column.row())
                        for column,k in zip(columns,ks)),zero_matrix(QQ,r,r))
    return pref,pref*gradient,pref*(gradient.column()*gradient.row()+hessian_log)

def multiply_linear_jet(jet,constant,linear):
    c,g,h = jet
    return (constant*c,constant*g+c*linear,
            constant*h+g.column()*linear.row()+linear.column()*g.row())

def check_contiguity(Q,max_degree):
    """Independent Gamma shift recurrences, including negative arguments."""
    count = 0
    qt = sum(Q.columns())
    for total in range(1,max_degree+1):
        for n in IntegerVectors(total,Q.nrows()):
            n = vector(ZZ,n)
            ks = n*Q
            k0 = n*qt
            for a in range(Q.nrows()):
                if n[a]==0:
                    continue
                previous = vector(ZZ,list(n))
                previous[a] -= 1
                left = gamma_coefficients(Q,n)
                right = gamma_coefficients(Q,previous)
                for i,qi in enumerate(Q.row(a)):
                    if qi>0:
                        for m in range(qi):
                            left = multiply_linear_jet(left,ks[i]-m,Q.column(i))
                    elif qi<0:
                        for m in range(1,-qi+1):
                            right = multiply_linear_jet(right,ks[i]+m,Q.column(i))
                for m in range(1,qt[a]+1):
                    right = multiply_linear_jet(right,k0-qt[a]+m,qt)
                assert left==right, ('Gamma contiguity',tuple(n),a)
                count += 1
    return count

def periods_and_gv(Q,kappas,degree):
    r = Q.nrows()
    series = PowerSeriesRing(QQ,names=tuple('z'+str(i+1) for i in range(r)),
                            default_prec=degree+1)
    z = series.gens()
    exponents = [tuple(n) for d in range(degree+1) for n in IntegerVectors(d,r)]
    zero = (0,)*r
    coefficient = lambda p,n: p[n[0] if r==1 else n]
    w0 = series(0).add_bigoh(degree+1)
    cg = [series(0).add_bigoh(degree+1) for _ in range(r)]
    contracted = [series(0).add_bigoh(degree+1) for _ in range(r)]
    for n in exponents:
        c,grad,hess = gamma_coefficients(Q,n)
        mon = prod(z[i]^n[i] for i in range(r))
        w0 += c*mon
        for a in range(r):
            cg[a] += grad[a]*mon
            contracted[a] += sum(kappas[a][b,c0]*hess[b,c0]
                                  for b in range(r) for c0 in range(r))*mon/2
    mirror_log = [c/w0 for c in cg]
    magnetic = [contracted[a]/w0-sum(kappas[a][b,c]*mirror_log[b]*mirror_log[c]
                                    for b in range(r) for c in range(r))/2
                for a in range(r)]
    remainder = list(magnetic)
    gv = {}
    for n in exponents:
        if n==zero:
            assert all(coefficient(p,zero)==0 for p in remainder)
            continue
        residual = vector(QQ,[coefficient(p,n) for p in remainder])
        j = next(i for i in range(r) if n[i])
        invariant = residual[j]/n[j]
        assert residual == invariant*vector(QQ,n), ('integrability',n,list(residual))
        assert invariant in ZZ, ('nonintegral GV',n,invariant)
        gv[n] = ZZ(invariant)
        if invariant:
            qn = prod(z[i]^n[i] for i in range(r))*sum(n[i]*mirror_log[i] for i in range(r)).exp()
            dilog = sum(qn^k/QQ(k)^2 for k in range(1,degree//sum(n)+1))
            remainder = [remainder[a]-invariant*n[a]*dilog for a in range(r)]
    assert all(p==0 for p in remainder)
    return {'ring':series,'w0':w0,'mirror_log':mirror_log,'magnetic':magnetic,'gv':gv}

# Independent one-parameter benchmark exercises mirror map and multicovers.
quintic = periods_and_gv(matrix(ZZ,[[1]*5]),[matrix(QQ,[[5]])],3)
assert [quintic['gv'][(i,)] for i in [1,2,3]] == [2875,609250,317206375]
print('Quintic benchmark: 2875, 609250, 317206375.')
print('Gamma coefficient recurrences through degree 3:',check_contiguity(glsm,3))

degree = int(os.environ.get('MIRROR_DEGREE','6'))
assert degree>=2
nef_kappa = [matrix(QQ,6,6,lambda b,c:trip(nef.column(a),nef.column(b),nef.column(c)))
            for a in range(6)]
data = periods_and_gv(glsm,nef_kappa,degree)
gv = data['gv']
unit_gv = [gv[tuple(1 if j==i else 0 for j in range(6))] for i in range(6)]
assert unit_gv == [1,1,-2,0,1,1]
assert gv[(1,0,0,0,1,0)] == 0
if degree>=3:
    assert gv[(1,0,0,0,1,1)] == 1
    assert nef.T*curve(4,9) == vector(ZZ,[1,0,0,0,1,1])
    assert nef.T*curve(3,9) == vector(ZZ,[0,0,0,0,0,1])
for n,v in gv.items():
    if all(n[i]==0 for i in [1,2,3,5]):
        assert v == (1 if sum(n)==1 else 0)
if degree>=6:
    assert data['w0'].add_bigoh(7) == 1+6*prod(data['ring'].gens()[i]^n
                                               for i,n in enumerate([0,0,1,1,2,2]))
    assert sum(v!=0 for n,v in gv.items() if sum(n)<=6) == 41
if degree>=8:
    assert sum(v!=0 for n,v in gv.items() if sum(n)<=8) == 80

# Independently assemble the printed flat-coordinate Li2 expansion.
flat = PowerSeriesRing(QQ,names=tuple('q'+str(i+1) for i in range(6)),default_prec=3)
q1,q2,q3,q4,q5,q6 = flat.gens()
leading = [flat(0).add_bigoh(3),flat(0).add_bigoh(3)]
for n,v in gv.items():
    if v and sum(n)<=2:
        mon = prod(flat.gen(i)^n[i] for i in range(6))
        covers = sum(mon^k/QQ(k)^2 for k in range(1,2//sum(n)+1))
        leading[0] += v*(-n[0]+n[1])*covers
        leading[1] += v*(n[0]-n[1]-n[4]+n[5])*covers
assert leading[0] == -q1+q2+(-q1^2+q2^2)/4+3*q2*q3-2*q2*q6
assert leading[1] == q1-q2-q5+q6+(q1^2-q2^2-q5^2+q6^2)/4-3*q2*q3+q3*q6
print('Integral nef basis (columns):\n',nef)
print('GLSM charge matrix, columns D0,...,D8,E:\n',glsm)
print('Total degree',degree,':',len(gv),'nonzero classes tested;',sum(v!=0 for v in gv.values()),'nonzero GV invariants.')
for n,v in gv.items():
    if sum(n)<=2:
        print('GV',n,'=',v)
print('Fundamental period:',data['w0'])
print('ALL PERIOD/GV CHECKS PASSED; no strong-coupling continuation asserted.')

def truncated_wall(gv,cutoff,a=0,b=0,eps=1,lam=1,label='F2',digits=50):
    """Numerical root from GV classes up to cutoff, with their full Li2 covers.

    Comparing cutoffs is an empirical convergence diagnostic, NOT a rigorous
    bound for the omitted curve classes or analytic continuation.
    """
    import mpmath as mp
    assert label in ['F2','FR']
    assert a>=0 and b>=0 and eps>0 and lam>0
    assert 1<=cutoff<=max(sum(n) for n in gv)
    with mp.workdps(digits):
        a,b,eps,lam = [mp.mpf(str(QQ(x).numerator()))/mp.mpf(str(QQ(x).denominator()))
                       for x in [a,b,eps,lam]]
        va = (3*a*a+6*a*b+2*b*b+42*a*eps+34*b*eps+121*eps*eps)/2
        vv = (11*a*a+10*a*b+2*b*b+90*a*eps+38*b*eps+176*eps*eps)/2
        f = 3*a+2*b+13*eps
        g = -(3*a+2*b+11*eps)/2
        h = (3*a+2*b+(12 if label=='F2' else 14)*eps)/2
        sp = (lam*lam*(va*h-vv*g)+5*g/12+h/2)/(lam*lam*f*vv+h-5*f/12)
        terms = [(n,int(v)) for n,v in gv.items() if v and sum(n)<=cutoff]
        areas = [eps,eps,b+eps,a+eps,eps,eps]
        def determinant(s,return_parts=False):
            psiA,psiV = mp.mpc(0),mp.mpc(0)
            for n,v in terms:
                arg = mp.exp(-2*mp.pi*(lam*sum(n[i]*areas[i] for i in range(6))+mp.j*s*n[0]))
                weight = v*mp.polylog(2,arg)/(4*mp.pi*mp.pi)
                psiA += (-n[0]+n[1])*weight
                psiV += (n[0]-n[1]-n[4]+n[5])*weight
            zg = lam*lam*va-s+mp.mpf('0.5')+mp.j*lam*(g+s*f)+psiA
            zh = lam*lam*vv-mp.mpf(5)/12+mp.j*lam*h+psiV
            return (zg,zh) if return_parts else (zg*mp.conj(zh)).imag
        scalar_determinant = lambda x: determinant(x)
        root = mp.findroot(scalar_determinant,(sp-mp.mpf('0.001'),sp+mp.mpf('0.001')))
        zg,zh = determinant(root,True)
        assert 0<root<1 and (zg*mp.conj(zh)).real>0
        assert mp.diff(scalar_determinant,root)>0
        assert abs(determinant(root))<mp.mpf(10)**(-digits+8)
        return tuple(mp.nstr(x,25) for x in [sp,root,root-sp])
