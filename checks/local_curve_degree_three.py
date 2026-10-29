"""Exact BP reconstruction through degree three.

Arithmetic verification for Paper 2; geometric interpretation is proved in Paper 2.
Run with SymPy 1.14.0. The independent quotient-algebra check is separate.
"""
from collections import Counter
from math import factorial
import sympy as S

r, w, T = S.symbols("r w T")

def partitions(n, maximum=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n, maximum or n), 0, -1):
        for tail in partitions(n-first, first):
            yield (first,) + tail

def z(part):
    return S.prod(k**m * factorial(m) for k,m in Counter(part).items())

def canon(count):
    return tuple(sorted((k for k,m in count.items() for _ in range(m)), reverse=True))

def clean(matrix):
    return matrix.applyfunc(S.cancel)

def build(d,t1=1,t2=1):
    parts=list(reversed(list(partitions(d))))
    loc={p:i for i,p in enumerate(parts)}
    n=len(parts)
    M=S.zeros(n)
    for j,part in enumerate(parts):
        count=Counter(part)
        M[j,j]=S.cancel(-S.Rational(1,2)*(t1+t2)*sum(
            k*m*(k*(r**k+1)/(r**k-1)-(r+1)/(r-1)) for k,m in count.items()))
        for k in range(1,d):
            for l in range(1,d-k+1):
                if count[k+l]:
                    target=count.copy()
                    target[k+l]-=1
                    target[k]+=1
                    target[l]+=1
                    dest=canon(target)
                    M[loc[dest],j] += -S.Rational(1,2)*t1*t2*(k+l)*count[k+l]*z(dest)/z(part)
                multiplicity=count[k]*(count[l]-(k==l))
                if multiplicity:
                    target=count.copy()
                    target[k]-=1
                    target[l]-=1
                    target[k+l]+=1
                    dest=canon(target)
                    M[loc[dest],j] += S.Rational(1,2)*k*l*multiplicity*z(dest)/z(part)
    M=clean(M)
    identity=S.eye(n)[:,0]
    powers=[S.eye(n)]
    for _ in range(n-1):
        powers.append(clean(powers[-1]*M))
    krylov=S.Matrix.hstack(*(P*identity for P in powers))
    krylov_inv=clean(krylov.inv())
    multiplication=[]
    for j in range(n):
        L=sum((krylov_inv[i,j]*powers[i] for i in range(n)),S.zeros(n))
        multiplication.append(clean(L))
    D=S.diag(*(S.Integer(-1)**len(p)/(z(p)*S.Integer(t1*t2)**len(p)) for p in parts))
    for i,L in enumerate(multiplication):
        assert clean(L*identity-S.eye(n)[:,i]) == S.zeros(n,1)
        assert clean(D*L-(D*L).T) == S.zeros(n)
    G=clean(sum((multiplication[i]**2/D[i,i] for i in range(n)),S.zeros(n)))
    cap1=S.Matrix([(-1)**d*S.Integer(t1)**len(p)*S.prod(1/(1-r**(-k)) for k in p) for p in parts])
    cap2=S.Matrix([(-1)**d*S.Integer(t2)**len(p)*S.prod(1/(1-r**(-k)) for k in p) for p in parts])
    A=clean(sum((cap1[i]*multiplication[i] for i in range(n)),S.zeros(n)))
    Abar=clean(sum((cap2[i]*multiplication[i] for i in range(n)),S.zeros(n)))
    B=clean(r**d*G*A.inv()*Abar.inv())
    return parts,M,D,G,A,Abar,B

def in_w(value):
    num,den=S.fraction(S.cancel(value))
    modulus=S.Poly(r**2-(2-w)*r+1,r,domain=S.QQ.frac_field(w))
    num=S.Poly(num,r,domain=S.QQ.frac_field(w)).rem(modulus)
    den=S.Poly(den,r,domain=S.QQ.frac_field(w)).rem(modulus)
    result=(num*den.invert(modulus)).rem(modulus)
    assert result.degree() <= 0, result
    return S.factor(result.as_expr())

if __name__ == "__main__":
    for d in (1,2,3):
        for weights in ((1,1),(1,-1)):
            print("Starting",d,weights,flush=True)
            parts,M,D,G,A,Abar,B=build(d,*weights)
            print("Basis",parts,flush=True)
            if d==3 and weights==(1,1):
                print("M =",M,flush=True)
                print("B =",B,flush=True)
            print("trace B =",in_w(S.trace(B)),flush=True)
            print("characteristic polynomial =",[in_w(c) for c in B.charpoly(T).all_coeffs()],flush=True)
