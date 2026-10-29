"""Exact source X9 genus-zero numerical realization on a downward-closed box.
No universal or higher-genus claim. Uses the rational Frobenius Gamma jet and
master equation implemented in canonical check_x9_mirror_periods.sage.
"""
from sage.all import *
from itertools import product, permutations
from functools import lru_cache
from pathlib import Path
import argparse,ast,re,json,hashlib,time
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-root',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
PIN={'check_x9_mirror_periods.sage':'3fdc272a99e7dd506db0b9da66c694eb8fcbee32c0d96e715e5258a5de4949c6','check_x9_prepotential.sage':'c5c25b261f92db49ef4e78d8625e451cac2b01b60c9021b7fa3f96b156ef662b'}
for name,digest in PIN.items():assert hashlib.sha256((args.source_root/'checks'/name).read_bytes()).hexdigest()==digest
Q=matrix(ZZ,[[0,0,-1,0,-1,1,0,1,0,0],[0,0,0,1,1,0,0,0,-1,-1],[0,1,0,0,0,0,0,2,-2,1],[1,0,-1,0,0,0,0,0,0,1],[0,0,1,1,0,0,-1,-1,0,0],[0,0,0,-1,0,0,1,0,1,-1]])
N=matrix(ZZ,[[0,0,1,1,2,2],[1,1,2,4,5,5],[0,1,2,4,4,5],[0,-1,-1,-3,-3,-3],[0,-1,-1,-2,-2,-2],[0,0,1,2,2,2]])
src=(args.source_root/'checks/check_x9_prepotential.sage').read_text()
raw=ast.literal_eval(re.search(r'EXPECTED_KAPPA = (\{.*?\})\n',src,re.S).group(1))
labels=['D2','D5','D6','D7','D8','E']
base={tuple(sorted(labels.index(v) for v in key)):ZZ(val) for key,val in raw.items()}
def trip(a,b,c):return sum(a[i]*b[j]*c[k]*base.get(tuple(sorted((i,j,k))),0) for i in range(6) for j in range(6) for k in range(6))
assert Q.matrix_from_columns([2,5,6,7,8,9])*N == identity_matrix(ZZ,6)
kappa=[matrix(QQ,6,6,lambda b,c:trip(N.column(a),N.column(b),N.column(c))) for a in range(6)]
@lru_cache(None)
def harmonic(n,p=1):return sum((QQ(1)/QQ(k)**p for k in range(1,int(n)+1)),QQ(0))
def gamma_factory(Q):
 r=Q.nrows();cols=Q.columns();qt=sum(cols)
 @lru_cache(None)
 def gamma(n):
  v=vector(ZZ,n);ks=[v*c for c in cols];k0=v*qt;assert k0>=0
  neg=[i for i,k in enumerate(ks) if k<0]
  if len(neg)>2:return QQ(0),zero_vector(QQ,r),zero_matrix(QQ,r,r)
  pref=QQ(factorial(k0))/prod(factorial(k) for k in ks if k>=0)
  for i in neg:pref*=(-1)**int(-ks[i]-1)*factorial(-ks[i]-1)
  if len(neg)==2:
   u,w=[cols[i] for i in neg]
   return QQ(0),zero_vector(QQ,r),pref*(u.column()*w.row()+w.column()*u.row())
  g=qt*harmonic(k0)-sum((c*harmonic(k) for c,k in zip(cols,ks) if k>=0),zero_vector(QQ,r))
  if len(neg)==1:
   u=cols[neg[0]];g-=u*harmonic(-ks[neg[0]]-1)
   return QQ(0),pref*u,pref*(u.column()*g.row()+g.column()*u.row())
  h=-harmonic(k0,2)*(qt.column()*qt.row())+sum((harmonic(k,2)*(c.column()*c.row()) for c,k in zip(cols,ks)),zero_matrix(QQ,r,r))
  return pref,pref*g,pref*(g.column()*g.row()+h)
 return gamma
class Box:
 def __init__(self,caps):
  self.caps=tuple(caps);self.r=len(caps);self.zero=(0,)*self.r
  self.keys=sorted(product(*(range(c+1) for c in caps)),key=lambda n:(sum(n),n));self.top=sum(caps)
  self.pairs={n:[] for n in self.keys}
  for n in self.keys:
   for m in product(*(range(c-a+1) for c,a in zip(caps,n))):self.pairs[n].append((m,tuple(a+b for a,b in zip(n,m))))
 def clean(self,p):return {n:QQ(v) for n,v in p.items() if v}
 def add(self,p,q,c=1):
  out=dict(p)
  for n,v in q.items():out[n]=out.get(n,0)+c*v
  return self.clean(out)
 def scale(self,p,c):return self.clean({n:c*v for n,v in p.items()})
 def mul(self,p,q):
  if len(p)>len(q):p,q=q,p
  out={}
  for n,v in p.items():
   for m,k in self.pairs[n]:
    w=q.get(m,0)
    if w:out[k]=out.get(k,0)+v*w
  return self.clean(out)
 def exp(self,p):
  assert p.get(self.zero,0)==0
  out={self.zero:QQ(1)};term=dict(out)
  for k in range(1,self.top+1):
   term=self.scale(self.mul(term,p),QQ(1)/k)
   if not term:break
   out=self.add(out,term)
  return out
 def inv(self,p):
  assert p.get(self.zero,0)==1
  h=self.add(p,{self.zero:QQ(1)},-1)
  out={self.zero:QQ(1)};term=dict(out)
  for k in range(1,self.top+1):
   term=self.scale(self.mul(term,h),-1)
   if not term:break
   out=self.add(out,term)
  assert self.mul(p,out)=={self.zero:QQ(1)}
  return out

def compute(Q,kappas,caps,check_jets=False):
 A=Box(caps);r=A.r;gamma=gamma_factory(Q);w0={};cg=[{} for i in range(r)];ch=[{} for i in range(r)]
 for n in A.keys:
  c,g,h=gamma(n)
  if c:w0[n]=c
  for a in range(r):
   if g[a]:cg[a][n]=g[a]
   v=sum(kappas[a][b,k]*h[b,k] for b in range(r) for k in range(r))/2
   if v:ch[a][n]=v
 inv=A.inv(w0);g=[A.mul(v,inv) for v in cg]
 gg={(b,c):A.mul(g[b],g[c]) for b in range(r) for c in range(b,r)}
 magnetic=[]
 for a in range(r):
  p=A.mul(ch[a],inv)
  for b in range(r):
   for c in range(b,r):p=A.add(p,gg[b,c],-kappas[a][b,c]*(QQ(1)/2 if b==c else 1))
  magnetic.append(p)
 # Complete mirror-map substitution in the quotient by (z_i^(caps_i+1)).
 qi=[]
 for i in range(r):
  e=tuple(int(j==i) for j in range(r))
  qi.append(A.mul({e:QQ(1)},A.exp(g[i])) if caps[i] else {})
 qmon={A.zero:{A.zero:QQ(1)}}
 for n in A.keys[1:]:
  i=next(i for i,x in enumerate(n) if x)
  prev=list(n);prev[i]-=1
  qmon[n]=A.mul(qmon[tuple(prev)],qi[i])
  assert qmon[n].get(n,0)==1
 # Triangular inversion of the six master equations; no degree outside the box can contribute.
 remainder=[dict(p) for p in magnetic];gv={};dilog={}
 for n in A.keys[1:]:
  residual=vector(QQ,[p.get(n,0) for p in remainder]);j=next(i for i,x in enumerate(n) if x)
  invariant=residual[j]/n[j]
  assert residual==invariant*vector(QQ,n),('integrability',n,list(residual))
  assert invariant in ZZ,('nonintegral',n,invariant)
  gv[n]=ZZ(invariant)
  if invariant:
   maxk=min(c//x for c,x in zip(caps,n) if x)
   p={}
   for k in range(1,maxk+1):p=A.add(p,qmon[tuple(k*x for x in n)],QQ(1)/k**2)
   dilog[n]=p
   remainder=[A.add(remainder[a],p,-invariant*n[a]) for a in range(r)]
 assert all(not p for p in remainder)
 recurrences=0
 if check_jets:
  def muljet(jet,constant,linear):
   c,g,h=jet
   return constant*c,constant*g+c*linear,constant*h+g.column()*linear.row()+linear.column()*g.row()
  qt=sum(Q.columns())
  for n in A.keys[1:]:
   nn=vector(ZZ,n);ks=nn*Q;k0=nn*qt
   for a in range(r):
    if not n[a]:continue
    prev=list(n);prev[a]-=1;left=gamma(n);right=gamma(tuple(prev))
    for i,qi0 in enumerate(Q.row(a)):
     if qi0>0:
      for m in range(qi0):left=muljet(left,ks[i]-m,Q.column(i))
     elif qi0<0:
      for m in range(1,-qi0+1):right=muljet(right,ks[i]+m,Q.column(i))
    for m in range(1,qt[a]+1):right=muljet(right,k0-qt[a]+m,qt)
    assert left==right,('Gamma recurrence',n,a)
    recurrences+=1
 return A,gv,{'w0':w0,'mirror_log':g,'magnetic':magnetic,'gamma_recurrences':recurrences}
start=time.time()
_,bench,_=compute(matrix(ZZ,[[1]*5]),[matrix(QQ,[[5]])],(3,))
assert [bench[(i,)] for i in [1,2,3]]==[2875,609250,317206375]
print('PASS: independent quintic benchmark.',flush=True)
caps=(2,1,1,0,5,5)
A,gv,periods=compute(Q,kappa,caps,True)
print('PASS:',len(A.keys),'box points;',periods['gamma_recurrences'],'Gamma recurrences; all six integrability equations and integral GV coefficients.',flush=True)
candidates=sorted((a,1,1,0,b,c) for a in range(3) for b in range(max(0,a-1),6) for c in range(max(a,b),min(5,2+b,3+b-a)+1))
assert len(candidates)==37
rows=[{'degree':list(n),'genus_zero_bps':int(gv[n]),'genus_zero_gw':str(gv[n]),'primitive':gcd(n)==1} for n in candidates]
assert all(x['primitive'] for x in rows)
other=[{'degree':list(n),'genus_zero_bps':int(v)} for n,v in gv.items() if n[1:4]==(1,1,0) and n not in candidates and v]
print('CANDIDATE_NONZERO',[(x['degree'],x['genus_zero_bps']) for x in rows if x['genus_zero_bps']],flush=True)
print('CANDIDATE_SUM',sum(x['genus_zero_bps'] for x in rows),flush=True)
print('OUTSIDE_CANDIDATE_NONZERO_IN_BOX',other,flush=True)
def gw(n):
 return sum((QQ(gv[tuple(x//k for x in n)])/k**3 for k in divisors(gcd(n))),QQ(0))
def degree_row(n):
 return {'degree':list(n),'genus_zero_bps':int(gv[n]),'genus_zero_gw':str(gw(n)),'primitive':bool(gcd(n)==1)}
sectors={f'{b},{c}':[degree_row(n) for n in sorted(gv) if n[1:4]==(b,c,0)] for b,c in [(0,0),(1,0),(0,1),(1,1)]}
for key,sector in sectors.items():
 print('SECTOR',key,'nonzero',[(row['degree'],row['genus_zero_bps']) for row in sector if row['genus_zero_bps']],flush=True)
def serial(p):return [{'degree':list(n),'coefficient':str(v)} for n,v in sorted(p.items())]
data={'scope':'genus-zero numerical realization only; not a universal or higher-genus coefficient','source_snapshot_head':'d0d815886953bdc051f0c6c28247aa184a896518','method':'Exact rational Gamma 2-jets, mirror-map substitution, six Li_2 master equations, and multicover extraction in a componentwise downward-closed box. No smoothing invariant is an input.','enumerative_input':'Canonical Higgs toric-hypersurface mirror identification; Demirtas et al. arXiv:2303.00757v2 equation (4.30), section 4.5, and Appendix A.','scope_of_finiteness':'The 37-candidate bound is imported from the separate common-blowup certificate. This script proves only exact arithmetic in the displayed box.','source_sha256':PIN,'caps':list(caps),'box_points':len(A.keys),'max_total_degree':sum(caps),'glsm':[list(map(int,row)) for row in Q],'nef_divisor_matrix':[list(map(int,row)) for row in N],'nef_kappa':[[list(map(str,row)) for row in k] for k in kappa],'gamma_recurrences':periods['gamma_recurrences'],'quintic_benchmark':[2875,609250,317206375],'all_box_bps':[{'degree':list(n),'value':int(v)} for n,v in sorted(gv.items())],'all_box_bps_and_gw':[degree_row(n) for n in sorted(gv)],'bidegree_sectors':sectors,'w0':serial(periods['w0']),'mirror_log':[serial(p) for p in periods['mirror_log']],'magnetic':[serial(p) for p in periods['magnetic']],'candidate_results':rows,'candidate_sum':sum(x['genus_zero_bps'] for x in rows),'nonzero_outside_candidates_with_target_beta':other}
args.output.write_text(json.dumps(data,indent=2)+'\n')
print('OUTPUT',args.output,'seconds',float(time.time()-start),flush=True)
