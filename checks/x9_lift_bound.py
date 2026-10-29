"""Exact outer lift certificate; SageMath 10.9, no period computations.
Run: sage -python x9_lift_bound.py --output certificate.json
All raw inputs are embedded; --source optionally verifies the canonical source hash.
Mathematical interpretation: W is the blowup of the four stated disjoint
(-1,-1) curves on the Delta-regular X9 hypersurface. D10 misses W.
This proves containment in the candidate lists, not effectivity of candidates.
"""
from sage.all import *
from itertools import combinations
from pathlib import Path
import json, hashlib, argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True,help='Exact JSON certificate destination')
parser.add_argument('--source',type=Path,help='Optional canonical check_x9_prepotential.sage; verify pinned SHA256, never execute it')
args=parser.parse_args()
SOURCE_SHA256='c5c25b261f92db49ef4e78d8625e451cac2b01b60c9021b7fa3f96b156ef662b'
if args.source is not None:
 assert hashlib.sha256(args.source.read_bytes()).hexdigest()==SOURCE_SHA256
HEIGHT_BASE=[10,11,3,9,2,7,1,4,8,-5,-3]
HEIGHTS=[QQ(h)+QQ(i+1)/QQ(1009) for i,h in enumerate(HEIGHT_BASE)]
OLD_RAYS=[(1,0,0,0),(0,1,0,0),(-1,-1,0,0),(0,0,1,0),(0,0,0,1),(-4,-2,-1,0),(-4,-2,0,-1),(3,1,1,1),(2,1,1,1),(-2,-1,0,0),(-1,0,0,0)]
OLD_CONES=[(0,1,3,6),(0,1,3,7),(0,1,4,5),(0,1,4,7),(0,1,5,6),(0,2,4,5),(0,2,4,7),(0,2,5,6),(0,2,6,7),(0,3,6,7),(1,3,6,10),(1,3,7,8),(1,3,8,10),(1,4,5,10),(1,4,7,8),(1,4,8,10),(1,5,6,10),(2,4,5,9),(2,4,7,9),(2,5,6,9),(2,6,7,9),(3,6,7,9),(3,6,9,10),(3,7,8,9),(3,8,9,10),(4,5,9,10),(4,7,8,9),(4,8,9,10),(5,6,9,10)]
# Reconstruct the 29-cone source fan from the printed rational heights.
poly=Polyhedron(vertices=OLD_RAYS[:9])
reconstructed=set()
for ineq in poly.inequalities():
 normal=-vector(ineq.A());assert ineq.b()==1
 on=[i for i,v in enumerate(OLD_RAYS) if normal*vector(v)==1]
 if len(on)==4:
  reconstructed.add(tuple(sorted(on)));continue
 drop=next(j for j in range(4) if normal[j]!=0)
 keep=[j for j in range(4) if j!=drop]
 lifted=[tuple([OLD_RAYS[i][j] for j in keep]+[HEIGHTS[i]]) for i in on]
 lift_poly=Polyhedron(vertices=lifted)
 for facet in lift_poly.faces(3):
  constraints=[e for e in facet.ambient_Hrepresentation() if e.is_inequality()]
  assert len(constraints)==1
  if -vector(constraints[0].A())[-1]<0:
   indices=tuple(sorted(on[lifted.index(tuple(v))] for v in facet.vertices()))
   assert len(indices)==4;reconstructed.add(indices)
assert sorted(reconstructed)==sorted(OLD_CONES)
CENTERS=[(2,4),(6,7),(4,9),(3,9)]
BASIS=[2,5,6,7,8,9]
N=matrix(ZZ,[[0,0,1,1,2,2],[1,1,2,4,5,5],[0,1,2,4,4,5],[0,-1,-1,-3,-3,-3],[0,-1,-1,-2,-2,-2],[0,0,1,2,2,2]])
# D_i gamma = column_i(Q).n, columns 0..9; D10 gamma=0.
Q=matrix(ZZ,[[0,0,-1,0,-1,1,0,1,0,0],[0,0,0,1,1,0,0,0,-1,-1],[0,1,0,0,0,0,0,2,-2,1],[1,0,-1,0,0,0,0,0,0,1],[0,0,1,1,0,0,-1,-1,0,0],[0,0,0,-1,0,0,1,0,1,-1]])
assert Q*matrix(ZZ,OLD_RAYS[:10])==0
assert Q.matrix_from_columns(BASIS)*N==identity_matrix(ZZ,6)
assert N.det()==-1
r=[vector(ZZ,x) for x in OLD_RAYS]; cones=list(OLD_CONES)
f=Fan(cones=[list(c) for c in cones],rays=r,check=True)
assert f.is_complete() and f.is_smooth()
for center in CENTERS:
 idx=len(r);r.append(r[center[0]]+r[center[1]])
 assert gcd(r[-1])==1
 hit=[c for c in cones if set(center)<=set(c)];assert hit
 cones=sorted(set([c for c in cones if not set(center)<=set(c)]+[tuple(sorted((set(c)-{omit})|{idx})) for c in hit for omit in center]))
 f=Fan(cones=[list(c) for c in cones],rays=r,check=True)
 assert f.is_complete() and f.is_smooth()
assert len(r)==15 and len(cones)==47
walls={}
for c in cones:
 for wall in combinations(c,3):walls.setdefault(wall,[]).append(next(i for i in c if i not in wall))
wall_rows=[]
for wall,ab in walls.items():
 assert len(ab)==2
 a,b=ab;co=matrix(QQ,[r[i] for i in wall]).T.solve_right(-r[a]-r[b])
 z=zero_vector(QQ,len(r));z[a]=z[b]=1
 for i,t in zip(wall,co):z[i]=t
 assert sum((z[i]*r[i] for i in range(len(r))),zero_vector(QQ,4))==0
 wall_rows.append(z)
assert len(wall_rows)==94
# Existence of a divisor positive on every toric wall verifies projectivity.
ample_space=Polyhedron(ieqs=[[-1]+list(w) for w in wall_rows],backend='ppl')
assert not ample_space.is_empty()
ample=ample_space.representative_point()
assert all(w*ample>=1 for w in wall_rows)
M=Polyhedron(rays=wall_rows,backend='ppl')
eqs=[]
for i in range(10,15):
 z=[0]*16;z[i+1]=1;eqs.append(z)
M0=M.intersection(Polyhedron(eqns=eqs,ambient_dim=15,backend='ppl'))
def deg(v):return N.T*vector(QQ,[v[i] for i in BASIS])
R=Polyhedron(rays=[deg(v.vector()) for v in M0.rays()],backend='ppl')
EXPECTED_RAYS=[(0,0,0,0,1,1),(0,0,0,1,0,0),(0,0,1,0,0,0),(0,0,1,0,0,1),(0,1,0,0,0,0),(0,1,0,0,0,1),(1,0,1,0,1,1),(1,1,0,0,0,1),(1,1,0,0,1,1)]
assert R==Polyhedron(rays=EXPECTED_RAYS,backend='ppl')
# This outer ambient cone alone is insufficient: e5+e6 remains.
assert R.contains(vector(QQ,[0,0,0,0,1,1]))
def pullback(d):
 d=list(map(QQ,d))
 for a,b in CENTERS:d.append(d[a]+d[b])
 return vector(QQ,d)
F=pullback([1]*11)
PF=Polyhedron(ieqs=[[F[i]]+list(r[i]) for i in range(15)],backend='ppl')
fvertices=[v.vector() for v in PF.vertices()]
allc=sorted(set(c for top in cones for k in range(1,5) for c in combinations(top,k)),key=lambda c:(-len(c),c))
meet=[]; empty_faces=[]
for c in allc:
 face=[v for v in fvertices if all(F[i]+r[i]*v==0 for i in c)]
 assert face
 if len(face)==1:empty_faces.append((c,face[0]))
 elif not any(set(c)<=set(d) for d in meet):meet.append(c)
assert len(meet)==41 and all(10 not in c for c in meet)
# D=phi*(2H2+3H3+H4-H6)-sum of the four exceptional divisors.
bcoef=N*vector(ZZ,[0,2,3,1,0,-1]);old=zero_vector(ZZ,11)
for i,x in zip(BASIS,bcoef):old[i]=x
D=pullback(old)
for i in range(11,15):D[i]-=1
assert list(D)==[0,0,2,0,0,7,7,-5,-5,3,0,1,1,2,2]
# All sections are global on the SAME open U obtained by removing D10.
PD=Polyhedron(ieqs=[[D[i]]+list(r[i]) for i in range(15) if i!=10],backend='ppl')
dvertices=[v.vector() for v in PD.vertices()]
witness=[]
for c in meet:
 v=next(v for v in dvertices if all(D[i]+r[i]*v==0 for i in c))
 assert all(D[i]+r[i]*v>=0 for i in range(15) if i!=10)
 witness.append((c,v))
den=lcm([x.denominator() for c,v in witness for x in v])
assert den==2
# Genuine integral character sections of O_U(2D).
for c,v in witness:
 m=vector(ZZ,den*v)
 assert all(den*D[i]+r[i]*m>=0 for i in range(15) if i!=10)
 assert all(den*D[i]+r[i]*m==0 for i in c)
# Thus D|W is nef; on jgamma this gives n6<=2n2+3n3+n4.
result={}
for d,expected in [(1,37),(2,183)]:
 eq=[]
 for i,k in [(1,d),(2,d),(3,0)]:
  z=[-k]+[0]*6;z[i+1]=1;eq.append(z)
 P=R.intersection(Polyhedron(eqns=eq,ieqs=[[0,0,2,3,1,0,-1]],ambient_dim=6,backend='ppl'))
 assert P.is_compact()
 cand=sorted(tuple(map(int,v)) for v in P.integral_points())
 elementary=sorted((a,d,d,0,b,c) for a in range(2*d+1) for b in range(max(0,a-d),5*d+1) for c in range(max(a,b),min(5*d,2*d+b,3*d+b-a)+1))
 assert cand==elementary and len(cand)==expected
 result[str(d)]={'count':len(cand),'vertices':[[str(x) for x in v.vector()] for v in P.vertices()],'candidates':cand}
def vec(v):return [str(x) for x in v]
data={'source_sha256':SOURCE_SHA256,'height_base':HEIGHT_BASE,'height_denominator':1009,'rational_heights':[str(x) for x in HEIGHTS],'HEAD':'d0d815886953bdc051f0c6c28247aa184a896518','old_rays':OLD_RAYS,'old_cones':OLD_CONES,'centers':CENTERS,'rays':[list(map(int,v)) for v in r],'cones':cones,'wall_relations':[vec(v) for v in wall_rows],'ample_witness':vec(ample),'source_nef_matrix':[list(map(int,row)) for row in N],'source_curve_divisor_matrix':[list(map(int,row)) for row in Q],'ambient_outer_rays':[vec(v.vector()) for v in R.rays()],'ambient_outer_inequalities':[vec(v.vector()) for v in R.inequalities()],'anticanonical_pullback':vec(F),'anticanonical_polytope_vertices':[vec(v) for v in fvertices],'maximal_possibly_met_cones':meet,'empty_orbit_vertex_faces':[{'cone':c,'vertex':vec(v)} for c,v in empty_faces],'nef_divisor':vec(D),'monomial_denominator':int(den),'global_monomial_witnesses':[{'cone':c,'character_2D':list(map(int,den*v)),'slack_2D':[str(den*D[i]+r[i]*(den*v)) for i in range(15)]} for c,v in witness],'lift_candidates':result}
out=args.output
out.write_text(json.dumps(data,indent=2)+'\n')
print('PASS: smooth projective 15-ray,47-cone ambient common blowup;94 wall relations.')
print('PASS: global integral sections of 2D generate on all41 possibly met maximal cones.')
print('PASS: complete outer candidate sets:37 over beta;183 over2beta.')
print('These are necessary-condition candidate sets, not asserted effective lift sets.')
print('JSON SHA256',hashlib.sha256(out.read_bytes()).hexdigest())
