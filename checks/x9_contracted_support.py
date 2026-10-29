"""Exact toric support data for the X9 source Euler-characteristic bound.
Run with SageMath: sage -python checks/x9_contracted_support.py --output result.json
This certifies fan and Newton-face data; Paper 1 supplies the scheme-theoretic proof.
"""
from sage.all import *
from itertools import combinations
from pathlib import Path
import json, argparse, ast, hashlib
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output",type=Path,required=True)
args=parser.parse_args()
INPUT_SHA256 = '5e74c32b0b40ae4fb8d0d1439455bdc943f51c398bb164aacd776ea66d11fa18'
RAYS = [(1, 0, 0, 0), (0, 1, 0, 0), (-1, -1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-4, -2, -1, 0), (-4, -2, 0, -1), (3, 1, 1, 1), (2, 1, 1, 1), (-2, -1, 0, 0), (-1, 0, 0, 0)]
CONES = [(0, 1, 3, 6), (0, 1, 3, 7), (0, 1, 4, 5), (0, 1, 4, 7), (0, 1, 5, 6), (0, 2, 4, 5), (0, 2, 4, 7), (0, 2, 5, 6), (0, 2, 6, 7), (0, 3, 6, 7), (1, 3, 6, 10), (1, 3, 7, 8), (1, 3, 8, 10), (1, 4, 5, 10), (1, 4, 7, 8), (1, 4, 8, 10), (1, 5, 6, 10), (2, 4, 5, 9), (2, 4, 7, 9), (2, 5, 6, 9), (2, 6, 7, 9), (3, 6, 7, 9), (3, 6, 9, 10), (3, 7, 8, 9), (3, 8, 9, 10), (4, 5, 9, 10), (4, 7, 8, 9), (4, 8, 9, 10), (5, 6, 9, 10)]
# Check the embedded data against the pinned source.
source=Path(__file__).with_name('x9_lift_bound.py')
assert hashlib.sha256(source.read_bytes()).hexdigest()==INPUT_SHA256
assignments={}
for node in ast.parse(source.read_text()).body:
 if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ('OLD_RAYS','OLD_CONES'):
  assignments[node.targets[0].id]=ast.literal_eval(node.value)
assert RAYS==assignments['OLD_RAYS'] and CONES==assignments['OLD_CONES']
r=[vector(ZZ,v) for v in RAYS]
fan=Fan(cones=[list(c) for c in CONES],rays=r,check=True)
assert fan.is_complete() and fan.is_smooth()
P=Polyhedron(ieqs=[[int(i==0)]+list(v) for i,v in enumerate(r)])
A=Polyhedron(ieqs=[[1]+list(v) for v in r])
verts=[v.vector() for v in P.vertices()];averts=[v.vector() for v in A.vertices()]
assert P.dim()==4 and all(x in ZZ for v in verts for x in v)
# The support function of D0 is convex on this fan, hence globally generated.
for c in CONES:
 m=matrix(QQ,[r[i] for i in c]).solve_right(vector(QQ,[-int(i==0) for i in c]))
 assert all(int(i==0)+r[i]*m>=0 for i in range(11))
allc=sorted({tuple(c) for top in CONES for k in range(1,5) for c in combinations(top,k)},key=lambda c:(len(c),c))
exc=[];records=[]
for c in allc:
 face=[v for v in verts if all(int(i==0)+r[i]*v==0 for i in c)]
 assert face
 image_cone_dim=4-Polyhedron(vertices=face).dim()
 if image_cone_dim>len(c):
  records.append((c,image_cone_dim))
  if not any(set(e)<=set(c) for e in exc):exc.append(c)
assert exc==[(8,),(9,),(10,),(2,4),(6,7)]
def antiface(c):return [v for v in averts if all(1+r[i]*v==0 for i in c)]
assert len(antiface((10,)))==1
assert len(antiface((6,7,9)))==len(antiface((2,4,9)))==1
assert not any({6,7,8}<=set(c) or {2,4,8}<=set(c) or {2,4,6,7}<=set(c) for c in CONES)
stars=sorted(tuple(i for i in c if i!=8) for c in CONES if 8 in c and 10 not in c)
assert stars==[(1,3,7),(1,4,7),(3,7,9),(4,7,9)]
project=lambda v:(v[1]-v[3],v[2]-v[3])
base={i:project(r[i]) for i in (1,3,4,9)}
assert base=={1:(1,0),3:(0,1),4:(-1,-1),9:(-1,0)}
face=antiface((8,));fiber_vertex=vector(ZZ,(1,-1,-1,-1))
assert fiber_vertex in face
assert [1+r[i]*fiber_vertex for i in (1,3,4,7,9,10)]==[0,0,0,1,0,0]
assert all(1+r[7]*m==0 for m in face if m!=fiber_vertex)
assert antiface((9,))==[vector(ZZ,(0,1,-1,-1)),fiber_vertex] or set(map(tuple,antiface((9,))))=={(0,1,-1,-1),(1,-1,-1,-1)}
# Positive-fiber loci support the ambient exceptional set. The statement
# about the hypersurface additionally uses Delta-regularity, its nonzero
# vertex coefficients, and the graph argument in the accompanying proof.
out={'input_sha256':INPUT_SHA256,'D0_polytope_vertices':[list(v) for v in verts],
 'minimal_contracted_orbit_cones':exc,'all_contracted_orbits':records,
 'D8_star_away_D10':stars,'D8_base_rays':base,
 'D8_anticanonical_vertices':[list(v) for v in face],
 'D8_fiber_vertex':list(fiber_vertex),
 'E_anticanonical_edge':[list(v) for v in antiface((9,))],
 'empty_intersection_faces':{str(c):[list(v) for v in antiface(c)] for c in [(10,),(6,7,9),(2,4,9)]}}
args.output.write_text(json.dumps(out,default=str,indent=2)+'\n')
print('PASS: contracted-locus support, disjointness, D8 graph fan and E primitive edge.')
