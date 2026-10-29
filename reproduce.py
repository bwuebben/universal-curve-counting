#!/usr/bin/env python3
"""Replay the exact finite calculations accompanying the two papers."""
import argparse, contextlib, hashlib, io, json, os, runpy, shutil, subprocess, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--sage',default='sage');args=parser.parse_args()
 import sympy as S
 assert S.__version__=='1.14.0',S.__version__
 os.chdir(ROOT)
 for line in (ROOT/'SHA256SUMS').read_text().splitlines():
  digest,path=line.split('  ',1);assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
 with tempfile.TemporaryDirectory(prefix='curve-counting-') as temp:
  tmp=Path(temp)
  def run(name, certificate=None, sage=False, stdout_json=False, extra=()):
   command=[args.sage,'-python'] if sage else [sys.executable]
   command+=['checks/'+name+'.py',*extra]
   output=tmp/(certificate or (name+'.txt'))
   if certificate and not stdout_json: command+=['--output',str(output)]
   print('Running',name,flush=True)
   p=subprocess.run(command,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   if p.returncode:raise RuntimeError(p.stdout+'\n'+p.stderr)
   if stdout_json:output.write_text(p.stdout)
   if certificate:assert json.loads(output.read_text())==json.loads((ROOT/'evidence'/certificate).read_text()),certificate
  run('x9_plane_geometry','x9-plane-fan-certificate.json',stdout_json=True)
  run('x9_source_chamber','x9-source-chamber-certificate.json',stdout_json=True)
  run('x9_lift_bound','x9-lift-bound-certificate.json',sage=True)
  run('x9_source_genus_zero_box','x9-source-genus-zero-box.json',sage=True,extra=('--source-root','vendor/global-higgs-branch'))
  run('x9_contracted_support','x9-contracted-support-certificate.json',sage=True)
  run('integrality_counterexample_geometry','integrality-counterexample-geometry.json',stdout_json=True)
  run('x9_universal_line','x9-universal-line-certificate.json')
  run('local_series_comparison')
  run('local_curve_degree_three')
  run('integrality_counterexample','integrality-counterexample-certificate.json')
  sys.path.insert(0,str(ROOT/'checks'))
  from local_curve_degree_three import build,in_w,T
  with contextlib.redirect_stdout(io.StringIO()):quotient=runpy.run_path(str(ROOT/'checks/degree_three_quotient_check.py'))
  B=build(3)[-1]
  coefficients=[in_w(x) for x in B.charpoly(T).all_coeffs()]
  assert all(S.simplify(a-b)==0 for a,b in zip(coefficients,quotient['co']))
  assert S.simplify(in_w(S.trace(B**2))-quotient['tr2'])==0
  print('PASS: matrix and independent quotient characteristic polynomials and traces agree.',flush=True)
  # Reconstruct the appendix in a temporary tree, leaving the manuscript intact.
  shutil.copytree(ROOT/'evidence',tmp/'evidence');(tmp/'paper/sections').mkdir(parents=True)
  shutil.copy2(ROOT/'paper/generate_data.py',tmp/'paper/generate_data.py')
  subprocess.run([sys.executable,str(tmp/'paper/generate_data.py')],cwd=tmp,check=True,stdout=subprocess.PIPE)
  assert (tmp/'paper/sections/data.tex').read_bytes()==(ROOT/'paper/sections/data.tex').read_bytes()
 print('PASS: all eight certificates, local identities, independent algebra, and printed tables agree.')
if __name__=='__main__':main()
