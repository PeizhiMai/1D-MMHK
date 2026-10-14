#!/usr/bin/env python3
"""Optional independent 44,100-state sparse ED/Krylov portable-backend anchor.

Run in the full Docker image. Writes only the requested new output directory.
This checks a canonical Sz=0 sector, not a substitute for the grand-ground
trace used by the full density workflow.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'response'))
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from reproducibility.workflow import solvers
from postprocess_response import reconstruct
from run_single_particle_scan import fraction_G

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve()
if out.exists():raise RuntimeError('Choose a new independent-check output directory')
out.mkdir(parents=True);start=time.monotonic()
L=10;U=8.;nu=nd=4
masks=lambda N:[sum(1<<j for j in v) for v in itertools.combinations(range(L),N)]
states=[u|(d<<L) for u in masks(nu) for d in masks(nd)]
where={s:i for i,s in enumerate(states)};rows=[];cols=[];values=[]
for col,s in enumerate(states):
    rows.append(col);cols.append(col);values.append(U*sum(((s>>j)&1)*((s>>(j+L))&1) for j in range(L)))
    for spin in [0,1]:
        for j in range(L):
            x,y=j+spin*L,(j+1)%L+spin*L
            for dest,src in [(x,y),(y,x)]:
                if not(s>>src)&1 or (s>>dest)&1:continue
                sign=(-1)**((s&((1<<src)-1)).bit_count());v=s^(1<<src)
                sign*=(-1)**((v&((1<<dest)-1)).bit_count());v|=1<<dest
                rows.append(where[v]);cols.append(col);values.append(-sign)
H=coo_matrix((values,(rows,cols)),shape=(len(states),len(states))).tocsr()
assert (H-H.T).nnz==0
E,V=eigsh(H,k=4,ncv=40,which='SA',tol=1e-11,maxiter=4000,v0=np.random.default_rng(42).normal(size=len(states)))
idx=np.argsort(E);E=E[idx];V=V[:,idx]
assert E[1]-E[0]>1e-8;g=V[:,0];residual=float(np.linalg.norm(H@g-E[0]*g))
assert residual<1e-8
run=out/'native'
with (out/'native.log').open('w') as log:
    subprocess.run([str(solvers()['response']),'--L','10','--U','8','--nup','4','--ndn','4','--q','1','--phase','0','--depth','300','--output',str(run)],check=True,
                   env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),stdout=log,stderr=subprocess.STDOUT)
omega=np.linspace(0,5,5001);eta=.1;actual=reconstruct(run,omega,eta)
np.testing.assert_allclose(actual['metadata']['E0'],E[0],atol=1e-8,rtol=0)
errors={};source_checks={}
for spin,ch in [(False,'charge'),(True,'spin')]:
    density=np.array([sum(np.exp(2j*np.pi*j/L)*(((s>>j)&1)+(-1 if spin else 1)*((s>>(j+L))&1)) for j in range(L))/np.sqrt(L) for s in states])
    source=density*g;weight=float(np.vdot(source,source).real);v=source/np.sqrt(weight)
    old=np.zeros_like(v);beta=0.;alpha=[];b=[0.]
    for i in range(600):
        r=H@v-beta*old;ai=float(np.vdot(v,r).real);r-=ai*v
        bn=float(np.linalg.norm(r));alpha.append(ai)
        if bn<1e-12:break
        old,v=v,r/bn;beta=bn
        if i<599:b.append(bn)
    cf=dict(steps=len(alpha),norm=np.sqrt(weight),a=np.asarray(alpha),b=np.asarray(b),E0=float(E[0]))
    spectrum=-fraction_G(cf,omega+1j*eta,True).imag/np.pi
    error=float(abs(spectrum-actual[ch][1]).max());errors[ch]=error
    np.testing.assert_allclose(spectrum,actual[ch][1],atol=1e-7,rtol=1e-6)
    source_checks[ch]=dict(weight=weight,first_moment=float(np.vdot(source,(H@source-E[0]*source)).real))
report=dict(passed=True,L=L,U=U,N=nu+nd,Nup=nu,Ndn=nd,phase_per_bond=0,dimension=len(states),nnz=H.nnz,
            eigenvalues=E.tolist(),independent_ground_residual=residual,native_E0=actual['metadata']['E0'],
            maximum_spectral_errors=errors,source_checks=source_checks,eta=eta,omega_range=[0,5],
            elapsed_seconds=time.monotonic()-start,scope='canonical sector anchor, independent unreduced sparse ED/Krylov; no historical arrays')
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
