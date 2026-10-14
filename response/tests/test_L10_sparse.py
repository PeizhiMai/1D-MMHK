#!/usr/bin/env python3
"""Independent unreduced 52920-state ED/Krylov check of the actual L=10 run."""
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import postprocess_response as post


def main():
    start=time.monotonic();L=10;U=8.;Nup=5;Ndn=4
    masks=lambda N:[sum(1<<j for j in inds) for inds in itertools.combinations(range(L),N)]
    states=[u|(d<<L) for u in masks(Nup) for d in masks(Ndn)]
    where={s:i for i,s in enumerate(states)};dim=len(states)
    row=[];col=[];data=[]
    for k,s in enumerate(states):
        row.append(k);col.append(k);data.append(U*sum(((s>>j)&1)*((s>>(j+L))&1) for j in range(L)))
        for sigma in [0,1]:
            for j in range(L):
                a,b=j+L*sigma,(j+1)%L+L*sigma
                for dest,src in [(a,b),(b,a)]:
                    if not(s>>src)&1 or(s>>dest)&1:continue
                    sign=(-1)**((s&((1<<src)-1)).bit_count());v=s^(1<<src)
                    sign*=(-1)**((v&((1<<dest)-1)).bit_count());v|=1<<dest
                    row.append(where[v]);col.append(k);data.append(-sign)
    H=coo_matrix((data,(row,col)),shape=(dim,dim)).tocsr()
    anti=H-H.T;assert anti.nnz==0 or np.max(abs(anti.data))<1e-14
    print("Independent matrix",dim,"states, nnz",H.nnz,flush=True)
    energies,V=eigsh(H,k=4,ncv=40,which="SA",tol=1e-11,maxiter=4000,v0=np.random.default_rng(42).normal(size=dim))
    order=np.argsort(energies);energies=energies[order];V=V[:,order]
    gs=np.flatnonzero(energies-energies[0]<1e-8)
    assert len(gs)==2,(energies,gs)
    print("Independent eigenvalues",energies,flush=True)
    omega=np.linspace(0,5,5001);eta=.2
    run=HERE/"results/zero_twist_benchmark/L10_U8_n0p45_phi0"
    actual=post.reconstruct(run,omega,eta)
    np.testing.assert_allclose(energies[0],actual["metadata"]["E0"],atol=1e-9,rtol=0)
    report={"dimension":dim,"nnz":int(H.nnz),"independent_lowest_eigenvalues":energies.tolist(),
            "ground_degeneracy":len(gs),"channels":{}}
    for spin,channel in [(False,"charge"),(True,"spin")]:
        density=np.array([sum(np.exp(2j*np.pi*j/L)*(((s>>j)&1)+(-1 if spin else 1)*((s>>(j+L))&1)) for j in range(L))/np.sqrt(L) for s in states])
        spectrum=np.zeros_like(omega)
        for g in gs:
            source=density*V[:,g];weight=np.vdot(source,source).real
            v=source/np.sqrt(weight);previous=np.zeros_like(v);beta=0.;a=[];b=[]
            # An ARPACK ground vector can mix both degenerate Q sectors, so
            # this unreduced source may need a deeper chain than each separate
            # translation-sector source used by the production calculation.
            for it in range(600):
                residual=H@v-beta*previous
                alpha=np.vdot(v,residual).real;residual-=alpha*v
                next_beta=np.linalg.norm(residual);a.append(alpha)
                if it:b.append(beta)
                if next_beta<1e-12:break
                previous,v=v,residual/next_beta;beta=next_beta
            from scipy.linalg import eigh_tridiagonal
            ep,z=eigh_tridiagonal(np.asarray(a),np.asarray(b),lapack_driver="stev")
            spectrum+=post.lorentzian(omega,ep-energies[g],weight*z[0]**2/len(gs),eta)
        np.testing.assert_allclose(spectrum,actual[channel][1],atol=1e-7,rtol=1e-6)
        diff=float(np.max(abs(spectrum-actual[channel][1])))
        report["channels"][channel]={"max_abs_error":diff,"relative_to_max":diff/float(np.max(spectrum)),
                                    "independent_first_peak":post.first_peak(omega,spectrum),
                                    "qbasis_first_peak":post.first_peak(omega,actual[channel][1])}
        print("PASS L10",channel,report["channels"][channel],flush=True)
    report["elapsed_seconds"]=time.monotonic()-start;report["passed"]=True
    output=HERE/"results/validation/L10_independent_sparse_report.json"
    output.write_text(json.dumps(report,indent=2)+"\n")
    print("PASS independent L10 sparse ED",output)


if __name__=="__main__":main()
