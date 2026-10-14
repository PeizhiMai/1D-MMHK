#!/usr/bin/env python3
"""Independent dense Fock-space ED checks (does not reuse qbasis operators)."""
import importlib.util
import itertools
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location("post",ROOT/"response/postprocess_response.py")
post=importlib.util.module_from_spec(spec);spec.loader.exec_module(post)


def mask_states(L,N):
    return [sum(1<<i for i in inds) for inds in itertools.combinations(range(L),N)]


def independent_ed(L,Nup,Ndn,U,phase,omega,eta,observable="density"):
    states=[u|(d<<L) for u in mask_states(L,Nup) for d in mask_states(L,Ndn)]
    where={s:i for i,s in enumerate(states)};dim=len(states)
    H=np.zeros((dim,dim),dtype=complex)
    for col,s in enumerate(states):
        H[col,col]=U*sum(((s>>j)&1)*((s>>(j+L))&1) for j in range(L))
        for spin in [0,1]:
            for j in range(L):
                a,b=j+spin*L,(j+1)%L+spin*L
                for dest,src,coef in [(a,b,-np.exp(1j*phase)),(b,a,-np.exp(-1j*phase))]:
                    if not (s>>src)&1 or (s>>dest)&1:continue
                    sign=(-1)**((s&((1<<src)-1)).bit_count())
                    v=s^(1<<src);sign*=(-1)**((v&((1<<dest)-1)).bit_count())
                    v|=1<<dest;H[where[v],col]+=coef*sign
    np.testing.assert_allclose(H,H.conj().T,atol=1e-14)
    energies,V=np.linalg.eigh(H);gs=np.flatnonzero(energies-energies[0]<1e-9)
    spectra={};weights={};elastic={}
    for spin,channel in [(False,"charge"),(True,"spin")]:
        spectra[channel]=np.zeros((L,len(omega)));weights[channel]=np.zeros(L);elastic[channel]=np.zeros(L)
        for q in range(L):
            sources=[]
            if observable=="density":
                density=np.array([sum(np.exp(2j*np.pi*q*j/L)*(((s>>j)&1)+(-1 if spin else 1)*((s>>(j+L))&1)) for j in range(L))/np.sqrt(L) for s in states])
                sources=[density[:,None]*V[:,gs]]
            else:
                for K in range(L):
                    A=np.zeros_like(H)
                    for col,s in enumerate(states):
                        for sigma in [0,1]:
                            for i in range(L):
                                for j in range(L):
                                    dest,src=i+sigma*L,j+sigma*L
                                    coef=np.exp(2j*np.pi*(-K*i+(K+q)*j)/L)/(L*np.sqrt(L))
                                    if spin and sigma:coef=-coef
                                    if dest==src:
                                        A[col,col]+=coef*((s>>src)&1);continue
                                    if not(s>>src)&1 or(s>>dest)&1:continue
                                    sign=(-1)**((s&((1<<src)-1)).bit_count());v=s^(1<<src)
                                    sign*=(-1)**((v&((1<<dest)-1)).bit_count());v|=1<<dest
                                    A[where[v],col]+=coef*sign
                    sources.append(A@V[:,gs])
            for source_matrix in sources:
                for ig,g in enumerate(gs):
                    W=abs(V.conj().T@source_matrix[:,ig])**2/len(gs)
                    E=energies-energies[g];positive=E>1e-8
                    spectra[channel][q]+=post.lorentzian(omega,E[positive],W[positive],eta)
                    weights[channel][q]+=W.sum();elastic[channel][q]+=W[abs(E)<=1e-8].sum()
    return {"E0":float(energies[0]),"ground_degeneracy":len(gs),"spectra":spectra,"weights":weights,"elastic":elastic}


def main():
    exe=ROOT/"response/bin/chain_Hubbard_tbc_response.x"
    env=dict(os.environ,OMP_NUM_THREADS="2",MKL_NUM_THREADS="1")
    omega=np.linspace(-1,16,1701);eta=0.2;reports=[]
    tests=[("half_U8",2,2,8.,0.,"density"),("doped_U8",2,1,8.,0.,"density"),
           ("dilute_U0",1,1,0.,0.,"density"),("twisted_U4",2,1,4.,0.137,"density"),
           ("pairs_half_U8",2,2,8.,0.,"orbital-pairs"),("pairs_twisted_U4",2,1,4.,.137,"orbital-pairs")]
    base=ROOT/"response/results/validation";base.mkdir(parents=True,exist_ok=True)
    for name,nu,nd,U,phase,observable in tests:
        with tempfile.TemporaryDirectory(prefix=name+"_",dir=base) as temp:
            run=Path(temp)
            with (run/"run.log").open("w") as log:
                subprocess.run([str(exe),"--L","4","--nup",str(nu),"--ndn",str(nd),
                                "--U",str(U),"--phase",str(phase),"--depth","120","--observable",observable,"--output",str(run)],
                               check=True,env=env,stdout=log,stderr=subprocess.STDOUT)
            actual=post.reconstruct(run,omega,eta);exact=independent_ed(4,nu,nd,U,phase,omega,eta,observable)
            np.testing.assert_allclose(actual["metadata"]["E0"],exact["E0"],atol=2e-9,rtol=0)
            errors={}
            for channel in ["charge","spin"]:
                np.testing.assert_allclose(actual[channel],exact["spectra"][channel],atol=2e-8,rtol=2e-8)
                np.testing.assert_allclose(actual[f"weight_{channel}"],exact["weights"][channel],atol=1e-8,rtol=1e-8)
                np.testing.assert_allclose(actual[f"elastic_{channel}"],exact["elastic"][channel],atol=1e-8,rtol=1e-8)
                errors[channel]=float(np.max(abs(actual[channel]-exact["spectra"][channel])))
            reports.append({"name":name,"passed":True,"E0":exact["E0"],
                            "ground_degeneracy":int(exact["ground_degeneracy"]),"max_spectral_errors":errors})
            print("PASS",name,reports[-1])
    for args in [["--filling","0.453125"],["--filling","0.51"],["--nup","-2","--ndn","-2"],["--q","10"]]:
        completed=subprocess.run([str(exe),*args],capture_output=True,env=env,text=True)
        assert completed.returncode != 0, args
    (base/"test_report.json").write_text(json.dumps(reports,indent=2)+"\n")
    print("PASS invalid-input rejection; all independent checks passed")


if __name__=="__main__":main()
