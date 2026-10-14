#!/usr/bin/env python3
"""Independent L<=6 Hubbard-TBC ED benchmark of the physical MMHK response.

No student modules are imported. The reference builds the boundary-bond gauge
directly and diagonalizes fixed (Nup,Ndn) Fock sectors without translation
reduction. It compares all neutral momenta and the complete spin ground trace.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import eigh, eigvalsh

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from mmhk_density_workflow import DEGENERACY_TOL, finite_ring_probe_weights
from postprocess_response import lorentzian


@lru_cache(None)
def sector_basis(L,Nu,Nd):
    masks=lambda n:[sum(1<<j for j in js) for js in itertools.combinations(range(L),n)]
    states=[u|(d<<L) for u in masks(Nu) for d in masks(Nd)]
    return states,{s:i for i,s in enumerate(states)}


def boundary_hamiltonian(L,Nu,Nd,U,theta):
    states,where=sector_basis(L,Nu,Nd);H=np.zeros((len(states),len(states)),complex)
    for column,state in enumerate(states):
        H[column,column]=U*sum(((state>>j)&1)*((state>>(L+j))&1) for j in range(L))
        for spin in [0,1]:
            bonds=[(j+1+spin*L,j+spin*L,-1.) for j in range(L-1)]
            bonds.append((spin*L,L-1+spin*L,-np.exp(1j*theta)))
            for a,b,hopping in bonds:
                for dest,src,coef in [(a,b,hopping),(b,a,np.conj(hopping))]:
                    if not(state>>src)&1 or(state>>dest)&1:continue
                    sign=(-1)**((state&((1<<src)-1)).bit_count());new=state^(1<<src)
                    sign*=(-1)**((new&((1<<dest)-1)).bit_count());new|=1<<dest
                    H[where[new],column]+=coef*sign
    np.testing.assert_allclose(H,H.conj().T,atol=1e-14)
    return H


def canonical_response(L,Nu,Nd,E,V,omega,eta,q_values):
    gs=np.flatnonzero(E-E[0]<DEGENERACY_TOL);states,_=sector_basis(L,Nu,Nd)
    result={ch:np.zeros((len(q_values),len(omega))) for ch in ['charge','spin']}
    for ch,sign in [('charge',1),('spin',-1)]:
        for iq,q in enumerate(q_values):
            diagonal=np.array([sum(np.exp(1j*q*j)*(((s>>j)&1)+sign*((s>>(L+j))&1)) for j in range(L))/np.sqrt(L) for s in states])
            amplitudes=V.conj().T@(diagonal[:,None]*V[:,gs])
            for column,g in enumerate(gs):
                weights=abs(amplitudes[:,column])**2/len(gs)
                gaps=E-E[g];positive=gaps>1e-8
                result[ch][iq]+=lorentzian(omega,gaps[positive],weights[positive],eta)
    return result,len(gs)


def benchmark(directory):
    start=time.monotonic();directory=Path(directory)
    summary=json.loads((directory/'summary.json').read_text())
    L=summary['L'];U=summary['U'];mu=summary['mu_selection']['mu']
    if L>6:raise ValueError('Independent reference execution is limited to at most 6 orbitals')
    with np.load(directory/'physical_density_response.npz',allow_pickle=False) as saved:
        saved={key:saved[key].copy() for key in saved.files}
    with np.load(directory/'ground_energy_table.npz',allow_pickle=False) as table:
        actual_table=table['E0_by_N'].copy()
    k=saved['k'];wk=saved['k_weights'];omega=saved['omega'];eta=float(saved['eta'])
    off_q=np.pi/5
    q_values=np.r_[saved['q'],off_q]
    expected={ch:np.zeros((len(k),len(q_values),len(omega))) for ch in ['charge','spin']}
    table=np.zeros_like(actual_table);means=[];degeneracies=[];max_ring_errors={'charge':0.,'spin':0.}
    for ik,kval in enumerate(k):
        theta=L*kval
        for N in range(2*L+1):
            table[ik,N]=eigvalsh(boundary_hamiltonian(L,(N+1)//2,N//2,U,theta),check_finite=False)[0]
        grand=table[ik]-mu*np.arange(2*L+1)
        numbers=np.flatnonzero(grand-grand.min()<DEGENERACY_TOL)
        components=[]
        for N in numbers:
            # Independently examine ALL nonnegative-Sz sectors, not the
            # production workflow's monotonic early-termination shortcut.
            for Nu in range((int(N)+1)//2,min(L,int(N))+1):
                Nd=int(N)-Nu
                if not 0<=Nd<=L:continue
                E,V=eigh(boundary_hamiltonian(L,Nu,Nd,U,theta),check_finite=False)
                if abs(E[0]-table[ik,N])>DEGENERACY_TOL:continue
                response,deg=canonical_response(L,Nu,Nd,E,V,omega,eta,q_values)
                components.append((int(N),1 if Nu==Nd else 2,deg,response))
        dimension=sum(mult*deg for _,mult,deg,_ in components)
        means.append(sum(N*mult*deg for N,mult,deg,_ in components)/dimension)
        degeneracies.append(dimension)
        for N,mult,deg,response in components:
            probability=mult*deg/dimension
            for ch in expected:expected[ch][ik]+=probability*response[ch]
        np.testing.assert_allclose(means[-1],saved['mean_N_per_k'][ik],atol=1e-10)
        assert dimension==summary['k_blocks'][ik]['grand_ground_dimension']
        for ch in expected:
            error=float(np.max(abs(expected[ch][ik,:L]-saved[f'S_{ch}_per_k'][ik])))
            max_ring_errors[ch]=max(max_ring_errors[ch],error)
            np.testing.assert_allclose(expected[ch][ik,:L],saved[f'S_{ch}_per_k'][ik],atol=1e-7,rtol=1e-6)
            predicted=finite_ring_probe_weights(off_q,L)@saved[f'S_{ch}_per_k'][ik]
            np.testing.assert_allclose(expected[ch][ik,L],predicted,atol=1e-7,rtol=1e-6)
        print(f'PASS independent boundary-gauge ring U={U:g} k={ik+1}/{len(k)}, N={means[-1]:g}',flush=True)
    energy_error=float(np.max(abs(table-actual_table)))
    np.testing.assert_allclose(table,actual_table,atol=2e-9,rtol=0)
    average_errors={};probe_errors={}
    for ch in expected:
        average=np.einsum('k,kqw->qw',wk,expected[ch])
        np.testing.assert_allclose(average[:L],saved[f'S_{ch}'],atol=1e-7,rtol=1e-6)
        average_errors[ch]=float(np.max(abs(average[:L]-saved[f'S_{ch}'])))
        probe_errors[ch]=float(np.max(abs(average[L]-finite_ring_probe_weights(off_q,L)@saved[f'S_{ch}'])))
    report={'passed':True,'L':L,'U':U,'mu':mu,'reference':'independent full-sector Hubbard ED in boundary-bond TBC gauge',
            'k_count':len(k),'total_fixed_N_ground_entries_compared':int(table.size),
            'max_ground_energy_error':energy_error,'max_per_ring_spectral_errors':max_ring_errors,
            'max_averaged_spectral_errors':average_errors,'actual_filling':float(wk@np.asarray(means)/(2*L)),
            'off_grid_finite_ring_probe_q':float(off_q),'off_grid_finite_ring_probe_errors':probe_errors,
            'off_grid_scope':'Finite-ring window probe ONLY; not generic-q physical MMHK inter-block susceptibility',
            'elapsed_seconds':time.monotonic()-start}
    (directory/'independent_hubbard_tbc_benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directories',nargs='+',type=Path);args=parser.parse_args()
    for directory in args.directories:benchmark(directory)
