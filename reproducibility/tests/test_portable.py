"""Independent physical checks; no historical outputs/binaries are needed."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'response'))
from reproducibility.workflow import momentum_plan, load_config, solvers
from postprocess_response import reconstruct, read_fraction
from run_single_particle_scan import fraction_G
spec=importlib.util.spec_from_file_location('independent',ROOT/'response/tests/test_response.py')
independent=importlib.util.module_from_spec(spec);spec.loader.exec_module(independent)

@pytest.fixture(scope='module')
def exe():
    try:return solvers()
    except FileNotFoundError:pytest.skip('Native checks require freshly built portable solvers')

def run_native(exe,p,args):
    p.mkdir(parents=True,exist_ok=True)
    with (p/'run.log').open('w') as f:
        subprocess.run([str(exe),*map(str,args),'--output',str(p)],check=True,
                       env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),stdout=f,stderr=subprocess.STDOUT)
    return json.loads((p/'metadata.json').read_text())

@pytest.mark.parametrize('L,K',[(4,11),(10,101),(12,101)])
def test_momenta(L,K):
    grid,groups=momentum_plan(L,K);seen=[]
    for phase,points in groups:
        for j,q in points:
            np.testing.assert_allclose(-2*np.cos(grid[j]),-2*np.cos(2*np.pi*q/L+phase),atol=2e-12)
            seen.append(j)
    assert sorted(seen)==list(range(K))
    if L==10 and K==101:assert len(groups)==11

@pytest.mark.parametrize('key,value',[('L',4.5),('twists',3),('physical_densities',[1]),('eta',0),('plot_y_max',20)])
def test_invalid_config(tmp_path,key,value):
    d=json.loads((ROOT/'configs/smoke.json').read_text());d[key]=value
    p=tmp_path/'config.json';p.write_text(json.dumps(d))
    with pytest.raises(ValueError):load_config(p)

@pytest.mark.parametrize('nu,nd,U,phase,observable',[
 (2,2,8.,0.,'density'),(2,1,8.,0.,'density'),(1,1,0.,0.,'density'),
 (2,1,4.,.137,'density'),(2,2,8.,0.,'orbital-pairs'),(2,1,4.,.137,'orbital-pairs')])
def test_independent_neutral(exe,tmp_path,nu,nd,U,phase,observable):
    p=tmp_path/'native';omega=np.linspace(-1,16,1701);eta=.2
    run_native(exe['response'],p,['--L',4,'--nup',nu,'--ndn',nd,'--U',U,'--phase',phase,'--depth',120,'--observable',observable])
    actual=reconstruct(p,omega,eta);exact=independent.independent_ed(4,nu,nd,U,phase,omega,eta,observable)
    np.testing.assert_allclose(actual['metadata']['E0'],exact['E0'],atol=2e-9,rtol=0)
    for ch in ['charge','spin']:
        np.testing.assert_allclose(actual[ch],exact['spectra'][ch],atol=2e-8,rtol=2e-8)
        np.testing.assert_allclose(actual[f'weight_{ch}'],exact['weights'][ch],atol=1e-8,rtol=1e-8)
        np.testing.assert_allclose(actual[f'elastic_{ch}'],exact['elastic'][ch],atol=1e-8,rtol=1e-8)


def dense_hamiltonian(L,nu,nd,U,phase):
    states=[u|(d<<L) for u in independent.mask_states(L,nu) for d in independent.mask_states(L,nd)]
    where={s:i for i,s in enumerate(states)};H=np.zeros((len(states),len(states)),complex)
    for col,s in enumerate(states):
        H[col,col]=U*sum(((s>>j)&1)*((s>>(j+L))&1) for j in range(L))
        for spin in [0,1]:
            for j in range(L):
                a,b=j+spin*L,(j+1)%L+spin*L
                for dest,src,coef in [(a,b,-np.exp(1j*phase)),(b,a,-np.exp(-1j*phase))]:
                    if not(s>>src)&1 or (s>>dest)&1:continue
                    sign=(-1)**((s&((1<<src)-1)).bit_count());v=s^(1<<src)
                    sign*=(-1)**((v&((1<<dest)-1)).bit_count());v|=1<<dest
                    H[where[v],col]+=coef*sign
    return states,H

@pytest.mark.parametrize('U,phase',[(0.,.137),(4.,.137),(8.,0.)])
def test_independent_charged(exe,tmp_path,U,phase):
    L=4;nu=nd=1
    states,H=dense_hamiltonian(L,nu,nd,U,phase);E,V=np.linalg.eigh(H)
    assert E[1]-E[0]>1e-8
    ground=run_native(exe['response'],tmp_path/'ground',['--L',L,'--U',U,'--nup',nu,'--ndn',nd,'--phase',phase,'--ground-only'])
    assert len(ground['ground_states'])==1
    Q=ground['ground_states'][0]['Q']
    p=tmp_path/'charged'
    m=run_native(exe['spectral'],p,['--L',L,'--U',U,'--nup',nu,'--ndn',nd,'--phase',phase,'--ground-q',Q,'--depth',160])
    np.testing.assert_allclose(m['E0'],E[0],atol=1e-9,rtol=0)
    z=np.linspace(-6,16,301)+.29+1j*.1
    for spin in [0,1]:
      for particle in [False,True]:
        nup=nu+(1 if particle else -1)*(spin==0);ndn=nd+(1 if particle else -1)*(spin==1)
        targets,Hp=dense_hamiltonian(L,nup,ndn,U,phase);where={s:i for i,s in enumerate(targets)}
        Ep,Vp=np.linalg.eigh(Hp)
        for q in range(L):
            source=np.zeros(len(targets),complex)
            for col,s in enumerate(states):
                for j in range(L):
                    site=j+spin*L
                    if bool((s>>site)&1)==particle:continue
                    sign=(-1)**((s&((1<<site)-1)).bit_count());v=s^(1<<site)
                    source[where[v]]+=V[col,0]*sign*np.exp(((-1 if particle else 1)*2j*np.pi*q*j/L))/np.sqrt(L)
            weights=abs(Vp.conj().T@source)**2
            exact=np.sum(weights[:,None]/(z[None,:]+(-1 if particle else 1)*(Ep[:,None]-E[0])),axis=0)
            branch='particle' if particle else 'hole'
            actual=fraction_G(read_fraction(p/f'{branch}_spin0_g0_q{q}.dat'),z,particle)
            np.testing.assert_allclose(actual,exact,atol=2e-8,rtol=2e-8)
    import csv
    rows=list(csv.DictReader((p/'sum_rules.csv').open()))
    assert max(float(r['weight_error']) for r in rows)<1e-8
    assert max(float(r['moment_error']) for r in rows)<1e-8
