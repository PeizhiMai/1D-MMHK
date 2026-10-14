#!/usr/bin/env python3
"""Direct uniform 101-point A(k,w), preserving the 50-twist campaign's mu.

The required folded grid has 20 phases. Shared phases reuse certified ground
ensembles/sources. At new phases, a rigorous hopping norm bound prunes number
sectors, followed by all-Q energy solves and full spin/ground projector traces.
No interpolation of A and no chemical-potential retuning.
"""
import csv,json,os,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
import run_single_particle_scan as base
from postprocess_response import read_fraction
from twist_averaged_velocity_scan import Runs
from replot_twist50_response import digest,tag,HERE

ROOT=base.ROOT;OLD=base.OUT;OUT=ROOT/'single_particle_eta0p1_k101_depth600'
base.OUT=OUT
KGRID=.2*np.pi*(np.arange(20)-10)/20
REPS=[0]+list(range(10,20))
OMEGA=base.OMEGA
EXE=base.EXE
ENERGY_EXE=HERE/'bin/chain_Hubbard_tbc_response.x'
ENERGY_SOURCE=HERE.parent/'create_code/chain_Hubbard_tbc_response.cc'

def compatible_source(U,old_index,nu,nd,Q,phase,E0):
    p=OLD/'canonical'/f'U{U}'/f'twist{old_index:03d}'/f'Nup{nu}_Ndn{nd}_Q{Q}'
    if not (p/'COMPLETE').exists():return None
    line=json.loads((p/'lineage.json').read_text());m=json.loads((p/'metadata.json').read_text())
    assert line['source_sha256']==digest(HERE.parent/'create_code/chain_Hubbard_tbc_spectral.cc')
    assert line['helper_sha256']==digest(ENERGY_SOURCE) and line['executable_sha256']==digest(EXE)
    assert (m['L'],m['U'],m['Nup'],m['Ndn'],m['ground_Q'],m['depth'])==(10,U,nu,nd,Q,600)
    assert abs(m['phase_per_bond']-phase)<1e-12 and abs(m['E0']-E0)<1e-8
    return p

def phase_task(U,ik,reference,table_k,table_E):
    k=float(KGRID[ik]);nearest=int(np.argmin(abs(table_k-k)));delta=float(abs(table_k[nearest]-k))
    shared=delta<1e-13;phase=-k
    runs=Runs(OUT/f'U{U}',U,10,300,1,ENERGY_EXE,digest(ENERGY_SOURCE),digest(ENERGY_EXE))
    paths={};times={};details={};memo={}
    def energy(nu,nd):
        key=nu,nd
        if key not in memo:
            _,memo[key]=runs.run(ik,k,nu,nd,True)
        return memo[key]
    def charged(nu,nd,Q,E0):
        key=U,ik,nu,nd,Q
        if key not in paths:
            p=compatible_source(U,nearest,nu,nd,Q,phase,E0) if shared else None
            if p is None:_,p=base.job_run(key,phase,E0,600)
            paths[key]=p
        return paths[key]
    if shared:
        for n in [.7,.8,.9]:
            t=reference[n]['twists'][nearest]
            components=[]
            for c in t['components']:
                z={key:c[key] for key in ['N','Nup','Ndn','spin_flip_multiplicity','ground_dimension','probability','ground_states']}
                for g in c['ground_states']:charged(c['Nup'],c['Ndn'],g['Q'],g['E0'])
                components.append(z)
            details[n]=dict(twist_index=ik,k=k,phase_per_bond=phase,components=components,N=t['N'])
    else:
        # Two spin species, L bonds: ||H(phi+d)-H(phi)|| <=4Lt|sin(d/2)|.
        # Every number sector has this same bound. Include a numerical margin.
        bound=40*np.sin(delta/2)+1e-7
        numbers=np.arange(21);candidate_by_n={};all_candidates=set()
        for n in [.7,.8,.9]:
            mu=reference[n]['mu_selection']['mu'];old_grand=table_E[nearest]-mu*numbers
            candidate_by_n[n]=np.flatnonzero(old_grand-old_grand.min()<=2*bound+1e-8).tolist()
            all_candidates.update(candidate_by_n[n])
        energies={N:energy((N+1)//2,N//2)['E0'] for N in sorted(all_candidates)}
        for n in [.7,.8,.9]:
            mu=reference[n]['mu_selection']['mu'];candidate=candidate_by_n[n]
            minimum=min(energies[N]-mu*N for N in candidate)
            excluded=[int(N) for N in numbers if N not in candidate]
            assert all(table_E[nearest,N]-bound-mu*N>minimum-1e-8 for N in excluded)
            active=[N for N in candidate if abs(energies[N]-mu*N-minimum)<1e-8]
            components=[]
            for N in active:
                for nu in range((N+1)//2,min(N,10)+1):
                    nd=N-nu
                    if nd<0 or nd>10:continue
                    m=energy(nu,nd)
                    if m['E0']>energies[N]+1e-8:break
                    assert abs(m['E0']-energies[N])<1e-8
                    gs=[]
                    for g in m['ground_states']:
                        p=charged(nu,nd,g['Q'],g['E0']);cm=json.loads((p/'metadata.json').read_text())
                        for item in cm['grounds']:gs.append(dict(Q=g['Q'],id=item['id'],E0=cm['E0'],residual=item['residual']))
                    components.append(dict(N=N,Nup=nu,Ndn=nd,spin_flip_multiplicity=1 if nu==nd else 2,
                                           ground_dimension=len(gs),ground_states=gs))
            dim=sum(c['spin_flip_multiplicity']*c['ground_dimension'] for c in components)
            for c in components:c['probability']=c['spin_flip_multiplicity']*c['ground_dimension']/dim
            mean_N=sum(c['probability']*c['N'] for c in components)
            details[n]=dict(twist_index=ik,k=k,phase_per_bond=phase,components=components,N=mean_N,
                nearest_old_phase_index=nearest,hopping_norm_bound=bound,candidate_number_sectors=candidate,
                actual_grand_ground_numbers=active,number_sectors_excluded_with_certified_bound=excluded)
    (OUT/f'U{U}'/f'phase{ik:02d}').mkdir(parents=True,exist_ok=True)
    (OUT/f'U{U}'/f'phase{ik:02d}'/'projectors.json').write_text(json.dumps(details,indent=2)+'\n')
    return ik,details,paths

def render_case(U,n,reference,phases,paths):
    dest=OUT/f'U{U}'/f'n{tag(n)}';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE').exists():return
    k=np.linspace(0,np.pi,101);mu=reference['mu_selection']['mu'];source_hashes={};occupation=0.
    norm_error=moment_error=duplicate_error=0.
    if U==0:
        eps=-2*np.cos(k)-mu;A=.1/np.pi/((OMEGA[None,:]-eps[:,None])**2+.1**2)
        sample_density=None
    else:
        values={};expected_occupation=0.
        for ik in REPS:
            t=phases[ik][n];bands=np.zeros((10,len(OMEGA)));occ=0.
            for c in t['components']:
                counts={Q:sum(g['Q']==Q for g in c['ground_states']) for Q in {g['Q'] for g in c['ground_states']}}
                for Q,count in counts.items():
                    p=paths[(U,ik,c['Nup'],c['Ndn'],Q)];cm=json.loads((p/'metadata.json').read_text())
                    assert cm['ground_dimension']==count
                    factor=c['probability']/c['ground_dimension'];spins=[0] if c['Nup']==c['Ndn'] else [0,1]
                    for gi in range(count):
                      for spin in spins:
                       for q in range(10):
                        for branch in ['particle','hole']:
                            f=p/f'{branch}_spin{spin}_g{gi}_q{q}.dat';source_hashes[str(f)]=digest(f);d=read_fraction(f)
                            G=base.fraction_G(d,OMEGA+mu+.1j,branch=='particle');bands[q]-=factor/len(spins)*G.imag/np.pi
                            if branch=='hole':occ+=factor/len(spins)*d['norm']**2
                            if d['steps']:
                                ev=eigh_tridiagonal(d['a'],d['b'][1:],eigvals_only=True,lapack_driver='sterf')
                                edge=ev.min()-d['E0']+(-mu if branch=='particle' else mu)
                                if edge<-1e-7:
                                    ev,v=eigh_tridiagonal(d['a'],d['b'][1:],lapack_driver='stemr');active=d['norm']**2*v[0]**2>1e-8
                                    if active.any():assert (ev[active]-d['E0']+(-mu if branch=='particle' else mu)).min()>-1e-7
                    with (p/'sum_rules.csv').open() as f:
                        for r in csv.DictReader(f):norm_error=max(norm_error,float(r['weight_error']));moment_error=max(moment_error,float(r['moment_error']))
            wk=.05 if ik in [0,10] else .1
            occupation+=wk*occ/10;expected_occupation+=wk*t['N']/20
            for q in range(10):
                physical=2*np.pi*q/10+t['k'];wrapped=(physical+np.pi)%(2*np.pi)-np.pi
                j=int(round(abs(wrapped)*100/np.pi));assert abs(abs(wrapped)-k[j])<1e-12
                values.setdefault(j,[]).append(bands[q])
        assert set(values)==set(range(101))
        A=np.array([np.mean(values[j],axis=0) for j in range(101)])
        duplicate_error=max([np.max(abs(v[0]-v[1])) for v in values.values() if len(v)>1] or [0.])
        # Selected depth tests give <0.1% of peak. Report numerical duplicate
        # discrepancies explicitly rather than claiming exact floating equality.
        assert duplicate_error<1e-3*max(A.max(),1.)
        assert abs(occupation-expected_occupation)<1e-8
        sample_density=2*occupation
    assert A.min()>=0 and np.isfinite(A).all()
    np.savez_compressed(dest/'Akw.npz',k=k,omega=OMEGA,A_per_spin=A,mu=mu,eta=.1,
        target_density=n,actual_density=reference['actual_density_Ne_over_L'],reference_density_twists=50,
        folded_sampling_twists=20,sample_density=np.nan if sample_density is None else sample_density)
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(4.8,3.6),layout='constrained')
    im=ax.imshow(A.T,origin='lower',extent=[0,np.pi,-6,16],aspect='auto',cmap='rainbow',vmin=0,vmax=1.5,interpolation='none')
    ax.axhline(0,color='white',ls='--',lw=.6);ax.set(xlabel=r'$k$',ylabel=r'$\omega/t$',xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'])
    ax.text(.03,.96,rf'$U/t={U},\ n_{{\rm ref}}={reference["actual_density_Ne_over_L"]:.3f}$',transform=ax.transAxes,va='top',color='white',fontsize=9)
    fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.07,pad=.15,extend='max',label=r'$A(k,\omega)$ per spin')
    fig.savefig(dest/'spectral_Akw.png',dpi=180);fig.savefig(dest/'spectral_Akw.svg');plt.close(fig)
    m=dict(U=U,target_n=n,reference_n_50_twists=reference['actual_density_Ne_over_L'],sample_n_20_folded_phases=sample_density,
        mu=mu,eta=.1,k_points=101,k_grid='j*pi/100, j=0..100',interpolated=False,
        definition='A=-Im(Gup+Gdn)/(2*pi), unit total weight per spin',source_norm_error=norm_error,
        source_first_moment_error=moment_error,duplicate_inversion_max_abs_error=float(duplicate_error),source_hashes=source_hashes,
        mu_policy='unchanged chemical potential from original50-phase density constraint; new grid is spectral sampling only')
    (dest/'metadata.json').write_text(json.dumps(m,indent=2)+'\n');(dest/'COMPLETE').write_text('success\n')
    print(f'SPECTRAL_101_COMPLETE U={U} n={n:g}',flush=True)

def main():
    OUT.mkdir(exist_ok=True);start=time.monotonic()
    hashes={str(f):digest(f) for f in [Path(__file__),HERE/'run_single_particle_scan.py',ENERGY_SOURCE,HERE.parent/'create_code/chain_Hubbard_tbc_spectral.cc',EXE,ENERGY_EXE]}
    config=dict(eta=.1,depth=600,k_count=101,k_grid='linspace(0,pi,101)',interpolate=False,
       mu_reference_twists=50,spectral_folded_twists=20,representatives=REPS,source_sha256=hashes)
    cfg=OUT/'config.json'
    if cfg.exists():assert json.loads(cfg.read_text())==config
    else:cfg.write_text(json.dumps(config,indent=2)+'\n')
    for U in range(13):
        reference={n:json.loads((ROOT/f'U{U}'/f'n{tag(n)}'/'metadata.json').read_text()) for n in [.7,.8,.9]}
        phases={};paths={}
        if U:
            with np.load(ROOT/f'U{U}'/'ground_energy_table.npz') as table:
                oldk=table['k'].copy();E=table['E0_by_N'].copy()
            with ThreadPoolExecutor(max_workers=8) as pool:
                futures=[pool.submit(phase_task,U,ik,reference,oldk,E) for ik in REPS]
                for f in as_completed(futures):
                    ik,d,p=f.result();phases[ik]=d;paths.update(p)
                    print(f'U={U}: folded phase {len(phases)}/11, elapsed={time.monotonic()-start:.1f}s',flush=True)
        for n in [.7,.8,.9]:render_case(U,n,reference[n],phases,paths)
    assert all(digest(f)==sha for f,sha in hashes.items())
    (OUT/'COMPLETE').write_text('39 direct101-point spectra completed\n');print('ALL_101_SPECTRA_COMPLETE',flush=True)
if __name__=='__main__':main()
