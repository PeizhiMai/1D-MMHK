#!/usr/bin/env python3
"""Single-particle A(k,w), with the completed common-mu response ensembles.

Use inversion-related TBC blocks once. All 251 distinct momenta on 0..pi
of the L*50 grid are covered, with no change of mu or ground-state weights.
"""
import argparse,csv,json,os,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from replot_twist50_response import digest,tag,HERE
from postprocess_response import read_fraction

ROOT=HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
OUT=ROOT/'single_particle_eta0p1_depth600'
EXE=HERE/'bin/chain_Hubbard_tbc_spectral.x'
REPS=[0]+list(range(25,50))
OMEGA=np.linspace(-6,16,2201)

def fraction_G(d,z,particle):
    if d['steps']==0:return np.zeros(z.shape,dtype=complex)
    a,b=d['a'],d['b'];shift=z+d['E0'] if particle else z-d['E0']
    sign=-1 if particle else 1
    tail=shift+sign*a[-1]
    for i in range(len(a)-2,-1,-1):tail=shift+sign*a[i]-b[i+1]**2/tail
    return d['norm']**2/tail

def job_run(key,phase,expected,depth):
    U,ik,nu,nd,Q=key
    output=OUT/'canonical'/f'U{U}'/f'twist{ik:03d}'/f'Nup{nu}_Ndn{nd}_Q{Q}'
    command=[str(EXE),'--L','10','--U',str(U),'--nup',str(nu),'--ndn',str(nd),
             '--phase',repr(phase),'--ground-q',str(Q),'--depth',str(depth),'--output',str(output)]
    lineage=dict(command=command,source_sha256=digest(HERE.parent/'create_code/chain_Hubbard_tbc_spectral.cc'),
                 helper_sha256=digest(HERE.parent/'create_code/chain_Hubbard_tbc_response.cc'),executable_sha256=digest(EXE))
    if (output/'COMPLETE').exists():
        assert json.loads((output/'lineage.json').read_text())==lineage
    else:
        output.mkdir(parents=True,exist_ok=True)
        if (output/'metadata.json').exists():raise RuntimeError(f'Interrupted source run: {output}')
        (output/'lineage.json').write_text(json.dumps(lineage,indent=2)+'\n')
        start=time.monotonic()
        with (output/'run.log').open('w') as f:
            subprocess.run(command,env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'),stdout=f,stderr=subprocess.STDOUT,check=True)
        (output/'timing.json').write_text(json.dumps(dict(seconds=time.monotonic()-start)))
    m=json.loads((output/'metadata.json').read_text())
    assert abs(m['E0']-expected)<1e-8
    assert max(g['residual'] for g in m['grounds'])<1e-8
    return key,output

def load_plan(U_list):
    cases={};jobs={};hashes={}
    for U in U_list:
      for n in [.7,.8,.9]:
        file=ROOT/f'U{U}'/f'n{tag(n)}'/'metadata.json';m=json.loads(file.read_text());hashes[str(file)]=digest(file);cases[(U,n)]=m
        if U==0:continue
        for ik in range(1,25):
            a,b=m['twists'][ik],m['twists'][50-ik]
            assert a['N']==b['N'] and a['ground_dimension']==b['ground_dimension']
            signature=lambda t:sorted((c['Nup'],c['Ndn'],c['ground_dimension'],round(c['probability'],12)) for c in t['components'])
            assert signature(a)==signature(b)
            for c,cc in zip(a['components'],b['components']):
                assert abs(min(g['E0'] for g in c['ground_states'])-min(g['E0'] for g in cc['ground_states']))<1e-8
        for ik in REPS:
          t=m['twists'][ik]
          for c in t['components']:
            for g in c['ground_states']:
              key=U,ik,c['Nup'],c['Ndn'],g['Q']
              jobs[key]=(t['phase_per_bond'],g['E0'])
    return cases,jobs,hashes

def plot_case(U,n,m,paths):
    dest=OUT/f'U{U}'/f'n{tag(n)}';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE').exists():return
    mu=m['mu_selection']['mu'];eta=.1;k=np.linspace(0,np.pi,251)
    source_hashes={};norm_error=moment_error=negative_pole=0.;occupation=0.
    if U==0:
        energies=-2*np.cos(k)-mu
        A=eta/np.pi/((OMEGA[None,:]-energies[:,None])**2+eta**2)
        occupation=m['actual_density_Ne_over_L']/2
    else:
        values={};weights={}
        for ik in REPS:
            t=m['twists'][ik];bands=np.zeros((10,len(OMEGA)));occ=0.
            for c in t['components']:
                counts={Q:sum(g['Q']==Q for g in c['ground_states']) for Q in {g['Q'] for g in c['ground_states']}}
                for Q,count in counts.items():
                    source=paths[(U,ik,c['Nup'],c['Ndn'],Q)]
                    cm=json.loads((source/'metadata.json').read_text());assert cm['ground_dimension']==count
                    factor=c['probability']/c['ground_dimension']
                    spins=[0] if c['Nup']==c['Ndn'] else [0,1]
                    for gi in range(count):
                      for spin in spins:
                        for q in range(10):
                          for branch in ['particle','hole']:
                            f=source/f'{branch}_spin{spin}_g{gi}_q{q}.dat';source_hashes[str(f)]=digest(f);d=read_fraction(f)
                            G=fraction_G(d,OMEGA+mu+1j*eta,branch=='particle')
                            bands[q]+=factor/len(spins)*(-G.imag/np.pi)
                            if branch=='hole':occ+=factor/len(spins)*d['norm']**2
                            if d['steps']:
                                ev=eigh_tridiagonal(d['a'],d['b'][1:],eigvals_only=True,lapack_driver='sterf')
                                edge=ev.min()-d['E0']-mu if branch=='particle' else ev.min()-d['E0']+mu
                                # Spurious dark Ritz eigenvalues need a residue check before rejecting.
                                if edge<-1e-7:
                                    ev,vec=eigh_tridiagonal(d['a'],d['b'][1:],lapack_driver='stemr')
                                    active=d['norm']**2*vec[0]**2>1e-8
                                    if active.any():
                                        e=(ev[active]-d['E0']-mu if branch=='particle' else ev[active]-d['E0']+mu).min()
                                        negative_pole=min(negative_pole,float(e));assert e>-1e-7,(U,n,ik,q,branch,e)
                    with (source/'sum_rules.csv').open() as f:
                        for r in csv.DictReader(f):
                            norm_error=max(norm_error,float(r['weight_error']));moment_error=max(moment_error,float(r['moment_error']))
            wk=.02 if ik in [0,25] else .04
            occupation+=wk*occ/10
            for q in range(10):
                physical=2*np.pi*q/10+t['k'];wrapped=(physical+np.pi)%(2*np.pi)-np.pi
                j=int(round(abs(wrapped)*250/np.pi));assert abs(abs(wrapped)-k[j])<1e-12
                values.setdefault(j,[]).append(bands[q]);weights.setdefault(j,[]).append(1.)
        assert set(values)==set(range(251))
        A=np.array([np.mean(values[j],axis=0) for j in range(251)])
        # Self-inverse boundary/zero-twist samples can contribute twice to one |k|.
        duplicates=max([np.max(np.abs(v[0]-v[1])) for v in values.values() if len(v)>1] or [0.])
        assert duplicates<1e-5,duplicates
        assert abs(occupation-m['actual_density_Ne_over_L']/2)<1e-8
    assert A.min()>-1e-10 and np.isfinite(A).all()
    np.savez_compressed(dest/'Akw.npz',k=k,omega=OMEGA,A_per_spin=A,mu=mu,eta=eta,
                        target_density=n,actual_density=m['actual_density_Ne_over_L'])
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(4.8,3.6),layout='constrained')
    im=ax.imshow(A.T,origin='lower',extent=[0,np.pi,OMEGA[0],OMEGA[-1]],aspect='auto',cmap='rainbow',vmin=0,vmax=1.5,interpolation='none')
    ax.axhline(0,color='white',lw=.6,ls='--');ax.set(xlabel=r'$k$',ylabel=r'$\omega/t$',xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'])
    ax.text(.03,.96,rf'$U/t={U},\ n={m["actual_density_Ne_over_L"]:.3f}$',transform=ax.transAxes,va='top',color='white',fontsize=9)
    fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.07,pad=.15,extend='max',label=r'$A(k,\omega)$ per spin')
    fig.savefig(dest/'spectral_Akw.png',dpi=180);fig.savefig(dest/'spectral_Akw.svg');plt.close(fig)
    (dest/'metadata.json').write_text(json.dumps(dict(U=U,target_n=n,actual_n=m['actual_density_Ne_over_L'],mu=mu,eta=.1,
        definition='A=-Im(Gup+Gdn)/(2*pi); integral over all frequency=1',
        physical_k='K-phi=2*pi*q/L+k_RBZ',particle_poles='Em(N+1)-E0-mu',hole_poles='E0-Em(N-1)-mu',
        ground_source=str(ROOT/f'U{U}'/f'n{tag(n)}'/'metadata.json'),source_hashes=source_hashes,
        distinct_positive_momenta=251,original_twists=50,representatives_used=26,inversion_reduction=True,
        pole_occupation_per_spin=occupation,source_norm_error=norm_error,source_first_moment_error=moment_error,
        minimum_significant_grand_excitation=negative_pole,colour_scale=[0,1.5],omega_window=[-6,16]),indent=2)+'\n')
    (dest/'COMPLETE').write_text('success\n');print(f'SPECTRAL_COMPLETE U={U} n={n:g}',flush=True)

def panel_plots():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for n in [.7,.8,.9]:
      if not all((OUT/f'U{U}'/f'n{tag(n)}'/'COMPLETE').exists() for U in range(13)):continue
      fig,axes=plt.subplots(4,4,figsize=(10.8,8),layout='constrained')
      for U,ax in zip(range(13),axes.flat):
        with np.load(OUT/f'U{U}'/f'n{tag(n)}'/'Akw.npz') as d:
          im=ax.imshow(d['A_per_spin'].T,origin='lower',extent=[0,np.pi,-6,16],aspect='auto',cmap='rainbow',vmin=0,vmax=1.5,interpolation='none')
          ax.axhline(0,color='white',lw=.5,ls='--');ax.set(xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'],yticks=[-5,0,5,10,15])
          ax.text(.04,.96,rf'$U/t={U}$',transform=ax.transAxes,va='top',fontsize=9,color='white')
      for ax in list(axes.flat)[13:]:ax.set_visible(False)
      fig.suptitle(rf'Single-particle $A(k,\omega)$: target $n={n:g}$, 10-MMHK, $\eta/t=0.1$',fontsize=13)
      fig.supxlabel(r'$k$');fig.supylabel(r'$\omega/t$')
      fig.colorbar(im,ax=list(axes.flat)[:13],orientation='horizontal',fraction=.035,pad=.03,extend='max',label=r'$A(k,\omega)$ per spin')
      fig.savefig(OUT/f'spectral_panel_n{tag(n)}.png',dpi=180);fig.savefig(OUT/f'spectral_panel_n{tag(n)}.svg');plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workers',type=int,default=8);p.add_argument('--U',type=int,nargs='+',default=list(range(13)));p.add_argument('--depth',type=int,default=600);a=p.parse_args()
    OUT.mkdir(exist_ok=True);cases,jobs,hashes=load_plan(a.U)
    config=dict(U_list=a.U,depth=a.depth,eta=.1,physical_n=[.7,.8,.9],original_twists=50,representatives=REPS,
      source_sha256=digest(HERE.parent/'create_code/chain_Hubbard_tbc_spectral.cc'),helper_sha256=digest(HERE.parent/'create_code/chain_Hubbard_tbc_response.cc'),exe_sha256=digest(EXE),runner_sha256=digest(__file__),ground_input_sha256=hashes)
    cfg=OUT/'config.json'
    if cfg.exists():assert json.loads(cfg.read_text())==config
    else:cfg.write_text(json.dumps(config,indent=2)+'\n')
    paths={};start=time.monotonic();print('Canonical jobs:',len(jobs),flush=True)
    for U in a.U:
      todo={key:value for key,value in jobs.items() if key[0]==U}
      with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(job_run,key,*value,a.depth) for key,value in todo.items()]
        for i,f in enumerate(as_completed(futures),1):
          key,where=f.result();paths[key]=where
          if i%10==0 or i==len(todo):print(f'U={U}: canonical {i}/{len(todo)}, elapsed={time.monotonic()-start:.1f}s',flush=True)
      for n in [.7,.8,.9]:plot_case(U,n,cases[(U,n)],paths)
    panel_plots()
    assert all(digest(f)==sha for f,sha in hashes.items())
    (OUT/'COMPLETE').write_text('All requested spectral cases completed\n')
    print('ALL_SPECTRAL_COMPLETE',flush=True)
if __name__=='__main__':main()
