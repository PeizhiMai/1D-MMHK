"""Fresh common-mu ground scans, physical neutral responses and charged G.

No historical runs, original workstation paths, book PDFs, or student code
are needed. The original audited helpers/operators are imported unchanged.
Only newly written provenance is normalized to paths relative to this run.
"""
import csv
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'response'))
sys.path.insert(0, str(ROOT/'skills/mmhk-1d-spin-charge/scripts'))
import numpy as np
from scipy.linalg import eigh_tridiagonal
from twist_averaged_velocity_scan import Runs, energy_table, ground_components, save_density
from mmhk_density_workflow import k_grid, choose_mu
from extract_zero_twist_velocities import merge_poles, write_csv
from run_response import free_case
from postprocess_response import read_fraction
from run_single_particle_scan import fraction_G


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tag(value):
    return f'{float(value):g}'.replace('.', 'p')


def dump(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def load_config(path):
    c=json.loads(Path(path).read_text())
    required=['L','physical_densities','U_list','twists','density_tolerance','response_depth',
              'spectral_depth','k_points','eta','omega_min','omega_max','omega_points',
              'plot_y_min','plot_y_max','ba_step','ba_nodes','workers','threads']
    if any(k not in c for k in required):raise ValueError('Incomplete configuration')
    for key in ['L','twists','response_depth','spectral_depth','k_points','omega_points','ba_nodes','workers','threads']:
        if isinstance(c[key],bool) or not isinstance(c[key],int):raise ValueError(f'{key} must be an integer')
    L=c['L'];nk=c['twists']
    if L<4 or L%2 or (nk!=1 and (nk<4 or nk%2)):raise ValueError('Require even L>=4 and twists=1 or even>=4')
    if any(not 0<n<1 for n in c['physical_densities']) or not c['physical_densities']:raise ValueError('Physical density must be 0<n<1')
    if len(set(c['physical_densities']))!=len(c['physical_densities']):raise ValueError('Duplicate densities')
    if not c['U_list'] or any(U<0 for U in c['U_list']) or len(set(c['U_list']))!=len(c['U_list']):raise ValueError('Invalid coupling list')
    if not all(np.isfinite(v) for v in [*c['U_list'],*c['physical_densities'],*map(c.get,['eta','density_tolerance','omega_min','omega_max','ba_step'])]):raise ValueError('Nonfinite configuration')
    if min(c['workers'],c['threads'])<1 or min(c['response_depth'],c['spectral_depth'])<3 or c['k_points']<3 or c['omega_points']<3:raise ValueError('Invalid resources/grid')
    if c['eta']<=0 or c['density_tolerance']<0 or c['omega_min']>=c['omega_max'] or c['ba_nodes']<50 or c['ba_step']<=0:raise ValueError('Invalid analysis settings')
    if not c['omega_min']<=c['plot_y_min']<c['plot_y_max']<=c['omega_max']:raise ValueError('Plot range must lie inside stored frequency range')
    dimension=math.comb(L,L//2)**2
    if dimension*(1+4*L)>=2**31:raise ValueError('Conservative CSR size bound exceeds the portable LP64 index range')
    if dimension>c.get('max_canonical_dimension',1_000_000):raise ValueError(f'Canonical dimension {dimension} exceeds laptop resource gate; explicitly increase max_canonical_dimension to opt in')
    return c


def solvers(directory=None):
    d=Path(directory or os.environ.get('MMHK_SOLVER_DIR',ROOT/'.build/solver/bin')).resolve()
    result={kind:d/f'chain_Hubbard_tbc_{kind}' for kind in ['response','spectral']}
    for p in result.values():
        if not p.is_file():raise FileNotFoundError(f'Build portable solvers first; missing {p}')
    return result


class PortableRuns(Runs):
    def run(self, ik, k, Nu, Nd, ground_only=True):
        kind='energy_scan' if ground_only else 'qmin_responses'
        p=self.root/kind/f'twist{ik:03d}'/f'Nup{Nu}_Ndn{Nd}'
        identity=dict(L=self.L,U=self.U,Nup=Nu,Ndn=Nd,phase_per_bond=float(-k),
                      depth=self.depth,ground_only=ground_only,q_indices=[] if ground_only else [1],
                      source_sha256=self.source_hash,executable_sha256=self.exe_hash)
        if (p/'COMPLETE').exists():
            if json.loads((p/'lineage.json').read_text())!=identity:raise RuntimeError(f'Incompatible cached run: {p}')
        else:
            if (p/'metadata.json').exists():raise RuntimeError(f'Interrupted native run; inspect before retrying: {p}')
            p.mkdir(parents=True,exist_ok=True);dump(p/'lineage.json',identity)
            command=[str(self.exe),'--L',str(self.L),'--U',str(self.U),'--nup',str(Nu),'--ndn',str(Nd),
                     '--phase',repr(float(-k)),'--depth',str(self.depth),'--observable','density','--output',str(p.resolve())]
            command+=['--ground-only'] if ground_only else ['--q','1']
            start=time.monotonic()
            with (p/'run.log').open('w') as f:subprocess.run(command,env=dict(self.env,OPENBLAS_NUM_THREADS='1'),stdout=f,stderr=subprocess.STDOUT,check=True)
            with self.lock:self.new_calls+=1;self.seconds+=time.monotonic()-start
        m=json.loads((p/'metadata.json').read_text())
        if (m['L'],m['U'],m['Nup'],m['Ndn'])!=(self.L,self.U,Nu,Nd) or abs(m['phase_per_bond']+k)>1e-12 or not m['ground_scan_all_Q']:raise RuntimeError('Wrong native response metadata')
        return p,m


def normalize_new_metadata(path, root):
    m=json.loads(path.read_text())
    for key in ['raw_fraction_sha256','source_hashes']:
        if key in m:m[key]={str(Path(p).relative_to(root)):s for p,s in m[key].items()}
    for check in m.get('source_moment_checks',[]):
        check['file']=str(Path(check['file']).relative_to(root))
    for block in m.get('twists',[]):
        for component in block['components']:
            if 'directory' in component:component['directory']=str(Path(component['directory']).relative_to(root))
    dump(path,m)


def initialize(c,out,exe):
    out=Path(out).resolve()
    archive=ROOT/'response/results'
    if out==archive or archive in out.parents:raise ValueError('Fresh calculations must not overwrite the historical archive; use runs/NAME')
    source_files=[ROOT/p for p in ['create_code/chain_Hubbard_tbc_response.cc','create_code/chain_Hubbard_tbc_spectral.cc','reproducibility/workflow.py','reproducibility/compat/mkl.h','cmake/patch_qbasis.py','response/twist_averaged_velocity_scan.py','response/mmhk_density_workflow.py','response/extract_zero_twist_velocities.py','response/postprocess_response.py','response/run_single_particle_scan.py','skills/mmhk-1d-spin-charge/scripts/run_response.py']]
    config=dict(c,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_files},
                executable_sha256={k:sha(v) for k,v in exe.items()},
                qbasis_commit='33324ce6081145ea9af981c64b2a2cc32c7faf0a',backend='OpenBLAS/LAPACKE/ARPACK LP64',
                ensemble='complete grand ground trace, common mu, no forced N or Sz',
                density_convention='Ne/L, half filling=1')
    p=out/'config.json'
    if p.exists() and json.loads(p.read_text())!=config:raise RuntimeError('Configuration/source/build differs; use a fresh output root')
    dump(p,config)
    return out


def responses(c,out,exe):
    k,w=k_grid(c['L'],c['twists'],'periodic');statistics=[];density_rows=[]
    with ThreadPoolExecutor(max_workers=c['workers']) as pool:
      for U in c['U_list']:
        root=out/f'U{U:g}';root.mkdir(exist_ok=True)
        config=json.loads((out/'config.json').read_text())
        config.update(solver_source_sha256=sha(ROOT/'create_code/chain_Hubbard_tbc_response.cc'),
                      solver_executable_sha256=sha(exe['response']))
        if U==0:
            for n in c['physical_densities']:
                dest=root/f'n{tag(n)}'
                if (dest/'COMPLETE').exists():continue
                stage=root/f'.postprocess_n{tag(n)}';stage.mkdir(exist_ok=True)
                # A true full free grand-ground projector is allowed. Arbitrary
                # shell weights are never selected to fabricate a target.
                try:
                    selection=free_case(stage,c['L'],k,n,config,merge_poles,write_csv,tag)
                    actual=selection['actual_density_Ne_over_L']
                except RuntimeError as error:
                    if 'target' not in str(error) and 'attainable' not in str(error):raise
                    eps=-2*np.cos(2*np.pi*np.arange(c['L'])[None,:]/c['L']+k[:,None])
                    levels=np.sort(np.repeat(eps,2,axis=1),axis=1)
                    E=np.column_stack([np.zeros(len(k)),np.cumsum(levels,axis=1)])
                    selection=choose_mu(E,w,n/2,c['L']);actual=2*selection['filling']
                    if abs(actual-n)>c['density_tolerance']+1e-12:raise RuntimeError(f'U=0 target {n} unattainable within explicit density tolerance') from error
                    free_case(stage,c['L'],k,actual,config,merge_poles,write_csv,tag)
                produced=stage/f'n{tag(actual)}'
                if dest.exists():shutil.rmtree(dest)
                shutil.copytree(produced,dest)
                shutil.copy2(stage/'ground_energy_table.npz',root/'ground_energy_table.npz')
                shutil.rmtree(stage)
                density_rows.append(dict(U=U,target_density=n,actual_density=actual,mu=selection['mu'],error=actual-n))
                m=json.loads((dest/'metadata.json').read_text());m['target_density_Ne_over_L']=n
                m['density_tolerance']=c['density_tolerance'];dump(dest/'metadata.json',m)
            continue
        runs=PortableRuns(root,U,c['L'],c['response_depth'],c['threads'],exe['response'],
                          sha(ROOT/'create_code/chain_Hubbard_tbc_response.cc'),sha(exe['response']))
        if (root/'ground_energy_table.npz').exists():
            with np.load(root/'ground_energy_table.npz') as d:E=d['E0_by_N'].copy();np.testing.assert_allclose(d['k'],k,atol=1e-14)
        elif len(k)==1:
            E=np.array([[runs.run(0,0.,(N+1)//2,N//2)[1]['E0'] for N in range(2*c['L']+1)]])
            np.savez_compressed(root/'ground_energy_table.npz',k=k,weights=w,E0_by_N=E)
        else:E=energy_table(runs,k,pool)
        selections={}
        for n in c['physical_densities']:
            s=choose_mu(E,w,n/2,c['L']);actual=2*s['filling']
            if abs(actual-n)>c['density_tolerance']+1e-12:raise RuntimeError(f'U={U}, target n={n}: nearest stable n={actual}, beyond density tolerance; increase twists or explicitly change tolerance')
            if not 0<actual<1:raise RuntimeError('Nearest stable density is outside the metallic BA comparison domain 0<n<1')
            selections[n]=s
            density_rows.append(dict(U=U,target_density=n,actual_density=actual,mu=s['mu'],error=actual-n))
        dump(root/'density_selections.json',selections)
        tasks=sorted({(i,N) for s in selections.values() for i,N in enumerate(s['numbers'])})
        futures={pool.submit(ground_components,runs,i,k[i],N,E[i,N]):(i,N) for i,N in tasks}
        components={futures[f]:f.result() for f in as_completed(futures)}
        for n,s in selections.items():
            dest=root/f'n{tag(n)}'
            if (dest/'COMPLETE').exists():continue
            actual=2*s['filling'];proxy=copy.copy(runs)
            proxy.root=root/f'.postprocess_n{tag(n)}';proxy.root.mkdir(exist_ok=True)
            save_density(proxy,actual,s,k,components)
            produced=proxy.root/f'n{tag(actual)}'
            if dest.exists():shutil.rmtree(dest)
            shutil.copytree(produced,dest);shutil.rmtree(proxy.root)
            m=json.loads((dest/'metadata.json').read_text());m.update(target_density_Ne_over_L=n,density_tolerance=c['density_tolerance'])
            dump(dest/'metadata.json',m);normalize_new_metadata(dest/'metadata.json',out)
        statistics.append(dict(U=U,native_calls=runs.new_calls,solver_seconds=runs.seconds))
        dump(out/'response_statistics.json',statistics)
    if density_rows:write_csv(out/'density_summary.csv',density_rows)
    (out/'RESPONSES_COMPLETE').write_text('all requested physical densities completed\n')


def momentum_plan(L,count):
    """Direct k=j*pi/(count-1), grouped into exact inversion-related blocks.

    Inversion is applied to the FULL unpolarized ground trace, not a selected
    member of a spin/momentum multiplet. No interpolation of A is performed.
    """
    groups={};grid=np.linspace(0,np.pi,count);spacing=2*np.pi/L
    for j,p in enumerate(grid):
        q=int(np.floor((p+spacing/2)/spacing+1e-12));k=p-q*spacing
        if abs(k)<1e-12:k=0.
        if k<0:q=-q;k=-k
        q%=L
        key=round(float(k),13)
        groups.setdefault(key,[]).append((j,q))
    return grid,sorted(groups.items())


def charged_job(c,out,exe,U,ik,k,Nu,Nd,Q,E0):
    p=out/'spectral/canonical'/f'U{U:g}'/f'phase{ik:03d}'/f'Nup{Nu}_Ndn{Nd}_Q{Q}'
    identity=dict(L=c['L'],U=U,phase=-k,Nup=Nu,Ndn=Nd,Q=Q,depth=c['spectral_depth'],
                  source_sha256=sha(ROOT/'create_code/chain_Hubbard_tbc_spectral.cc'),
                  helper_sha256=sha(ROOT/'create_code/chain_Hubbard_tbc_response.cc'),executable_sha256=sha(exe['spectral']))
    if (p/'COMPLETE').exists():
        if json.loads((p/'lineage.json').read_text())!=identity:raise RuntimeError('Incompatible charged cache')
    else:
        if (p/'metadata.json').exists():raise RuntimeError(f'Interrupted charged run needs inspection: {p}')
        p.mkdir(parents=True,exist_ok=True);dump(p/'lineage.json',identity)
        command=[str(exe['spectral']),'--L',str(c['L']),'--U',str(U),'--nup',str(Nu),'--ndn',str(Nd),
                 '--phase',repr(float(-k)),'--ground-q',str(Q),'--depth',str(c['spectral_depth']),'--output',str(p.resolve())]
        with (p/'run.log').open('w') as f:subprocess.run(command,env=dict(os.environ,OMP_NUM_THREADS=str(c['threads']),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'),check=True,stdout=f,stderr=subprocess.STDOUT)
    m=json.loads((p/'metadata.json').read_text())
    if abs(m['E0']-E0)>1e-8 or max(g['residual'] for g in m['grounds'])>1e-8:raise RuntimeError('Uncertified charged ground')
    return p,m


def spectral_phase(c,out,exe,U,ik,k,reference,oldk,E):
    L=c['L'];nearest=int(np.argmin(abs(oldk-k)));bound=4*L*abs(np.sin((oldk[nearest]-k)/2))+1e-7
    runs=PortableRuns(out/'spectral'/f'U{U:g}',U,L,c['response_depth'],c['threads'],exe['response'],
                      sha(ROOT/'create_code/chain_Hubbard_tbc_response.cc'),sha(exe['response']))
    memo={};sources={};states={};numbers=np.arange(2*L+1)
    def energy(Nu,Nd):
        if (Nu,Nd) not in memo:memo[Nu,Nd]=runs.run(ik,k,Nu,Nd)[1]
        return memo[Nu,Nd]
    candidates={}
    for n,m in reference.items():
        g=E[nearest]-m['mu_selection']['mu']*numbers
        candidates[n]=np.flatnonzero(g-g.min()<=2*bound+1e-8).tolist()
    ground_E={N:energy((N+1)//2,N//2)['E0'] for N in sorted(set(sum(candidates.values(),[])))}
    for n,m in reference.items():
        mu=m['mu_selection']['mu'];minimum=min(ground_E[N]-mu*N for N in candidates[n])
        excluded=[int(N) for N in numbers if N not in candidates[n]]
        if any(E[nearest,N]-bound-mu*N<=minimum-1e-8 for N in excluded):raise RuntimeError('Invalid number-sector exclusion')
        active=[N for N in candidates[n] if abs(ground_E[N]-mu*N-minimum)<1e-8];components=[]
        for N in active:
            for Nu in range((N+1)//2,min(N,L)+1):
                Nd=N-Nu
                if Nd<0 or Nd>L:continue
                e=energy(Nu,Nd)
                if e['E0']>ground_E[N]+1e-8:break
                if abs(e['E0']-ground_E[N])>1e-8:raise RuntimeError('Unexpected lower spin ground state')
                gs=[]
                for g in e['ground_states']:
                    key=(Nu,Nd,g['Q'])
                    if key not in sources:sources[key]=charged_job(c,out,exe,U,ik,k,*key,g['E0'])
                    p,cm=sources[key]
                    gs.extend(dict(Q=g['Q'],id=item['id'],directory=str(p.relative_to(out))) for item in cm['grounds'])
                components.append(dict(N=N,Nup=Nu,Ndn=Nd,spin_flip_multiplicity=1 if Nu==Nd else 2,
                                       ground_dimension=len(gs),ground_states=gs))
        dim=sum(x['spin_flip_multiplicity']*x['ground_dimension'] for x in components)
        for x in components:x['probability']=x['spin_flip_multiplicity']*x['ground_dimension']/dim
        states[str(n)]=dict(k=k,phase_per_bond=-k,components=components,
                           candidate_N=candidates[n],active_N=active,excluded_N=excluded,
                           hopping_norm_bound=bound,nearest_reference_twist=nearest)
    dump(out/'spectral'/f'U{U:g}'/f'phase{ik:03d}'/'projectors.json',states)
    return ik,k,states


def spectra(c,out,exe):
    grid,groups=momentum_plan(c['L'],c['k_points']);omega=np.linspace(c['omega_min'],c['omega_max'],c['omega_points'])
    for U in c['U_list']:
        reference={n:json.loads((out/f'U{U:g}'/f'n{tag(n)}'/'metadata.json').read_text()) for n in c['physical_densities']}
        phases={}
        if U:
            with np.load(out/f'U{U:g}'/'ground_energy_table.npz') as d:oldk=d['k'].copy();E=d['E0_by_N'].copy()
            with ThreadPoolExecutor(max_workers=c['workers']) as pool:
                tasks=[pool.submit(spectral_phase,c,out,exe,U,ik,k,reference,oldk,E) for ik,(k,_) in enumerate(groups)]
                for f in as_completed(tasks):
                    ik,k,s=f.result();phases[ik]=s;print(f'SPECTRAL GROUNDS U={U}: {len(phases)}/{len(groups)}',flush=True)
        for n,m in reference.items():
            p=out/'spectral'/f'U{U:g}'/f'n{tag(n)}';p.mkdir(parents=True,exist_ok=True)
            if (p/'COMPLETE').exists():continue
            mu=m['mu_selection']['mu'];A=np.zeros((len(grid),len(omega)));hashes={};plan=[];norm_error=moment_error=0.
            if U==0:A=c['eta']/np.pi/((omega[None,:]+mu+2*np.cos(grid[:,None]))**2+c['eta']**2)
            else:
              for ik,(_,points) in enumerate(groups):
                components=phases[ik][str(n)]['components']
                for component in components:
                    spins=[0] if component['Nup']==component['Ndn'] else [0,1]
                    factor=component['probability']/component['ground_dimension']/len(spins)
                    for g in component['ground_states']:
                        directory=out/g['directory']
                        with (directory/'sum_rules.csv').open() as f:
                            for r in csv.DictReader(f):
                                norm_error=max(norm_error,float(r['weight_error']));moment_error=max(moment_error,float(r['moment_error']))
                        for j,q in points:
                            for spin in spins:
                                for branch in ['particle','hole']:
                                    file=directory/f'{branch}_spin{spin}_g{g["id"]}_q{q}.dat'
                                    key=str(file.relative_to(out));hashes[key]=sha(file);fraction=read_fraction(file)
                                    A[j]-=factor*fraction_G(fraction,omega+mu+1j*c['eta'],branch=='particle').imag/np.pi
                                    plan.append(dict(j=j,file=key,factor=factor,particle=branch=='particle'))
                                    if fraction['steps']:
                                        ev=eigh_tridiagonal(fraction['a'],fraction['b'][1:],eigvals_only=True,lapack_driver='sterf')
                                        sign=-mu if branch=='particle' else mu
                                        if ev.min()-fraction['E0']+sign < -1e-7:
                                            ev,v=eigh_tridiagonal(fraction['a'],fraction['b'][1:],lapack_driver='stemr')
                                            visible=fraction['norm']**2*v[0]**2>1e-8
                                            if visible.any() and (ev[visible]-fraction['E0']+sign).min()<-1e-7:raise RuntimeError('Negative charged excitation with nonnegligible residue')
            if not np.isfinite(A).all() or (A<0).any() or norm_error>1e-8 or moment_error>1e-8:raise RuntimeError('Spectral positivity/sum-rule gate failed')
            np.savez_compressed(p/'Akw.npz',k=grid,omega=omega,A_per_spin=A,mu=mu,eta=c['eta'],
                                target_density=n,actual_density=m['actual_density_Ne_over_L'],reference_density_twists=c['twists'])
            dump(p/'metadata.json',dict(U=U,L=c['L'],target_n=n,reference_n=m['actual_density_Ne_over_L'],mu=mu,eta=c['eta'],
                  source_hashes=hashes,source_norm_error=norm_error,source_first_moment_error=moment_error,
                  k_points=len(grid),interpolated=False,mu_policy='reference density grid unchanged',
                  phase_count=len(groups),inversion_reuse='complete unpolarized ground trace'))
            dump(p/'spectral_plan.json',plan);(p/'COMPLETE').write_text('success\n')
            print(f'SPECTRAL COMPLETE U={U} n={n:g}',flush=True)
    (out/'SPECTRA_COMPLETE').write_text('all directly calculated momenta completed\n')
