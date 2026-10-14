"""Portable plot-only mode: numerical arrays/poles in, figures/tables out.

Never launches a many-body solver or inspects executables.
The BA reference uses integral equations, never a book PDF or a many-body run.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'response'))
sys.path.insert(0,str(ROOT/'response/essler_bethe_ansatz'))
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR',str(Path(os.environ.get('TMPDIR','/tmp'))/'mmhk-matplotlib'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fast_lorentzian_spectrum import spectrum
from postprocess_twist50_density_scan import peak_record
from hubbard_velocities import solve
from mmhk_density_workflow import choose_mu
from postprocess_response import read_fraction
from run_single_particle_scan import fraction_G

COLORS={'charge':'#c43c64','spin':'#234f91'}
STYLES={'charge':'-','spin':'--'}


def tag(n):return f'{float(n):g}'.replace('.','p')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_csv(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def save(fig,out,name):
    fig.savefig(out/f'{name}.png',dpi=180,facecolor='white')
    fig.savefig(out/f'{name}.svg',facecolor='white');plt.close(fig)


def required_integrity(input_root,files):
    """Verify selected archive inputs against RELATIVE repository provenance.

    Historical binary hashes/absolute workstation paths are never needed by
    plot-only mode. For a new campaign its own fraction/data gates apply.
    """
    manifest=ROOT/'repository_manifest.tsv';prefix='response/results/'+input_root.name+'/'
    expected={}
    if manifest.exists():
        with manifest.open(newline='') as f:
            expected={r['path'][len(prefix):]:r['sha256'] for r in csv.DictReader(f,delimiter='\t') if r['path'].startswith(prefix)}
    hashes={}
    for p in files:
        rel=str(p.relative_to(input_root));hashes[rel]=digest(p)
        if rel in expected and hashes[rel]!=expected[rel]:raise RuntimeError(f'Archived input hash mismatch: {rel}')
    if expected and any(rel not in expected for rel in hashes):raise RuntimeError('Archive input missing from the relative file manifest')
    return hashes,bool(expected)


def panel_axes(count,landscape=False):
    cols=min(5 if landscape else 4,count);rows=math.ceil(count/cols)
    if landscape and count==13:fig,axes=plt.subplots(3,5,figsize=(12.24,4.98),sharex=True,sharey=True,squeeze=False);fig.subplots_adjust(left=.052,right=.987,bottom=.10,top=.965,wspace=.16,hspace=.28)
    else:fig,axes=plt.subplots(rows,cols,figsize=(2.7*cols,2.05*rows+.9),sharex=True,sharey=True,squeeze=False);fig.subplots_adjust(left=.075,right=.985,bottom=.13,top=.89,wspace=.18,hspace=.25)
    for ax in list(axes.flat)[count:]:ax.set_visible(False)
    # Show the x ticks of every bottommost LIVE panel in each column.
    for col in range(cols):
        live=[axes[r,col] for r in range(rows) if r*cols+col<count]
        if live:live[-1].tick_params(labelbottom=True)
    return fig,axes


def spectral_panel(arrays,Ulist,n,c,out,landscape=False):
    fig,axes=panel_axes(len(Ulist),landscape)
    for U,ax in zip(Ulist,axes.flat):
        d=arrays[(U,n)];im=ax.pcolormesh(d['k'],d['omega'],d['A_per_spin'].T,shading='nearest',cmap='rainbow',vmin=0,vmax=1.5,rasterized=True)
        ax.set(xlim=(0,np.pi),ylim=(c['plot_y_min'],c['plot_y_max']),xticks=[0,np.pi/2,np.pi],xticklabels=['0','π/2','π'],yticks=[-5,0,5,10])
        ax.tick_params(labelsize=9);ax.axhline(0,color='white',ls='--',lw=.5)
        ax.text(.035,.94,f'U/t = {U:g}',transform=ax.transAxes,va='top',color='white',fontsize=10)
        actual=float(d['actual_density'])
        if abs(actual-n)>1e-10:ax.text(.965,.94,f'n̄ = {actual:.3f}',transform=ax.transAxes,va='top',ha='right',color='white',fontsize=8)
    fig.supylabel('ω/t',x=.01);fig.supxlabel('k',y=.015)
    if landscape and len(Ulist)==13:
        ca=fig.add_axes([.654,.22,.309,.032])
        cb=fig.colorbar(im,cax=ca,orientation='horizontal',extend='max',ticks=[0,.5,1,1.5]);cb.set_label('A(k,ω) per spin')
    else:
        fig.suptitle(f'{c["L"]}-MMHK, target n={n:g}, η/t={c["eta"]:g}',fontsize=12)
        ca=fig.add_axes([.29,.07,.48,.015]);fig.colorbar(im,cax=ca,orientation='horizontal',extend='max',label='A(k,ω) per spin')
    save(fig,out,f'spectral_{"slide" if landscape else "panel"}_n{tag(n)}')


def rebroaden_spectral(root,p,d,eta):
    if abs(float(d['eta'])-eta)<1e-12:return d
    plan=p/'spectral_plan.json'
    if not plan.exists():raise ValueError('This archived spectral array is eta=0.1. Use eta=0.1, or a fresh campaign with portable raw spectral plans for rebroadened spectra.')
    A=np.zeros_like(d['A_per_spin']);mu=float(d['mu']);omega=d['omega']
    entries=json.loads(plan.read_text())
    if not entries:
        A=eta/np.pi/((omega[None,:]+mu+2*np.cos(d['k'][:,None]))**2+eta**2)
    hashes=json.loads((p/'metadata.json').read_text())['source_hashes']
    for entry in entries:
        relative=Path(entry['file'])
        if relative.is_absolute() or '..' in relative.parts:raise RuntimeError('Unsafe raw-fraction plan path')
        file=root/relative
        if hashes.get(str(relative))!=digest(file):raise RuntimeError('Raw spectral fraction hash mismatch')
        fraction=read_fraction(file)
        A[entry['j']]-=entry['factor']*fraction_G(fraction,omega+mu+1j*eta,entry['particle']).imag/np.pi
    return dict(d,A_per_spin=A,eta=np.array(eta))


PEAK_KEYS=['first_omega','dominant_omega','first_velocity','dominant_velocity',
           'lowest_significant_pole','lowest_pole_residue','lowest_pole_fraction',
           'integrated_weight','first_moment','lowest_pole_velocity']
def safe_peak(E,W,eta,q):
    if len(E) and W.sum()>0:
        try:return peak_record(E,W,eta,q)
        except (RuntimeError,IndexError,ValueError):pass
    # Missing/absent peaks are null, never a fabricated zero velocity.
    return dict.fromkeys(PEAK_KEYS)


def plot_all(root,out,c,slides=True):
    root=Path(root).resolve();out=Path(out).resolve()
    if out==root or root in out.parents and 'response/results' in str(out):raise ValueError('Write plots outside the immutable historical archive')
    out.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    specroot=root/'spectral'
    if not specroot.exists():specroot=root/'single_particle_eta0p1_k101_depth600'
    files=[];cache={};arrays={};metadata={}
    for U in c['U_list']:
        for n in c['physical_densities']:
            p=root/f'U{U:g}'/f'n{tag(n)}';s=specroot/f'U{U:g}'/f'n{tag(n)}'
            files.extend([p/'metadata.json',p/'averaged_qmin_poles.npz',s/'metadata.json',s/'Akw.npz'])
            m=json.loads((p/'metadata.json').read_text());metadata[(U,n)]=m
            if m['L']!=c['L'] or m['twist_count']!=c['twists']:raise RuntimeError('Wrong input mixing size/twist count')
            if abs(m['actual_density_Ne_over_L']-n)>c['density_tolerance']+1e-10:raise RuntimeError('Input density outside declared tolerance')
            with np.load(p/'averaged_qmin_poles.npz',allow_pickle=False) as d:
                for ch in COLORS:
                    E,W=d[f'{ch}_energies'].copy(),d[f'{ch}_weights'].copy()
                    if not np.isfinite(E).all() or not np.isfinite(W).all() or (E<0).any() or (W<0).any():raise RuntimeError('Invalid neutral poles')
                    cache[(U,n,ch)]=(E,W)
            with np.load(s/'Akw.npz',allow_pickle=False) as d:a={k:d[k].copy() for k in d.files}
            np.testing.assert_allclose(a['k'],np.linspace(0,np.pi,c['k_points']),atol=1e-12)
            if a['A_per_spin'].shape!=(c['k_points'],len(a['omega'])) or not np.isfinite(a['A_per_spin']).all() or (a['A_per_spin']<0).any():raise RuntimeError('Invalid spectral array')
            np.testing.assert_allclose(a['mu'],m['mu_selection']['mu'],atol=1e-12)
            np.testing.assert_allclose(a['actual_density'],m['actual_density_Ne_over_L'],atol=1e-12)
            arrays[(U,n)]=rebroaden_spectral(root,s,a,c['eta'])
    input_hashes,archive_checked=required_integrity(root,files)
    signature=dict(parameters=c,input_sha256=input_hashes)
    manifest=out/'plot_inputs.json'
    if manifest.exists() and json.loads(manifest.read_text())!=signature:raise RuntimeError('Plot inputs/settings differ; use a new output directory')
    manifest.write_text(json.dumps(signature,indent=2)+'\n')
    Ugrid=np.unique(np.r_[np.arange(min(c['U_list']),max(c['U_list'])+c['ba_step']/2,c['ba_step']),c['U_list']])
    ba=[solve(float(U),float(n),c['ba_nodes']) for n in c['physical_densities'] for U in Ugrid]
    write_csv(out/'BA_velocities_fine_grid.csv',ba)
    rows=[];low=np.linspace(0,2,2001);full=np.linspace(0,65,6501);curves={};population=[];sensitivity=[]
    for U in c['U_list']:
      for n in c['physical_densities']:
        m=metadata[(U,n)];actual=m['actual_density_Ne_over_L'];nominal=solve(U,n,c['ba_nodes']);at_actual=solve(U,actual,c['ba_nodes'])
        row=dict(U_over_t=U,target_density=n,actual_density=actual,mu=m['mu_selection']['mu'],q_min=2*np.pi/c['L'],eta=c['eta'])
        for ch in COLORS:
            E,W=cache[(U,n,ch)]
            r=safe_peak(E,W,c['eta'],row['q_min'])
            row.update({f'{ch}_{k}':v for k,v in r.items()})
            row[f'{ch}_BA_velocity_nominal_n']=nominal[f'v_{ch}'];row[f'{ch}_BA_velocity_actual_n']=at_actual[f'v_{ch}']
            curves[(U,n,ch)]={'low':spectrum(low,E,W,c['eta']) if len(E) else np.zeros_like(low),
                             'full':spectrum(full,E,W,c['eta']) if len(E) else np.zeros_like(full)}
            for eta in sorted(set([.05,.1,.2,.3,c['eta']])):
                try:r0=peak_record(E,W,eta,row['q_min']) if len(E) else {}
                except RuntimeError:r0={}
                sensitivity.append(dict(U=U,density=n,channel=ch,eta=eta,first_omega=r0.get('first_omega'),dominant_omega=r0.get('dominant_omega')))
        rows.append(row)
        for t in m['twists']:population.append(dict(U=U,target_density=n,k=t['k'],N=t['N'],ground_dimension=t['ground_dimension']))
        np.savez_compressed(out/f'U{U:g}_n{tag(n)}_responses.npz',omega_low=low,omega_full=full,
                            charge_low=curves[(U,n,'charge')]['low'],spin_low=curves[(U,n,'spin')]['low'],
                            charge_full=curves[(U,n,'charge')]['full'],spin_full=curves[(U,n,'spin')]['full'],eta=c['eta'])
    write_csv(out/'velocity_comparison.csv',rows);write_csv(out/'broadening_sensitivity.csv',sensitivity);write_csv(out/'selected_ground_states.csv',population)
    ymax=max(.05,.05*np.ceil(1.12*max(x['low'].max() for x in curves.values())/.05))
    for n in c['physical_densities']:
      for mode in ['low','full']:
        fig,axes=panel_axes(len(c['U_list']))
        for U,ax in zip(c['U_list'],axes.flat):
            for ch in COLORS:ax.plot(low if mode=='low' else full,curves[(U,n,ch)][mode],color=COLORS[ch],ls=STYLES[ch],lw=1.5)
            ax.set(xlim=(0,2 if mode=='low' else 65),ylim=(0,ymax) if mode=='low' else (1e-6,max(ymax,1)),yscale='linear' if mode=='low' else 'log')
            ax.text(.965,.94,f'U/t = {U:g}',transform=ax.transAxes,va='top',ha='right')
            actual=metadata[(U,n)]['actual_density_Ne_over_L']
            if abs(actual-n)>1e-10:ax.text(.965,.75,f'n̄ = {actual:.3f}',transform=ax.transAxes,va='top',ha='right',fontsize=8)
            ax.tick_params(labelsize=8)
        fig.suptitle(f'{c["L"]}-MMHK, target n={n:g}, {c["twists"]} twists, η/t={c["eta"]:g}',fontsize=12)
        fig.supxlabel('ω/t',y=.02);fig.supylabel('S(q₁,ω)',x=.012)
        from matplotlib.lines import Line2D
        fig.legend(handles=[Line2D([],[],color=COLORS[ch],ls=STYLES[ch],label=ch if ch=='charge' else 'spin (2Sz)') for ch in COLORS],loc='upper center',bbox_to_anchor=(.5,.95),ncol=2,frameon=False)
        save(fig,out,f'response_panel_n{tag(n)}'+('_full_range' if mode=='full' else ''))
      fig,axes=panel_axes(len(c['U_list']),landscape=True)
      for U,ax in zip(c['U_list'],axes.flat):
          for ch in COLORS:ax.plot(low,curves[(U,n,ch)]['low'],color=COLORS[ch],ls=STYLES[ch],lw=1.3)
          ax.set(xlim=(0,2),ylim=(0,ymax),xticks=[0,1,2]);ax.tick_params(labelsize=8)
          ax.text(.95,.94,f'U/t = {U:g}',transform=ax.transAxes,ha='right',va='top',fontsize=9)
          actual=metadata[(U,n)]['actual_density_Ne_over_L']
          if abs(actual-n)>1e-10:ax.text(.95,.73,f'n̄ = {actual:.3f}',transform=ax.transAxes,ha='right',fontsize=8)
      fig.supxlabel('ω/t',y=.015);fig.supylabel('S(q₁,ω)',x=.01)
      if len(c['U_list'])==13:
          from matplotlib.lines import Line2D
          axes.flat[13].set_visible(True);axes.flat[13].axis('off')
          axes.flat[13].legend(handles=[Line2D([],[],color=COLORS[ch],ls=STYLES[ch],label=ch if ch=='charge' else 'spin (2Sz)') for ch in COLORS],loc='center',frameon=False)
      save(fig,out,f'response_slide_n{tag(n)}')
      spectral_panel(arrays,c['U_list'],n,c,out)
      spectral_panel(arrays,c['U_list'],n,c,out,landscape=True)
      for U in c['U_list']:
        d=arrays[(U,n)];fig,ax=plt.subplots(figsize=(4.8,3.6),layout='constrained')
        im=ax.pcolormesh(d['k'],d['omega'],d['A_per_spin'].T,shading='nearest',cmap='rainbow',vmin=0,vmax=1.5,rasterized=True)
        ax.set(xlim=(0,np.pi),ylim=(c['plot_y_min'],c['plot_y_max']),xlabel='k',ylabel='ω/t',xticks=[0,np.pi/2,np.pi],xticklabels=['0','π/2','π']);ax.axhline(0,color='white',ls='--',lw=.5)
        fig.colorbar(im,ax=ax,orientation='horizontal',extend='max',label='A(k,ω) per spin')
        save(fig,out,f'Akw_U{U:g}_n{tag(n)}')
    fig,axes=plt.subplots(1,len(c['physical_densities']),figsize=(3.8*len(c['physical_densities']),4.5),sharey=True,squeeze=False)
    for ax,n in zip(axes.flat,c['physical_densities']):
        b=[r for r in ba if r['density']==n];r=[r for r in rows if r['target_density']==n]
        for ch in COLORS:
            ax.plot([z['U_over_t'] for z in b],[z[f'v_{ch}'] for z in b],color=COLORS[ch],ls=STYLES[ch])
            ax.plot([z['U_over_t'] for z in r],[z[f'{ch}_dominant_velocity'] for z in r],ls='none',marker='o' if ch=='charge' else 's',color=COLORS[ch])
            changed=[z for z in r if abs(z['actual_density']-n)>1e-10]
            ax.plot([z['U_over_t'] for z in changed],[z[f'{ch}_BA_velocity_actual_n'] for z in changed],ls='none',marker='D',mfc='none',color=COLORS[ch])
        ax.set(title=f'n = {n:g}',xlabel='U/t',ylim=(0,2.55));ax.grid(alpha=.15)
    axes.flat[0].set_ylabel('v/(t a₀/ℏ)');fig.suptitle('Spin and charge velocities')
    fig.subplots_adjust(bottom=.24,top=.82);fig.text(.5,.06,'BA: lines. MMHK dominant peak/q: markers. Finite-q estimators, q₁ = 2π/L.\nDiamonds: BA at the actual density when it differs from the target.',ha='center',fontsize=9)
    fig.savefig(out/'spin_charge_velocity_vs_U.png',dpi=180,facecolor='white')
    fig.savefig(out/'spin_charge_velocity_vs_U.svg',facecolor='white')
    fig.texts[-1].set_visible(False);fig.subplots_adjust(bottom=.15)
    save(fig,out,'spin_charge_velocity_slide')
    density_rows=[];stairs={}
    for U in c['U_list']:
        p=root/f'U{U:g}'/'ground_energy_table.npz'
        if not p.exists():continue
        with np.load(p) as d:E,w=d['E0_by_N'],d['weights']
        ps=choose_mu(E,w,c['physical_densities'][0]/2,c['L'])['plateaus'];stairs[U]=ps
        for n in c['physical_densities']:
            s=choose_mu(E,w,n/2,c['L']);density_rows.append(dict(U=U,target_density=n,nearest_stable_density=2*s['filling'],exact_stable=s['exact_target_attained'],mu=s['mu'],mu_left=s['mu_left'],mu_right=s['mu_right']))
    write_csv(out/'density_attainability.csv',density_rows)
    if stairs:
        fig,axes=panel_axes(len(stairs))
        for (U,ps),ax in zip(stairs.items(),axes.flat):
            ax.stairs([2*p['filling'] for p in ps],[ps[0]['mu_left']]+[p['mu_right'] for p in ps],baseline=None)
            ax.set(title=f'U/t = {U:g}',ylim=(0,2));ax.tick_params(labelsize=8)
        fig.supxlabel('μ/t',y=.02);fig.supylabel('n = Ne/L',x=.01);fig.suptitle(f'{c["L"]}-MMHK: {c["twists"]}-twist stable density plateaus')
        save(fig,out,'density_vs_mu')
    if 7 in stairs:
        ps=stairs[7];fig,axes=plt.subplots(1,len(c['physical_densities']),figsize=(3.6*len(c['physical_densities']),3.7),squeeze=False)
        for ax,n in zip(axes.flat,c['physical_densities']):
            r=next(r for r in density_rows if r['U']==7 and r['target_density']==n);width=max(.005,r['mu_right']-r['mu_left'])
            ax.stairs([2*p['filling'] for p in ps],[ps[0]['mu_left']]+[p['mu_right'] for p in ps],baseline=None)
            ax.axhline(n,color='#777777',ls='--');ax.plot(r['mu'],r['nearest_stable_density'],'o')
            ax.set(title=f'Target n={n:g}',xlim=(r['mu_left']-.35*width,r['mu_right']+.35*width),ylim=(n-.018,n+.018),xlabel='μ/t',ylabel='n')
        fig.suptitle('U/t=7: stable density plateaus');fig.tight_layout();save(fig,out,'U7_density_vs_mu')
    presentation=dict(config=c,rows=rows,ba=ba,ymax=float(ymax),responses=[dict(U=U,n=n,channel=ch,x=low.tolist(),y=d['low'].tolist()) for (U,n,ch),d in curves.items()])
    (out/'presentation_data.json').write_text(json.dumps(presentation,allow_nan=False))
    report=dict(cases=len(rows),source_integrity_checked=archive_checked,input_sha256=input_hashes,
                many_body_solver_launched=False,BA_book_required=False,eta=c['eta'],k_points=c['k_points'],
                velocity_definition='dominant omega(q_min)/q_min; finite-q estimator',
                response_definition='coherent charge and longitudinal 2Sz, positive inelastic Lehmann spectrum',
                spectral_definition='per-spin -Im(Gup+Gdn)/(2*pi)',plot_y_limits=[c['plot_y_min'],c['plot_y_max']])
    (out/'plot_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
