#!/usr/bin/env python3
"""Audit completed G data and make 13-panel A(k,w) figures for each density."""
import csv,json,os
from pathlib import Path
import numpy as np
from replot_twist50_response import digest,tag,HERE,NOTE

ROOT=HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
OUT=ROOT/'single_particle_eta0p1_depth600'

def main():
    assert (OUT/'COMPLETE').exists(),'Spectral campaign is still incomplete'
    cfg=json.loads((OUT/'config.json').read_text());assert cfg['depth']==600 and cfg['eta']==.1
    assert digest(HERE/'run_single_particle_scan.py')==cfg['runner_sha256']
    assert digest(HERE.parent/'create_code/chain_Hubbard_tbc_spectral.cc')==cfg['source_sha256']
    assert digest(HERE/'bin/chain_Hubbard_tbc_spectral.x')==cfg['exe_sha256']
    assert all(digest(f)==sha for f,sha in cfg['ground_input_sha256'].items())
    all_raw={};rows=[];data={}
    for U in range(13):
      for n in [.7,.8,.9]:
        directory=OUT/f'U{U}'/f'n{tag(n)}';assert (directory/'COMPLETE').exists()
        m=json.loads((directory/'metadata.json').read_text())
        with np.load(directory/'Akw.npz') as d:data[(U,n)]={k:d[k].copy() for k in d.files}
        d=data[(U,n)];assert d['A_per_spin'].shape==(251,2201)
        assert np.isfinite(d['A_per_spin']).all() and d['A_per_spin'].min()>=0
        ground=json.loads((ROOT/f'U{U}'/f'n{tag(n)}'/'metadata.json').read_text())
        assert abs(float(d['mu'])-ground['mu_selection']['mu'])<1e-12
        assert abs(m['actual_n']-ground['actual_density_Ne_over_L'])<1e-12
        assert abs(2*m['pole_occupation_per_spin']-m['actual_n'])<1e-8
        all_raw.update(m['source_hashes'])
        rows.append(dict(U_over_t=U,target_n=n,actual_n=m['actual_n'],mu=m['mu'],eta=.1,
            momentum_count=251,frequency_count=2201,source_norm_error=m['source_norm_error'],
            source_first_moment_error=m['source_first_moment_error'],pole_density=2*m['pole_occupation_per_spin'],
            A_max=float(d['A_per_spin'].max()),colour_clip_fraction=float(np.mean(d['A_per_spin']>1.5))))
    assert all(digest(f)==sha for f,sha in all_raw.items())
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for n in [.7,.8,.9]:
        fig,axes=plt.subplots(4,4,figsize=(10.8,8),sharex=True,sharey=True)
        fig.subplots_adjust(left=.075,right=.985,bottom=.18,top=.9,hspace=.2,wspace=.16)
        for U,ax in zip(range(13),axes.flat):
            d=data[(U,n)];im=ax.imshow(d['A_per_spin'].T,origin='lower',extent=[0,np.pi,-6,16],
                aspect='auto',cmap='rainbow',vmin=0,vmax=1.5,interpolation='none')
            ax.axhline(0,color='white',lw=.5,ls='--')
            ax.set(xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'],yticks=[-5,0,5,10,15])
            ax.tick_params(labelsize=8)
            ax.text(.04,.95,rf'$U/t={U}$',transform=ax.transAxes,va='top',fontsize=9,color='white')
            if abs(float(d['actual_density'])-n)>1e-10:
                ax.text(.96,.95,rf'$\bar n={float(d["actual_density"]):.3f}$',transform=ax.transAxes,
                        va='top',ha='right',fontsize=8,color='white')
        for ax in list(axes.flat)[13:]:ax.set_visible(False)
        for ax in axes[2,1:]:ax.tick_params(labelbottom=True)
        fig.suptitle(rf'Single-particle $A(k,\omega)$: target $n={n:g}$, 10-MMHK, $\eta/t=0.1$',fontsize=13,y=.98)
        fig.supxlabel(r'$k$',y=.105);fig.supylabel(r'$\omega/t$',x=.015)
        cax=fig.add_axes([.23,.045,.54,.018]);fig.colorbar(im,cax=cax,orientation='horizontal',extend='max',label=r'$A(k,\omega)$ per spin')
        fig.savefig(OUT/f'spectral_panel_n{tag(n)}.png',dpi=180);fig.savefig(OUT/f'spectral_panel_n{tag(n)}.svg');plt.close(fig)
    with (OUT/'case_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    anchors=[HERE/'results/spectral_benchmark_L10_U8_N8_Q0_depth900/depth_validation.json',
             HERE/'results/spectral_benchmark_L10_U12_N7_nonzero_depth900/depth_validation.json']
    anchor_checks=[json.loads(p.read_text()) for p in anchors]
    v=dict(passed=True,cases=39,momenta_per_case=251,raw_fractions_verified=len(all_raw),
        maximum_source_norm_error=max(r['source_norm_error'] for r in rows),
        maximum_source_first_moment_error=max(r['source_first_moment_error'] for r in rows),
        maximum_pole_density_error=max(abs(r['pole_density']-r['actual_n']) for r in rows),
        independent_L4_ED=json.loads((HERE/'results/spectral_validation_L4_U4/independent_ED_validation.json').read_text()),
        depth600_vs900_anchors=anchor_checks,full_all_case_depth_convergence_claimed=False,
        input_hashes_verified=True,note_sha256=digest(NOTE))
    (OUT/'validation.json').write_text(json.dumps(v,indent=2)+'\n')
    (OUT/'README.txt').write_text(
        'Single-particle spectra for L=10 1D-MMHK, U/t=0..12, targets n=Ne/L=.7,.8,.9.\n'
        'Use the identical common mu and full spin/translation ground projectors of the50-twist response campaign.\n'
        'A(k,w)=-Im[Gup(k,w)+Gdn(k,w)]/(2*pi), per-spin convention, total spectral weight1.\n'
        'Particle frequencies Em(N+1)-E0-mu; hole frequencies E0-Em(N-1)-mu.\n'
        'Do not use the half-filled U/2 shift away from half filling.\n'
        'One physical k selects one MMHK block/mode. Average ground degeneracy and spin at fixed k,\n'
        'rather than averaging distinct physical momenta into one spectrum.\n'
        'Physical momentum p=K-phi=2*pi*q/L+k_RBZ. This equals the external-twist q=0 measurement\n'
        'after the exact discrete large-gauge transformation phi->phi-K.\n'
        'All251 unique nonnegative points of the500-site momentum grid, from0 to pi, are included.\n'
        'Inversion symmetry reduces50 boundary phases to26 representatives. Validate reflected\n'
        'ground populations and energies; zero/boundary duplicate spectral values must agree.\n'
        'U0 uses the exact free Green function, independent of the degenerate free occupation trace.\n'
        'Positive U uses new qbasis charged Lanczos sources, depth600, OMP/MKL1, eight concurrent jobs.\n'
        'Eta/t=.1. Plot window omega/t=[-6,16], step.01. All raw fractions remain available.\n'
        'The rainbow palette matches the half-filled notebook style. Colour range0..1.5 per spin,\n'
        'with an upper saturation arrow; NPZ intensities are unmodified. The old half-filled\n'
        'notebooks plot -2ImGup, which is2*pi times this per-spin A when spin symmetry holds.\n'
        'All density exceptions are preserved: n0.9 uses.896 at U7,.902 at U9,.898 at U10.\n'
        'Each U/density folder contains Akw.npz, PNG/SVG, source hashes, mu and validation metadata.\n'
        'spectral_panel_n*: all13U panels per density. case_summary.csv and validation.json provide audits.\n'
        'Independent L4 full-Fock ED tests all momenta and both particle/hole branches.\n'
        '600/900-depth anchors at U8,N8,zero phase and U12,N7,nonzero phase test both channels.\n'
        'These selected convergence tests do not constitute a global all-case error bound.\n'
        'The student solver and open LaTeX note were not modified or executed.\n')
    print('FINAL_SPECTRAL_AUDIT_COMPLETE',OUT,flush=True)

if __name__=='__main__':main()
