#!/usr/bin/env python3
"""Validate direct101-point spectra and assemble all13U panels per density."""
import csv,json,os
from pathlib import Path
import numpy as np
from replot_twist50_response import digest,tag,HERE,NOTE
ROOT=HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
OUT=ROOT/'single_particle_eta0p1_k101_depth600'
Y_MIN=-6
Y_MAX=10

def main():
    assert (OUT/'COMPLETE').exists(),'The101-point calculation is still incomplete'
    cfg=json.loads((OUT/'config.json').read_text());assert cfg['k_count']==101 and not cfg['interpolate']
    assert all(digest(p)==sha for p,sha in cfg['source_sha256'].items())
    before=json.loads((ROOT/'single_particle_eta0p1_depth600/config.json').read_text())
    assert all(digest(p)==sha for p,sha in before['ground_input_sha256'].items())
    original_analysis=json.loads((ROOT/'analysis/metadata.json').read_text())
    assert all(digest(p)==sha for p,sha in original_analysis['input_sha256'].items())
    arrays={};records=[];raw={}
    for U in range(13):
      for n in [.7,.8,.9]:
        p=OUT/f'U{U}'/f'n{tag(n)}';assert (p/'COMPLETE').exists()
        m=json.loads((p/'metadata.json').read_text());assert m['k_points']==101 and not m['interpolated']
        ref=json.loads((ROOT/f'U{U}'/f'n{tag(n)}'/'metadata.json').read_text())
        assert abs(m['mu']-ref['mu_selection']['mu'])<1e-12
        assert abs(m['reference_n_50_twists']-ref['actual_density_Ne_over_L'])<1e-12
        with np.load(p/'Akw.npz') as d:arrays[(U,n)]={k:d[k].copy() for k in d.files}
        d=arrays[(U,n)];np.testing.assert_array_equal(d['k'],np.linspace(0,np.pi,101))
        assert d['A_per_spin'].shape==(101,2201) and np.isfinite(d['A_per_spin']).all() and d['A_per_spin'].min()>=0
        raw.update(m['source_hashes'])
        records.append(dict(U=U,target_n=n,reference_n_50_twists=m['reference_n_50_twists'],
            sampled_n_20_folded_phases=m['sample_n_20_folded_phases'],mu=m['mu'],eta=.1,k_points=101,
            source_norm_error=m['source_norm_error'],source_first_moment_error=m['source_first_moment_error'],
            duplicate_inversion_max_abs_error=m['duplicate_inversion_max_abs_error'],
            A_max=float(d['A_per_spin'].max()),fraction_of_colour_bins_saturated=float(np.mean(d['A_per_spin']>1.5))))
    assert all(digest(p)==sha for p,sha in raw.items())
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for n in [.7,.8,.9]:
        fig,axes=plt.subplots(4,4,figsize=(10.8,8),sharex=True,sharey=True)
        fig.subplots_adjust(left=.075,right=.985,bottom=.19,top=.9,hspace=.2,wspace=.16)
        for U,ax in zip(range(13),axes.flat):
            d=arrays[(U,n)];im=ax.pcolormesh(d['k'],d['omega'],d['A_per_spin'].T,
                shading='nearest',cmap='rainbow',vmin=0,vmax=1.5,rasterized=True)
            ax.set_xlim(0,np.pi);ax.set_ylim(Y_MIN,Y_MAX)
            ax.axhline(0,color='white',lw=.5,ls='--')
            ax.set(xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'],yticks=[-5,0,5,10]);ax.tick_params(labelsize=8)
            ax.text(.04,.95,rf'$U/t={U}$',transform=ax.transAxes,va='top',fontsize=9,color='white')
            nr=float(d['actual_density'])
            if abs(nr-n)>1e-10:ax.text(.96,.95,rf'$n_{{\rm ref}}={nr:.3f}$',transform=ax.transAxes,va='top',ha='right',fontsize=8,color='white')
        for ax in list(axes.flat)[13:]:ax.set_visible(False)
        for ax in axes[2,1:]:ax.tick_params(labelbottom=True)
        fig.suptitle(rf'Single-particle $A(k,\omega)$: target $n={n:g}$, 10-MMHK, $\eta/t=0.1$',fontsize=13,y=.98)
        fig.supxlabel(r'$k$',y=.12);fig.supylabel(r'$\omega/t$',x=.015)
        cax=fig.add_axes([.23,.085,.54,.018]);fig.colorbar(im,cax=cax,orientation='horizontal',extend='max',label=r'$A(k,\omega)$ per spin')
        fig.text(.5,.012,'101 directly calculated momentum points. Chemical potential retains the 50-twist density constraint.',ha='center',fontsize=8)
        fig.savefig(OUT/f'spectral_panel_n{tag(n)}.png',dpi=180);fig.savefig(OUT/f'spectral_panel_n{tag(n)}.svg');plt.close(fig)
    # Keep discrete sample centres at their actual coordinates in individual
    # maps too. The edge half-cells are clipped to the requested 0..pi domain.
    for (U,n),d in arrays.items():
        p=OUT/f'U{U}'/f'n{tag(n)}';fig,ax=plt.subplots(figsize=(4.8,3.6),layout='constrained')
        im=ax.pcolormesh(d['k'],d['omega'],d['A_per_spin'].T,shading='nearest',
                        cmap='rainbow',vmin=0,vmax=1.5,rasterized=True)
        ax.set(xlim=(0,np.pi),ylim=(Y_MIN,Y_MAX),xlabel=r'$k$',ylabel=r'$\omega/t$',
               xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'],yticks=[-5,0,5,10])
        ax.axhline(0,color='white',ls='--',lw=.6)
        ax.text(.03,.96,rf'$U/t={U},\ n_{{\rm ref}}={float(d["actual_density"]):.3f}$',
                transform=ax.transAxes,va='top',color='white',fontsize=9)
        fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.07,pad=.15,extend='max',label=r'$A(k,\omega)$ per spin')
        fig.savefig(p/'spectral_Akw.png',dpi=180,bbox_inches='tight');fig.savefig(p/'spectral_Akw.svg',bbox_inches='tight');plt.close(fig)
    with (OUT/'case_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    validation=dict(passed=True,cases=39,momentum_points=101,interpolated=False,
        plot_y_limits_omega_over_t=[Y_MIN,Y_MAX],
        raw_fraction_files_rehashed=len(raw),source_files_verified=True,
        maximum_source_norm_error=max(r['source_norm_error'] for r in records),
        maximum_source_first_moment_error=max(r['source_first_moment_error'] for r in records),
        maximum_duplicate_inversion_error=max(r['duplicate_inversion_max_abs_error'] for r in records),
        independent_L4_ED=json.loads((HERE/'results/spectral_validation_L4_U4/independent_ED_validation.json').read_text()),
        depth_anchors=[json.loads((HERE/'results'/p/'depth_validation.json').read_text()) for p in
            ['spectral_benchmark_L10_U8_N8_Q0_depth900','spectral_benchmark_L10_U12_N7_nonzero_depth900',
             'spectral_benchmark_L10_U12_N9_nonzero_depth900']],
        all_case_depth_convergence_claimed=False,open_note_sha256=digest(NOTE))
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    (OUT/'README.txt').write_text(
        'Direct single-particle A(k,w) for L10, U/t0..12, targets physical n=Ne/L=.7,.8,.9.\n'
        'Exactly101 uniform physical momenta k=j*pi/100, j=0..100. No interpolation from the251-point run.\n'
        'Eta/t=.1; stored omega/t=-6..16 with step.01; Lanczos depth600.\n'
        f'Spectral colormaps display omega/t={Y_MIN}..{Y_MAX}; stored NPZ data are unchanged.\n'
        'Per-spin A=-Im(Gup+Gdn)/(2*pi) has unit total spectral weight.\n'
        'Addition poles Em(N+1)-E0-mu; removal poles E0-Em(N-1)-mu.\n'
        'Reference density and common mu remain those of the original50-twist response campaign.\n'
        'The101 k grid selects20 folded phases, with11 inversion representatives. It changes\n'
        'spectral sampling only. Never retune mu using this coarser set of phases.\n'
        'case_summary.csv reports both reference_n_50_twists and sampled_n_20_folded_phases,\n'
        'which can differ. Figure n_ref labels refer to the50-twist density constraint.\n'
        'At off-grid phases, solve actual grand ground states at the same mu. Do not impose N or Sz.\n'
        'The full number-sector scan is rigorously pruned using the hopping norm bound\n'
        '||H(phi+d)-H(phi)||<=4Lt|sin(d/2)|. Solve candidate sectors over ALL translation Q,\n'
        'check the lower bounds of excluded N, then trace every ground spin multiplet/projector.\n'
        'Shared phases reuse only hash-matched native600-step charged sources.\n'
        'Physical momentum K-phi=2*pi*q/L+k_RBZ equals the external q=0 twist measurement\n'
        'after a discrete large-gauge transformation. No averaging of distinct physical k.\n'
        'U0 uses the exact free retarded Green function, independent of free occupation degeneracy.\n'
        'Rainbow colour range0..1.5, with an upper saturation arrow. NPZ values are unmodified.\n'
        'Legacy half-filled notebooks plot -2ImGup=2*pi*A_per_spin under spin symmetry.\n'
        'spectral_panel_n*:13 U panels per density. U*/n*: individual PNG/SVG, Akw.npz and metadata.\n'
        'Selected600/900 depth checks and independent L4 full-Fock ED support the implementation.\n'
        'Selected depth tests do not establish a global all-case numerical error bound.\n'
        'Spin/charge response plots and the four-slide PPT retain50 twists and eta.1.\n'
        'The student solver and open LaTeX note remain unchanged.\n')
    print('FINAL_101_AUDIT_COMPLETE',OUT,flush=True)
if __name__=='__main__':main()
