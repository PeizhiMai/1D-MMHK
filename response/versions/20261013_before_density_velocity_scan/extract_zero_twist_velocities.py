#!/usr/bin/env python3
"""Finite-q peak-based velocity estimates from existing zero-twist L=10 data.

No Hamiltonian/ED solver is run. First and strongest low-energy peaks are
reported separately; neither is claimed to be a thermodynamic velocity.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.signal import find_peaks

from postprocess_response import read_fraction, poles_from_fraction, lorentzian

HERE=Path(__file__).resolve().parent
DEFAULT_INPUT=HERE/'results/zero_twist_benchmark'
DEFAULT_BA=HERE/'essler_bethe_ansatz/results_n0p9/velocities.csv'


def merge_poles(energies,weights,tolerance=1e-6):
    """Combine almost identical ghost/degenerate Ritz poles by summed residue."""
    order=np.argsort(energies);energies=np.asarray(energies)[order];weights=np.asarray(weights)[order]
    E=[];W=[];start=0
    while start<len(energies):
        stop=start+1
        while stop<len(energies) and energies[stop]-energies[start]<=tolerance:stop+=1
        total=float(weights[start:stop].sum())
        if total>0:
            E.append(float(energies[start:stop]@weights[start:stop]/total));W.append(total)
        start=stop
    return np.asarray(E),np.asarray(W)


def load_channel(directory,channel,q,provenance):
    m=json.loads((directory/'metadata.json').read_text())
    if not (directory/'COMPLETE').is_file():raise ValueError(f'Incomplete run {directory}')
    if m['L']!=10 or m['Nup']+m['Ndn']!=9 or m['t']!=1 or m['phase_per_bond']!=0:
        raise ValueError('Require zero-twist L=10, Ne=9, t=1')
    if m.get('observable','density')!='density':raise ValueError('Wrong observable')
    energies=[];weights=[];components=[]
    seen=set()
    for path in sorted(directory.glob(f'{channel}_gsQ*_q*.dat')):
        d=read_fraction(path)
        if int(d['q_index'])!=q:continue
        identity=(int(d['ground_Q']),int(d.get('ground_id',0)))
        if identity in seen:raise ValueError('Duplicate ground component')
        seen.add(identity)
        e,w=poles_from_fraction(d);w*=d['ensemble_weight']
        mask=(e>1e-8)&(w>0);e=e[mask];w=w[mask]
        energies.extend(e);weights.extend(w)
        components.append((identity,*merge_poles(e,w)))
        provenance[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    required={(g['Q'],g.get('id',0)) for g in m['ground_states']}
    if seen!=required:raise ValueError('Missing ground-state source')
    rawE=np.asarray(energies);rawW=np.asarray(weights)
    e,w=merge_poles(rawE,rawW)
    np.testing.assert_allclose(w.sum(),rawW.sum(),atol=1e-12,rtol=1e-12)
    np.testing.assert_allclose(e@w,rawE@rawW,atol=1e-11,rtol=1e-11)
    return e,w,components,m


def refined_peaks(energies,weights,eta,cutoff=5.,prominence_fraction=.01):
    omega=np.linspace(0,cutoff,int(cutoff/.0005)+1)
    spectrum=lorentzian(omega,energies,weights,eta)
    indices,_=find_peaks(spectrum,height=.001*spectrum.max(),
                        prominence=prominence_fraction*spectrum.max())
    peaks=[]
    for i in indices:
        def slope(x):
            z=x-energies
            return float(np.sum(-2*weights*eta*z/(np.pi*(z*z+eta*eta)**2)))
        lo,hi=omega[i-1],omega[i+1]
        if slope(lo)>0 and slope(hi)<0:
            peak=brentq(slope,lo,hi,xtol=1e-12)
        else:
            peak=minimize_scalar(lambda x:-float(lorentzian(np.array([x]),energies,weights,eta)[0]),
                                 bounds=(lo,hi),method='bounded',options={'xatol':1e-11}).x
        height=float(lorentzian(np.array([peak]),energies,weights,eta)[0])
        peaks.append({'omega':float(peak),'height':height})
    return peaks


def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def make_plots(rows,ba,cache,output,eta):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11})
    fig,axes=plt.subplots(1,2,figsize=(9.1,4.3),sharey=True)
    fig.subplots_adjust(left=.075,right=.985,bottom=.29,top=.79,wspace=.10)
    U=np.array([r['U_over_t'] for r in rows])
    for ax,channel,color in zip(axes,['charge','spin'],['#c43c64','#234f91']):
        ax.plot([r['U_over_t'] for r in ba if r['U_over_t']<=12],
                [r[f'v_{channel}'] for r in ba if r['U_over_t']<=12],color=color,lw=2,label='Thermodynamic Bethe ansatz')
        ax.plot(U,[r[f'{channel}_first_velocity'] for r in rows],'o',color='#222222',ms=5,label=r'$L=10$: first low-energy peak')
        if channel=='charge':
            ax.plot(U,[r[f'{channel}_dominant_velocity'] for r in rows],'D',mfc='none',mec='#c47616',ms=6,
                    label=r'$L=10$: strongest low-energy peak')
        ax.set(title=channel.capitalize(),xlabel=r'$U/t$',xlim=(0,12.5),ylim=(0,2.3))
        ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)
        ax.legend(fontsize=8,loc='upper right',frameon=False)
    axes[0].set_ylabel(r'$v/(t a_0/\hbar)$')
    fig.suptitle(r'Zero-twist Hubbard ring: $L=10$, $N_e=9$, $n=0.9$',fontsize=13,y=.97)
    fig.text(.5,.055,rf'Finite-q estimator: $\widehat v=\omega_{{\rm peak}}(q_1)/q_1$,  $q_1=\pi/5$,  $\eta/t={eta:g}$'
             '\nThese points are not a fitted '+r'$q\to0$'+' velocity; no finite-size extrapolation.',ha='center',fontsize=9)
    fig.savefig(output/'zero_twist_L10_velocities_vs_bethe_ansatz.png',dpi=180)
    fig.savefig(output/'zero_twist_L10_velocities_vs_bethe_ansatz.svg')
    plt.close(fig)
    fig,axes=plt.subplots(len(rows),2,figsize=(8.5,2.05*len(rows)),layout='constrained',squeeze=False)
    omega=np.linspace(0,2.0,4001);q1=2*np.pi/10
    for i,r in enumerate(rows):
        for j,channel,color in zip([0,1],['charge','spin'],['#c43c64','#234f91']):
            E,W=cache[(r['U_over_t'],channel,1)];s=lorentzian(omega,E,W,eta)
            ax=axes[i,j];ax.plot(omega,s/s.max(),color=color,lw=1.25)
            first=r[f'{channel}_first_omega'];strong=r[f'{channel}_dominant_omega']
            ax.axvline(first,color='#222222',ls=':',lw=1,label='First peak')
            if channel=='charge':ax.axvline(strong,color='#c47616',ls='--',lw=1,label='Strongest peak')
            ax.axvline(q1*r[f'{channel}_BA_velocity'],color='#888888',ls='-.',lw=1,label=r'$q_1 v_{\rm BA}$ (linear guide)')
            ax.set(title=f"U/t={r['U_over_t']:g}, {channel}",xlim=(0,2),ylim=(0,1.1))
            if j==0:ax.set_ylabel('Normalized intensity')
            if i==0:ax.legend(fontsize=7,loc='upper right')
            if i==len(rows)-1:ax.set_xlabel(r'$\omega/t$')
    fig.suptitle(rf'Peak-selection audit at $q_1=\pi/5$: $L=10$, $N_e=9$, $\theta=0$, $\eta/t={eta:g}$',fontsize=12)
    fig.savefig(output/'qmin_peak_selection_audit.png',dpi=170);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=DEFAULT_INPUT)
    p.add_argument('--bethe-ansatz',type=Path,default=DEFAULT_BA)
    p.add_argument('--output',type=Path,default=HERE/'results/zero_twist_velocity_comparison_n0p9')
    p.add_argument('--eta',type=float,default=.2)
    p.add_argument('--peak-prominence',type=float,default=.01)
    args=p.parse_args()
    if args.eta<=0 or args.peak_prominence<0:p.error('Invalid broadening or peak prominence')
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    ba=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(args.bethe_ansatz.open())]
    if any(abs(r['density']-.9)>1e-12 for r in ba):raise ValueError('Wrong BA density')
    BA={r['U_over_t']:r for r in ba}
    q1=2*np.pi/10;rows=[];sensitivity=[];pole_audit=[];momentum_audit=[];component_audit=[]
    cache={};provenance={}
    for U in [4.,6.,7.,8.,12.]:
        directory=args.input/f'L10_U{U:g}_n0p45_phi0'
        row={'U_over_t':U,'density_Ne_over_L':.9,'q_min':q1,'eta':args.eta}
        for channel in ['charge','spin']:
            for q in [1,2,3]:
                E,W,components,m=load_channel(directory,channel,q,provenance)
                if m['U']!=U:raise ValueError('Run coupling does not match its directory label')
                if {g['Q'] for g in m['ground_states']}!={2,8}:raise ValueError('Unexpected ground momenta for this dataset')
                cache[(U,channel,q)]=(E,W)
                peaks=refined_peaks(E,W,args.eta,prominence_fraction=args.peak_prominence)
                if not peaks:raise ValueError(f'No resolved peak U={U} {channel} q={q}')
                first=peaks[0];dominant=max(peaks,key=lambda x:x['height']);qphysical=q*q1
                momentum_audit.append({'U_over_t':U,'channel':channel,'q_index':q,'q':qphysical,
                                       'first_omega':first['omega'],'first_omega_over_q':first['omega']/qphysical,
                                       'dominant_omega':dominant['omega'],'dominant_omega_over_q':dominant['omega']/qphysical})
                if q!=1:continue
                row.update({f'{channel}_first_omega':first['omega'],f'{channel}_dominant_omega':dominant['omega'],
                            f'{channel}_first_velocity':first['omega']/q1,
                            f'{channel}_dominant_velocity':dominant['omega']/q1,
                            f'{channel}_first_to_dominant_height':first['height']/dominant['height'],
                            f'{channel}_BA_velocity':BA[U][f'v_{channel}']})
                row[f'{channel}_first_relative_error']=(row[f'{channel}_first_velocity']/row[f'{channel}_BA_velocity']-1)
                row[f'{channel}_dominant_relative_error']=(row[f'{channel}_dominant_velocity']/row[f'{channel}_BA_velocity']-1)
                for eta in [.05,.1,.2,.3]:
                    found=refined_peaks(E,W,eta,prominence_fraction=args.peak_prominence)
                    if not found:continue
                    strongest=max(found,key=lambda x:x['height'])
                    sensitivity.append({'U_over_t':U,'channel':channel,'eta':eta,'first_omega':found[0]['omega'],
                                        'first_velocity':found[0]['omega']/q1,
                                        'dominant_omega':strongest['omega'],'dominant_velocity':strongest['omega']/q1})
                for threshold in [1e-4,1e-3,1e-2]:
                    keep=W>threshold*W.sum();lowest=float(E[keep].min())
                    pole_audit.append({'U_over_t':U,'channel':channel,'fractional_weight_threshold':threshold,
                                       'lowest_omega':lowest,'lowest_omega_over_q':lowest/q1})
                for identity,ce,cw in components:
                    found=refined_peaks(ce,cw,args.eta,prominence_fraction=args.peak_prominence)
                    if not found:continue
                    strongest=max(found,key=lambda x:x['height'])
                    component_audit.append({'U_over_t':U,'channel':channel,'ground_Q':identity[0],
                                            'ground_id':identity[1],'first_omega':found[0]['omega'],
                                            'dominant_omega':strongest['omega']})
        rows.append(row)
        print(f"U={U:g}: charge(first,strong,BA)={row['charge_first_velocity']:.6f}, {row['charge_dominant_velocity']:.6f}, {row['charge_BA_velocity']:.6f}; spin(first,BA)={row['spin_first_velocity']:.6f}, {row['spin_BA_velocity']:.6f}",flush=True)
    write_csv(output/'velocity_comparison.csv',rows)
    write_csv(output/'broadening_sensitivity.csv',sensitivity)
    write_csv(output/'pole_threshold_audit.csv',pole_audit)
    write_csv(output/'momentum_dependence.csv',momentum_audit)
    write_csv(output/'ground_component_peaks.csv',component_audit)
    make_plots(rows,ba,cache,output,args.eta)
    kF=np.pi*.9/2
    metadata={'L':10,'Nup':5,'Ndn':4,'density_Ne_over_L':.9,'earlier_filling_x':.45,
              'phase_per_bond':0,'ground_momenta_averaged':[2,8],
              'q_min':q1,'k_F_thermodynamic':kF,'folded_4kF_magnitude':2*np.pi-4*kF,
              'qmin_equals_folded_4kF':bool(np.isclose(q1,2*np.pi-4*kF)),
              'estimator':'first significant local maximum / q_min; strongest local maximum reported separately',
              'search_window_omega':[0,5],'eta':args.eta,'prominence_fraction':args.peak_prominence,
              'pole_merge_tolerance':1e-6,'no_new_solver_runs':True,
              'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'raw_fraction_sha256':provenance,
              'bethe_ansatz_csv':str(args.bethe_ansatz.resolve()),
              'bethe_ansatz_sha256':hashlib.sha256(args.bethe_ansatz.read_bytes()).hexdigest(),
              'qualification':'Finite-q peak estimates, not a controlled thermodynamic velocity extraction. One hole, odd-electron ground doublet, coarse q spacing and a folded 4kF charge harmonic.'}
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Saved comparison:',output)


if __name__=='__main__':main()
