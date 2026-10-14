#!/usr/bin/env python3
"""Rebroaden existing poles and re-extract peaks without a neutral solver run."""
import argparse, csv, json, os
from pathlib import Path
import numpy as np
from replot_twist50_response import digest, tag, write_csv, velocity_plot, response_plots, HERE, NOTE
from postprocess_twist50_density_scan import peak_record
from fast_lorentzian_spectrum import spectrum

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--eta',type=float,default=.1)
    a=p.parse_args()
    if a.eta<=0:p.error('eta must be positive')
    root=HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'; old=root/'analysis'
    fine=root/'analysis_ba_step0p2';out=root/f'analysis_ba_step0p2_eta{tag(a.eta)}';out.mkdir(exist_ok=True)
    config=json.loads((root/'config.json').read_text());m=json.loads((old/'metadata.json').read_text())
    before=digest(NOTE); inputs={}
    with (old/'velocity_comparison.csv').open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    with (fine/'BA_velocities_fine_grid.csv').open() as f:ba=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    assert len(rows)==39 and len(ba)==183
    spectra={}; audit=[]
    for r in rows:
        U,n=int(r['U_over_t']),r['target_density'];path=root/f'U{U}'/f'n{tag(n)}'/'averaged_qmin_poles.npz'
        inputs[str(path)]=digest(path);assert inputs[str(path)]==m['input_sha256'][str(path)]
        r['eta']=a.eta
        with np.load(path) as d:
            for ch in ['charge','spin']:
                E,W=d[ch+'_energies'],d[ch+'_weights'];prior=r[ch+'_dominant_velocity']
                rec=peak_record(E,W,a.eta,r['q_min']);r.update({f'{ch}_{k}':v for k,v in rec.items()})
                r[ch+'_relative_error_actual_BA']=rec['dominant_velocity']/r[ch+'_BA_velocity_actual_n']-1
                lo=np.linspace(0,2,2001);full=np.linspace(0,65,6501)
                s=dict(omega_low=lo,S_low=spectrum(lo,E,W,a.eta),omega_full=full,S_full=spectrum(full,E,W,a.eta),eta=np.array(a.eta))
                spectra[(U,n,ch)]=s;np.savez_compressed(out/f'U{U}_n{tag(n)}_{ch}_spectra.npz',**s)
                audit.append(dict(U=U,target_density=n,channel=ch,eta=a.eta,previous_eta0p2_velocity=prior,**rec))
        print(f'U={U}, n={n:g}: eta={a.eta:g} peaks updated',flush=True)
    write_csv(out/'velocity_comparison.csv',rows);write_csv(out/'peak_audit.csv',audit)
    write_csv(out/'BA_velocities_fine_grid.csv',ba)
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    velocity_plot(plt,rows,ba,config['physical_densities'],out,.2,eta=a.eta)
    panels=response_plots(plt,rows,spectra,config['physical_densities'],config['U_list'],out,eta=a.eta)
    panels+=response_plots(plt,rows,spectra,config['physical_densities'],config['U_list'],out,full=True,eta=a.eta)
    # Complete literal values for editable PowerPoint charts, without decimation.
    literal=dict(eta=a.eta,rows=rows,ba=ba,ymax=max(s['S_low'].max() for s in spectra.values())*1.12,
        responses=[dict(U=U,n=n,channel=ch,x=s['omega_low'].tolist(),y=s['S_low'].tolist()) for (U,n,ch),s in spectra.items()])
    (out/'presentation_data.json').write_text(json.dumps(literal))
    assert digest(NOTE)==before and all(digest(f)==sha for f,sha in inputs.items())
    (out/'metadata.json').write_text(json.dumps(dict(eta=a.eta,q_min=np.pi/5,L=10,twists=50,
        response_source='original averaged poles; all 39 cases rebroadened; peak positions re-extracted',
        BA_source=str(fine/'BA_velocities_fine_grid.csv'),BA_sha256=digest(fine/'BA_velocities_fine_grid.csv'),
        input_sha256=inputs,source_sha256=digest(__file__),open_note_unchanged=True,neutral_solver_executed=False,
        first_dominant_differ_count=sum(abs(r[ch+'_first_omega']-r[ch+'_dominant_omega'])>1e-5 for r in rows for ch in ['charge','spin']),
        panels=panels),indent=2)+'\n')
    print('OUTPUT:',out,flush=True)
if __name__=='__main__':main()
