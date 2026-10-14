#!/usr/bin/env python3
"""Regression-check a newly calculated campaign against the frozen study.

Different orbital/twist grids are not asserted equivalent. This audit is for
L10, 50 twists, densities .7/.8/.9, depths 300/600 and eta=.1, k=101.
The fresh workflow never imports old results as a solver input.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'response'))
from fast_lorentzian_spectrum import spectrum
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
p.add_argument('--couplings',nargs='+',type=int,default=list(range(13)))
a=p.parse_args();new=a.input.resolve();old=ROOT/'response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
c=json.loads((new/'config.json').read_text())
assert c['L']==10 and c['twists']==50 and c['k_points']==101 and c['eta']==.1
assert c['response_depth']==300 and c['spectral_depth']==600
rows=[];omega=np.linspace(0,5,5001)
for U in a.couplings:
    with np.load(new/f'U{U}/ground_energy_table.npz') as x,np.load(old/f'U{U}/ground_energy_table.npz') as y:
        error=float(abs(x['E0_by_N']-y['E0_by_N']).max())
        np.testing.assert_allclose(x['E0_by_N'],y['E0_by_N'],atol=1e-8,rtol=0)
    for n in c['physical_densities']:
        tag=str(n).replace('.','p');relative=f'U{U}/n{tag}'
        assert (new/relative/'COMPLETE').is_file() and (new/'spectral'/relative/'COMPLETE').is_file()
        m=json.loads((new/relative/'metadata.json').read_text());r=json.loads((old/relative/'metadata.json').read_text())
        np.testing.assert_allclose(m['actual_density_Ne_over_L'],r['actual_density_Ne_over_L'],atol=1e-10,rtol=0)
        np.testing.assert_allclose(m['mu_selection']['mu'],r['mu_selection']['mu'],atol=1e-8,rtol=0)
        errors={}
        with np.load(new/relative/'averaged_qmin_poles.npz') as x,np.load(old/relative/'averaged_qmin_poles.npz') as y:
            for ch in ['charge','spin']:
                fresh=spectrum(omega,x[ch+'_energies'],x[ch+'_weights'],.1)
                golden=spectrum(omega,y[ch+'_energies'],y[ch+'_weights'],.1)
                errors[ch]=float(abs(fresh-golden).max())
                assert errors[ch]/max(float(golden.max()),1e-12)<1e-6,(U,n,ch,errors[ch])
        with np.load(new/'spectral'/relative/'Akw.npz') as x,np.load(old/'single_particle_eta0p1_k101_depth600'/relative/'Akw.npz') as y:
            for key in ['k','omega','mu','actual_density']:np.testing.assert_allclose(x[key],y[key],atol=1e-8,rtol=0)
            spectral_error=float(abs(x['A_per_spin']-y['A_per_spin']).max())
            # Finite-depth Lanczos spectra are not bitwise backend/inversion invariant.
            # Existing depth-600/900 anchors and +/-phase duplicates already exhibit
            # small differences; gate against 1e-4 of the global spectral maximum,
            # not a relative error on a vanishing high-energy tail.
            normalized_error=spectral_error/max(float(y['A_per_spin'].max()),1e-12)
            assert normalized_error<1e-4,(U,n,spectral_error,normalized_error)
        rows.append(dict(U=U,density=n,actual_density=m['actual_density_Ne_over_L'],mu=m['mu_selection']['mu'],
                         maximum_ground_energy_error=error,maximum_neutral_spectrum_errors=errors,
                         maximum_single_particle_spectrum_error=spectral_error,maximum_single_particle_error_relative_to_global_peak=normalized_error))
report=dict(passed=True,neutral_regression_tolerance_relative_to_global_peak=1e-6,spectral_regression_tolerance_relative_to_global_peak=1e-4,
            tolerance_basis='finite-depth Lanczos; frozen archive documents depth-600/900 differences and 3.92e-5 absolute duplicate-inversion differences; tight independent ED and sum-rule gates are separate',cases=len(rows),couplings=a.couplings,backend=c['backend'],
            executable_sha256=c['executable_sha256'],source_sha256=c['source_sha256'],
            numerical_configuration={key:c[key] for key in ['L','twists','physical_densities','response_depth','spectral_depth','k_points','eta','density_tolerance']},
            tests=rows,scope='fresh common-mu calculations versus frozen numerical archive; no archive used during solve')
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
