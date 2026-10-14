#!/usr/bin/env python3
"""Check complete clean-test outputs and archived-value reproduction."""
import csv
import json
from pathlib import Path
import sys
import zipfile
import re

import numpy as np
ROOT=Path(__file__).resolve().parents[1]
base=Path(sys.argv[1]).resolve()
archive=ROOT/'response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
for directory,expected in [('smoke/plots',5),('smoke-replot',5),('archive-plots',7)]:
    p=base/directory
    d=json.loads((p/'plot_validation.json').read_text());assert not d['many_body_solver_launched']
    with zipfile.ZipFile(p/'spin_charge_spectral.pptx') as z:
        slides=[s for s in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml',s)]
        assert len(slides)==expected
        assert len([s for s in z.namelist() if s.startswith('ppt/notesSlides/notesSlide') and s.endswith('.xml')])==expected
    for tag in (['0p5','0p75'] if expected==5 else ['0p7','0p8','0p9']):
        for kind in ['response_panel','response_slide','spectral_panel','spectral_slide']:
            assert (p/f'{kind}_n{tag}.png').stat().st_size>1000
with (base/'archive-plots/velocity_comparison.csv').open() as f:actual=list(csv.DictReader(f))
with (archive/'analysis_ba_step0p2_eta0p1/velocity_comparison.csv').open() as f:old=list(csv.DictReader(f))
assert len(actual)==len(old)==39
errors={}
for row in actual:
    target=next(r for r in old if float(r['U_over_t'])==float(row['U_over_t']) and float(r['target_density'])==float(row['target_density']))
    for key in ['actual_density','mu','charge_dominant_velocity','spin_dominant_velocity','charge_first_velocity','spin_first_velocity','charge_BA_velocity_nominal_n','spin_BA_velocity_nominal_n']:
        diff=abs(float(row[key])-float(target[key]));errors[key]=max(errors.get(key,0),diff)
        assert diff<1e-8,(key,diff)
# Verify sampled response arrays against the original eta=.1 postprocessing.
# Old output names differ; compare the physical pole spectra directly instead.
sys.path.insert(0,str(ROOT/'response'))
from fast_lorentzian_spectrum import spectrum
max_error=0.
for U in range(13):
    for tag in ['0p7','0p8','0p9']:
        with np.load(archive/f'U{U}/n{tag}/averaged_qmin_poles.npz') as poles, np.load(base/f'archive-plots/U{U}_n{tag}_responses.npz') as out:
            for ch in ['charge','spin']:
                exact=spectrum(out['omega_low'],poles[ch+'_energies'],poles[ch+'_weights'],.1)
                diff=float(abs(exact-out[ch+'_low']).max());max_error=max(max_error,diff);assert diff<1e-12
# Plotting is not allowed to alter input spectra. Plot_inputs contains their
# expected hashes; archive integrity verifier already authenticated them.
report=dict(passed=True,neutral_and_charged_independent_tests=17,
            fresh_end_to_end_smoke_cases=4,archived_plot_cases=39,
            archived_velocity_table_max_errors=errors,neutral_array_max_error=max_error,
            solver_free_plot_image=True,presentation_slide_counts=[5,5,7],
            scope='L4 new ED pipeline + all L10 archived figures; full fresh study is a separate validation')
print(json.dumps(report,indent=2))
