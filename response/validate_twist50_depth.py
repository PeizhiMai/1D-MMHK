#!/usr/bin/env python3
"""Compare selected canonical response anchors at Lanczos depths 300/600."""
import json
from pathlib import Path
import time

import numpy as np

from twist_averaged_velocity_scan import fractions, digest
from postprocess_twist_velocity_scan import spectrum
from fast_lorentzian_spectrum import peaks

ROOT = Path(__file__).resolve().parent/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'


def main():
    cases = json.loads((ROOT/'validation/depth600/cases.json').read_text())
    rows, hashes = [], {}
    for c in cases:
        if c.get('selected_density_anchor') is False:
            continue
        shallow, deep = Path(c['reference_directory']), Path(c['directory'])
        while not (shallow/'COMPLETE').exists():
            time.sleep(5)
        a, b = [json.loads((p/'metadata.json').read_text()) for p in [shallow, deep]]
        assert (a['L'], a['U'], a['Nup'], a['Ndn']) == (b['L'], b['U'], b['Nup'], b['Ndn'])
        assert abs(a['phase_per_bond']-b['phase_per_bond']) < 1e-14 and abs(a['E0']-b['E0']) < 1e-8
        assert a['ground_subspace_complete'] and b['ground_subspace_complete']
        for ch in ['charge', 'spin']:
            EA, WA = fractions(dict(directory=str(shallow), metadata=a), ch, hashes, [])
            EB, WB = fractions(dict(directory=str(deep), metadata=b), ch, hashes, [])
            for eta in [.05, .2]:
                omega = np.linspace(0, 5, 5001)
                sa, sb = spectrum(omega, EA, WA, eta), spectrum(omega, EB, WB, eta)
                error = float(np.max(abs(sa-sb))); relative = error/max(float(sb.max()), 1e-30)
                found = []
                for energies, weights in [(EA, WA), (EB, WB)]:
                    try:
                        found.append(peaks(energies, weights, eta))
                    except RuntimeError as exception:
                        if str(exception) != 'No resolved positive-frequency local maximum':
                            raise
                        found.append([])
                pa, pb = found
                resolution_agreement = bool(pa) == bool(pb)
                peak_error = (abs(max(pa, key=lambda p: p['height'])['omega']-
                                  max(pb, key=lambda p: p['height'])['omega']) if pa and pb else None)
                row = dict(U=c['U'], N=c['N'], twist_index=c['twist_index'], channel=ch, eta=eta,
                           maximum_absolute_error=error, relative_to_peak_error=relative,
                           dominant_peak_frequency_difference=peak_error,
                           resolved_peak_agreement=resolution_agreement,
                           source_weight_difference=float(abs(WA.sum()-WB.sum())),
                           first_moment_difference=float(abs(EA@WA-EB@WB)))
                rows.append(row); print(row, flush=True)
    primary = all(r['relative_to_peak_error'] < 1e-6 and r['resolved_peak_agreement'] and
                  (r['dominant_peak_frequency_difference'] is None or r['dominant_peak_frequency_difference'] < 1e-6)
                  for r in rows if r['eta'] == .2)
    strict_small_eta = all(r['relative_to_peak_error'] < 1e-6 for r in rows if r['eta'] == .05)
    report = dict(passed_primary_eta0p2=primary, strict_1ppm_small_eta_test_passed=strict_small_eta,
                  window=[0, 5], comparisons=rows, source_fraction_sha256=hashes,
                  validator_sha256=digest(__file__), not_a_global_depth_convergence_proof=True)
    (ROOT/'validation/depth_convergence.json').write_text(json.dumps(report, indent=2)+'\n')
    assert primary, 'Primary eta=.2 spectrum or dominant-peak convergence needs refinement'
    print('PRIMARY ETA=.2 DEPTH VALIDATION PASSED; small-eta discrepancies retained in report', flush=True)


if __name__ == '__main__':
    main()
