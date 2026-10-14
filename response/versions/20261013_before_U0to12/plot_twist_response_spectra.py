#!/usr/bin/env python3
"""Show all saved n=0.8 q_min responses and audit peak versus spectral edge.

Read-only postprocessing of the twenty-twist common-mu ground traces. No
Hamiltonian solver, student code, or LaTeX editor is used. Intensities are not
independently normalized between charge and longitudinal 2Sz channels.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
from scipy.signal import find_peaks

from extract_zero_twist_velocities import write_csv
from postprocess_twist_velocity_scan import spectrum, peaks
from postprocess_response import read_fraction, poles_from_fraction

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit_weak_group_origins(root, bands, cache):
    """Locate the weak charge group without assigning a spinon branch label."""
    rows, hashes = [], {}
    for band in bands:
        U = band['U_over_t']; cut = band['separating_gap_midpoint']
        m = json.loads((root/f'U{U:g}'/'n0p8/metadata.json').read_text())
        total = cache[(U, 'charge')][1].sum(); sum_low = 0.; singlet_low = 0.
        for block in m['twists']:
            low = 0.
            highest_M = max((c['Nup']-c['Ndn'])/2 for c in block['components'])
            for c in block['components']:
                for path in Path(c['directory']).glob('charge_gsQ*_q1.dat'):
                    hashes[str(path.resolve())] = sha(path)
                    assert hashes[str(path.resolve())] == m['raw_fraction_sha256'][str(path.resolve())]
                    d = read_fraction(path); E, W = poles_from_fraction(d)
                    factor = c['probability']*d['ensemble_weight']/20
                    low += factor*W[(E > 1e-8) & (E < cut)].sum()
            sum_low += low
            if highest_M == 0:
                singlet_low += low
            rows.append({'U_over_t': U, 'theta_over_pi': block['theta_boundary']/np.pi,
                         'highest_ground_spin_projection': highest_M,
                         'low_charge_group_weight': float(low),
                         'low_charge_group_percent_of_averaged_charge_weight': float(100*low/total)})
        np.testing.assert_allclose(100*sum_low/total, band['low_group_percent_of_total_weight'], atol=1e-10, rtol=1e-10)
        assert singlet_low/total < 1e-9
        band['singlet_ground_twists_fraction_of_total_charge_weight_in_low_group'] = float(singlet_low/total)
        band['weak_group_origin'] = 'Triplet-ground twists only within numerical precision; this does not establish a thermodynamic spinon identity.'
    return rows, hashes


def make_plots(couplings, cache, rows, output, eta, full_max):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size': 9})
    low = np.linspace(0, 2, 4001)
    wide = np.linspace(0, full_max, round(full_max/.005)+1)
    fig, axes = plt.subplots(len(couplings), 2, figsize=(9.3, 9.7), squeeze=False)
    fig.subplots_adjust(left=.085, right=.98, bottom=.115, top=.875, hspace=.27, wspace=.24)
    spectra = {'omega_low': low, 'omega_full': wide, 'eta': eta, 'q': np.pi/5}
    for i, U in enumerate(couplings):
        for ch, color, style in [('charge', '#c43c64', '-'), ('spin', '#234f91', '--')]:
            E, W = cache[(U, ch)]; record = next(r for r in rows if r['U_over_t'] == U and r['channel'] == ch)
            a, b = spectrum(low, E, W, eta), spectrum(wide, E, W, eta)
            spectra[f'U{U:g}_{ch}_low'] = a; spectra[f'U{U:g}_{ch}_full'] = b
            axes[i, 0].plot(low, a, color=color, ls=style, lw=1.5)
            axes[i, 1].semilogy(wide, b, color=color, ls=style, lw=1.3)
            peak = record['dominant_peak_omega']
            height = float(spectrum(np.array([peak]), E, W, eta)[0])
            axes[i, 0].scatter([peak], [height], color=color, s=16, zorder=5)
            axes[i, 0].axvline(record['lowest_significant_pole_omega'], color=color, ls=':', lw=.8, alpha=.65)
            # Wide-range evaluation must not hide a larger high-frequency maximum.
            assert b.max() <= height*(1+1e-6)
        for j in [0, 1]:
            ax = axes[i, j]
            ax.text(.96, .84, rf'$U/t={U:g}$', transform=ax.transAxes, ha='right', fontsize=10)
            ax.grid(alpha=.14); ax.spines[['top', 'right']].set_visible(False)
            if i == len(couplings)-1:
                ax.set_xlabel(r'$\omega/t$')
        axes[i, 0].set(xlim=(0, 2), ylim=(0, axes[i, 0].get_ylim()[1]*1.07), ylabel=r'$S_\nu(q_1,\omega)$')
        axes[i, 1].set(xlim=(0, full_max), ylim=(1e-6, 1), ylabel=r'$S_\nu(q_1,\omega)$')
    axes[0, 0].set_title('Low-energy window (linear scale)', fontsize=10)
    axes[0, 1].set_title('Entire computed pole-energy range (log scale)', fontsize=10)
    fig.suptitle(rf'$n=0.8$, $L=10$, 20 twists, $q_1=\pi/5$, $\eta/t={eta:g}$', fontsize=12, y=.975)
    legend = [Line2D([0], [0], color='#c43c64', lw=1.5, label='Charge'),
              Line2D([0], [0], color='#234f91', ls='--', lw=1.5, label=r'Spin ($2S^z$)'),
              Line2D([0], [0], marker='o', ls='none', color='#333333', label='First local maximum = dominant peak'),
              Line2D([0], [0], color='#333333', ls=':', label='Lowest significant pole')]
    fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(.53, .947), ncol=2, frameon=False, fontsize=8)
    fig.text(.5, .025, 'Absolute intensities; no independent channel rescaling. Dotted poles require residue > '+r'$10^{-6}$'
             +' of total weight.\nPositive-frequency structure factors; weak high-energy tails shown at the saved Lanczos depth 300.',
             ha='center', fontsize=8)
    fig.savefig(output/'all_U_spin_charge_spectra.png', dpi=180)
    fig.savefig(output/'all_U_spin_charge_spectra.svg')
    plt.close(fig)
    np.savez_compressed(output/'all_U_spin_charge_spectra.npz', **spectra)
    # Discrete weak poles should not be erased just because eta=0.2 hides them.
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1))
    fig.subplots_adjust(left=.085, right=.98, bottom=.24, top=.82, wspace=.25)
    x = np.linspace(.2, 1.6, 5601)
    for ax, ch in zip(axes, ['charge', 'spin']):
        E, W = cache[(8., ch)]
        for eta_i, color, style in [(eta, '#222222', '-'), (.05, '#234f91', '--'), (.005, '#999999', '-')]:
            s = spectrum(x, E, W, eta_i)
            ax.plot(x, s, color=color, ls=style, lw=1 if eta_i == .005 else 1.6,
                    label=rf'$\eta/t={eta_i:g}$')
        r = next(r for r in rows if r['U_over_t'] == 8 and r['channel'] == ch)
        ax.axvline(r['lowest_significant_pole_omega'], color='#c47616', ls=':', lw=1, label='Lowest significant pole')
        ax.set(xlabel=r'$\omega/t$', ylabel=r'$S_\nu(q_1,\omega)$', xlim=(.2, 1.6), title=ch.capitalize())
        ax.set_yscale('log'); ax.set_ylim(1e-4, 2)
        ax.grid(alpha=.13); ax.spines[['top', 'right']].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc='upper center', bbox_to_anchor=(.53, .925), ncol=4, frameon=False)
    fig.suptitle(r'$U/t=8$: lowest poles versus broadened response maxima', fontsize=12, y=.965)
    fig.text(.5, .04, r'Smaller $\eta$ exposes discrete twist poles; it is a resolution check, not a thermodynamic limit.'
             '\nA pole onset, the first local maximum and the maximum of an averaged band are different observables.', ha='center', fontsize=8)
    fig.savefig(output/'U8_peak_vs_edge_broadening.png', dpi=180)
    fig.savefig(output/'U8_peak_vs_edge_broadening.svg')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=HERE/'results/twist20_common_mu_velocity_scan')
    p.add_argument('--eta', type=float, default=.2)
    args = p.parse_args()
    if args.eta <= 0:
        p.error('Require eta>0')
    os.environ.setdefault('MPLCONFIGDIR', str(HERE/'.mplconfig'))
    output = args.input/'analysis/n0p8_spectra_and_peak_audit'; output.mkdir(exist_ok=True)
    config = json.loads((args.input/'config.json').read_text()); couplings = config['U_list']
    cache, rows, bands, threshold_audit, hashes, resolutions, matching = {}, [], [], [], {}, [], []
    full_max = 0
    for U in couplings:
        directory = args.input/f'U{U:g}'/'n0p8'
        assert (directory/'COMPLETE').exists()
        metadata = json.loads((directory/'metadata.json').read_text())
        assert metadata['actual_density_Ne_over_L'] == .8 and metadata['twist_count'] == 20
        path = directory/'averaged_qmin_poles.npz'; hashes[str(path.resolve())] = sha(path)
        with np.load(path) as data:
            for ch in ['charge', 'spin']:
                E, W = data[f'{ch}_energies'].copy(), data[f'{ch}_weights'].copy()
                assert np.all(E > 0) and np.all(W > 0)
                cache[(U, ch)] = E, W; full_max = max(full_max, E.max())
                total = W.sum()
                found = peaks(E, W, args.eta); first, dominant = found[0], max(found, key=lambda p: p['height'])
                x = np.linspace(0, 5, 10001); s = spectrum(x, E, W, args.eta)
                ids, _ = find_peaks(s)  # No height/prominence threshold for this check.
                assert abs(x[ids[0]]-first['omega']) < .0005
                assert abs(first['omega']-dominant['omega']) < 1e-10
                good = W > 1e-6*total; j = np.flatnonzero(good)[0]
                rows.append({'U_over_t': U, 'channel': ch, 'eta': args.eta,
                             'first_local_maximum_omega': first['omega'], 'dominant_peak_omega': dominant['omega'],
                             'lowest_significant_pole_omega': float(E[j]),
                             'lowest_pole_percent_of_total_weight': float(100*W[j]/total),
                             'total_inelastic_Lehmann_weight': float(total),
                             'first_equals_dominant': True, 'unfiltered_first_grid_omega': float(x[ids[0]])})
                for threshold in [1e-8, 1e-6, 1e-5, 1e-4, 1e-3]:
                    keep = np.flatnonzero(W > threshold*total); j = keep[0]
                    threshold_audit.append({'U_over_t': U, 'channel': ch, 'residue_fraction_threshold': threshold,
                                            'lowest_omega': float(E[j]), 'lowest_pole_fraction': float(W[j]/total)})
                for eta in [.005, .02, .05, .1, .2, .3]:
                    candidates = peaks(E, W, eta)
                    strongest = max(candidates, key=lambda p: p['height'])
                    resolutions.append({'U_over_t': U, 'channel': ch, 'eta': eta,
                                        'first_significant_peak': candidates[0]['omega'],
                                        'dominant_significant_peak': strongest['omega']})
                if ch == 'charge':
                    ids = np.flatnonzero((W > 1e-6*total) & (E < dominant['omega']))
                    gap_index = np.argmax(np.diff(E[ids])); lower, upper = E[ids[gap_index]], E[ids[gap_index+1]]
                    cut = .5*(lower+upper)
                    spin_E, spin_W = data['spin_energies'], data['spin_weights']
                    spin_matches = set()
                    weak_ids = np.flatnonzero((E < cut) & (W > 1e-6*total))
                    for j in weak_ids:
                        nearest = int(np.argmin(abs(spin_E-E[j])))
                        error = float(abs(spin_E[nearest]-E[j]))
                        matched = error < 1e-6
                        if matched:
                            spin_matches.add(nearest)
                        matching.append({'U_over_t': U, 'charge_pole_omega': float(E[j]),
                                         'charge_weight_percent': float(100*W[j]/total),
                                         'nearest_spin_pole_omega': float(spin_E[nearest]),
                                         'frequency_error': error, 'matched_within_1e_minus_6': matched,
                                         'spin_weight_percent': float(100*spin_W[nearest]/spin_W.sum())})
                    bands.append({'U_over_t': U, 'low_group_last_omega': float(lower),
                                  'main_band_first_omega': float(upper), 'separating_gap_midpoint': float(cut),
                                  'low_group_percent_of_total_weight': float(100*W[E < cut].sum()/total),
                                  'main_band_to_omega2_percent_of_total_weight': float(100*W[(E >= cut) & (E < 2)].sum()/total),
                                  'matched_spin_poles_percent_of_spin_weight': float(100*spin_W[list(spin_matches)].sum()/spin_W.sum()),
                                  'all_significant_low_charge_poles_match_spin_energies': len(spin_matches)==len(weak_ids),
                                  'definition': 'largest gap between significant poles below dominant peak; does not establish an excitation-branch identity'})
    full_max = float(np.ceil(full_max/5)*5)
    origins, raw_hashes = audit_weak_group_origins(args.input, bands, cache)
    make_plots(couplings, cache, rows, output, args.eta, full_max)
    write_csv(output/'peak_vs_lowest_pole.csv', rows)
    write_csv(output/'charge_low_group_weights.csv', bands)
    write_csv(output/'lowest_pole_threshold_audit.csv', threshold_audit)
    write_csv(output/'peak_broadening_audit.csv', resolutions)
    write_csv(output/'weak_charge_poles_spin_comparison.csv', matching)
    write_csv(output/'weak_charge_group_twist_origins.csv', origins)
    metadata = {'input_pole_sha256': hashes, 'source_sha256': sha(__file__),
                'response_evaluator_sha256': sha(HERE/'postprocess_twist_velocity_scan.py'),
                'eta': args.eta, 'q': np.pi/5, 'positive_frequency_range': [0, full_max],
                'intensity_normalization': 'absolute; physical sources have 1/sqrt(L), spin is longitudinal 2Sz',
                'first_unfiltered_local_maximum_equals_dominant_for_all_channels_and_U': True,
                'weak_group_raw_fraction_sha256': raw_hashes,
                'lowest_significant_pole_default_residue_fraction': 1e-6,
                'no_new_solver_runs': True,
                'interpretation': 'Largest response feature is a mode/ridge proxy, not an automatic sound-velocity definition; lowest nonzero pole is a separate spectral-edge observable.'}
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    shutil.copy2(__file__, output/'plotter_source_snapshot.py')
    print('Peak/edge records:', json.dumps(rows, indent=2), flush=True)
    print('Charge group weights:', json.dumps(bands, indent=2), flush=True)
    print('Plots and audit:', output, flush=True)


if __name__ == '__main__':
    main()
