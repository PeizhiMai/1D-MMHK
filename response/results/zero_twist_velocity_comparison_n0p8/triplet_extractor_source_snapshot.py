#!/usr/bin/env python3
"""Ne=8 L=10 peak estimates, including the full zero-field ground triplet.

Uses saved physical-density fractions only. The isotropic triplet average of
longitudinal spectra is (S_M0+2*S_M1)/3, since M=-1 is a spin-flipped M=+1.
It is NOT a reduced-BZ sum, a thermal spectrum, or a q->0 velocity fit.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np

from extract_zero_twist_velocities import load_channel, merge_poles, refined_peaks, write_csv
from postprocess_response import lorentzian

HERE = Path(__file__).resolve().parent
BASE0 = HERE/'results/zero_twist_benchmark_n0p8'
BASE1 = HERE/'results/zero_twist_benchmark_n0p8_Sz1'
OUT = HERE/'results/zero_twist_velocity_comparison_n0p8'
BA_FILE = HERE/'essler_bethe_ansatz/results_n0p8/velocities.csv'
NOTE = HERE.parent/'note/spin_charge_response_definitions.tex'
COUPLINGS = [4., 6., 7., 8., 12.]
Q1 = 2*np.pi/10


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def selected(E, W, eta):
    found = refined_peaks(E, W, eta)
    if not found:
        raise ValueError('No resolved local maximum')
    return found[0], max(found, key=lambda peak: peak['height'])


def plot(rows, ba, cache):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(9.1, 4.3), sharey=True)
    fig.subplots_adjust(left=.075, right=.985, bottom=.30, top=.78, wspace=.12)
    U = [r['U_over_t'] for r in rows]
    for ax, channel, color in zip(axes, ['charge', 'spin'], ['#c43c64', '#234f91']):
        valid = [r for r in ba if r['U_over_t'] <= 12]
        ax.plot([r['U_over_t'] for r in valid], [r[f'v_{channel}'] for r in valid],
                color=color, lw=2, label='Thermodynamic Bethe ansatz')
        ax.plot(U, [r[f'{channel}_triplet_dominant_velocity'] for r in rows], 'o',
                color='#222222', ms=5, label=r'$L=10$: equal-triplet average')
        if channel == 'spin':
            ax.plot(U, [r['spin_M0_dominant_velocity'] for r in rows], 's', mfc='none',
                    mec='#c47616', ms=6, label=r'$L=10$: fixed $S^z=0$')
        ax.set(title=channel.capitalize(), xlabel=r'$U/t$', xlim=(0, 12.5), ylim=(0, 2.3))
        ax.grid(alpha=.18); ax.spines[['top', 'right']].set_visible(False)
        ax.legend(fontsize=8, loc='upper right', frameon=False)
    axes[0].set_ylabel(r'$\widehat v/(ta_0/\hbar)$')
    fig.suptitle(r'Zero-twist Hubbard ring: $L=10$, $N_e=8$, $n=0.8$', fontsize=13, y=.97)
    fig.text(.5, .065, r'Dominant-peak estimator: $\widehat v=\omega_{\rm peak}(q_1)/q_1$, $q_1=\pi/5$, $\eta/t=0.2$'
             '\nFinite-q estimates, not a '+r'$q\to0$'+' fit. Ground triplet: '+r'$(S_{M=0}+2S_{M=1})/3$.',
             ha='center', fontsize=9)
    fig.savefig(OUT/'zero_twist_L10_n0p8_ensemble_velocities_vs_bethe_ansatz.png', dpi=180)
    fig.savefig(OUT/'zero_twist_L10_n0p8_ensemble_velocities_vs_bethe_ansatz.svg')
    plt.close(fig)
    # U=8 diagnostic shows the choice of ensemble and weaker charge poles.
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), layout='constrained')
    omega = np.linspace(0, 1.7, 3401)
    for ax, channel, color in zip(axes, ['charge', 'spin'], ['#c43c64', '#234f91']):
        E, W, E0, W0, E1, W1 = cache[(8., channel)]
        values = [lorentzian(omega, E, W, .2), lorentzian(omega, E0, W0, .2),
                  lorentzian(omega, E1, W1, .2)]
        norm = max(s.max() for s in values)
        for s, style, label in zip(values, ['-', '--', ':'],
                                   ['Equal-triplet average', r'Fixed $M=0$', r'Fixed $M=1$']):
            ax.plot(omega, s/norm, color=color, ls=style, lw=1.6, label=label)
        r = next(row for row in rows if row['U_over_t'] == 8)
        ax.axvline(Q1*r[f'{channel}_BA_velocity'], color='#777777', ls='-.', lw=1,
                   label=r'$q_1v_{\rm BA}$ (linear guide)')
        good = (E < 1.7) & (W > 1e-4*W.sum())
        # Pole heights relative to the largest residue are not spectral heights.
        ax.vlines(E[good], 0, .22*W[good]/W[good].max(), color='#222222', lw=1,
                  label='Poles: relative residue, scaled')
        ax.set(title=channel.capitalize(), xlabel=r'$\omega/t$', xlim=(0, 1.7), ylim=(0, 1.12))
        ax.grid(alpha=.15); ax.spines[['top', 'right']].set_visible(False)
        ax.legend(fontsize=7, frameon=False, loc='upper right')
    axes[0].set_ylabel('Intensity / common panel maximum')
    fig.suptitle(r'$U/t=8$, $q_1=\pi/5$, $n=0.8$: ground-ensemble diagnostic', fontsize=12)
    fig.savefig(OUT/'U8_ground_ensemble_and_poles.png', dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(HERE/'.mplconfig'))
    ba = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(BA_FILE.open())]
    assert all(abs(r['density']-.8) < 1e-12 for r in ba)
    BA = {r['U_over_t']: r for r in ba}
    raw_hashes, cache, rows, broadening, pole_rows, validation = {}, {}, [], [], [], []
    manifests = [json.loads((base/'manifest.json').read_text()) for base in [BASE0, BASE1]]
    assert manifests[0]['source_sha256'] == manifests[1]['source_sha256']
    assert manifests[0]['executable_sha256'] == manifests[1]['executable_sha256']
    grid = np.linspace(0, 5, 5001)
    for U in COUPLINGS:
        d0, d1 = [base/f'L10_U{U:g}_n0p4_phi0' for base in [BASE0, BASE1]]
        control = BASE0/f'controls/U{U:g}_Nup6_Ndn2_ground'
        assert (control/'COMPLETE').exists()
        m2 = json.loads((control/'metadata.json').read_text())
        assert (m2['L'], m2['Nup'], m2['Ndn'], m2['U'], m2['phase_per_bond']) == (10, 6, 2, U, 0)
        row = {'U_over_t': U, 'density': .8, 'q_min': Q1, 'eta': .2}
        for channel in ['charge', 'spin']:
            E0, W0, _, m0 = load_channel(d0, channel, 1, raw_hashes, .8)
            E1, W1, _, m1 = load_channel(d1, channel, 1, raw_hashes, .8)
            assert (m0['Nup'], m0['Ndn'], m1['Nup'], m1['Ndn']) == (4, 4, 5, 3)
            assert m0['ground_scan_all_Q'] and m1['ground_scan_all_Q'] and m2['ground_scan_all_Q']
            assert len(m0['ground_states']) == len(m1['ground_states']) == 1
            assert abs(m0['E0']-m1['E0']) < 1e-9
            assert m2['E0']-m0['E0'] > 1e-6  # Rules out S>=2 in the ground state.
            E, W = merge_poles(np.r_[E0, E1], np.r_[W0/3, 2*W1/3])
            cache[(U, channel)] = E, W, E0, W0, E1, W1
            raw0 = lorentzian(grid, E0, W0, .2)
            raw1 = lorentzian(grid, E1, W1, .2)
            scalar_error = float(np.max(abs(raw0-raw1)))
            if channel == 'charge':
                assert scalar_error < 1e-6  # Spin-scalar density is multiplet independent.
            direct = (raw0+2*raw1)/3
            merged = lorentzian(grid, E, W, .2)
            assert np.max(abs(direct-merged)) < 1e-8
            first, dominant = selected(E, W, .2)
            _, dom0 = selected(E0, W0, .2)
            row.update({f'{channel}_triplet_first_omega': first['omega'],
                        f'{channel}_triplet_dominant_omega': dominant['omega'],
                        f'{channel}_triplet_first_velocity': first['omega']/Q1,
                        f'{channel}_triplet_dominant_velocity': dominant['omega']/Q1,
                        f'{channel}_M0_dominant_velocity': dom0['omega']/Q1,
                        f'{channel}_BA_velocity': BA[U][f'v_{channel}'],
                        f'{channel}_triplet_relative_error': dominant['omega']/Q1/BA[U][f'v_{channel}']-1})
            for eta in [.02, .05, .1, .2, .3]:
                first, dominant = selected(E, W, eta)
                broadening.append({'U_over_t': U, 'channel': channel, 'eta': eta,
                                    'first_omega': first['omega'], 'dominant_omega': dominant['omega'],
                                    'first_velocity': first['omega']/Q1,
                                    'dominant_velocity': dominant['omega']/Q1})
            for energy, weight in zip(E, W):
                if weight > 1e-8*W.sum():
                    pole_rows.append({'U_over_t': U, 'channel': channel, 'omega': float(energy),
                                      'weight': float(weight), 'fraction_of_inelastic_weight': float(weight/W.sum())})
            validation.append({'U': U, 'channel': channel, 'M0_M1_ground_energy_error': abs(m0['E0']-m1['E0']),
                               'M2_ground_energy_above_M0': m2['E0']-m0['E0'],
                               'M0_M1_spectrum_max_difference': scalar_error,
                               'direct_vs_merged_average_max_error': float(np.max(abs(direct-merged)))})
        rows.append(row)
        print(f"U={U:g}: triplet charge={row['charge_triplet_dominant_velocity']:.6f}, BA={row['charge_BA_velocity']:.6f}; "
              f"triplet spin={row['spin_triplet_dominant_velocity']:.6f}, M0={row['spin_M0_dominant_velocity']:.6f}, BA={row['spin_BA_velocity']:.6f}", flush=True)
    write_csv(OUT/'triplet_velocity_comparison.csv', rows)
    write_csv(OUT/'triplet_broadening_sensitivity.csv', broadening)
    write_csv(OUT/'triplet_merged_poles.csv', pole_rows)
    plot(rows, ba, cache)
    metadata = {'L': 10, 'Ne': 8, 'physical_density': .8, 'native_filling': .4,
                'phase_per_bond': 0, 'q_min': Q1, 'ground_total_spin': 1,
                'ensemble': '(M=0 longitudinal spectrum + 2*M=1 longitudinal spectrum)/3',
                'negative_M_component': 'Obtained by exact global spin-flip symmetry, not a new solver run',
                'estimator': 'dominant positive-frequency local maximum / q_min',
                'eta': .2, 'peak_window': [0, 5], 'source_sha256': sha(__file__),
                'helper_sha256': sha(HERE/'extract_zero_twist_velocities.py'),
                'BA_csv_sha256': sha(BA_FILE), 'raw_fraction_sha256': raw_hashes,
                'note_sha256': sha(NOTE), 'ground_spin_validation': validation,
                'caveat': 'Finite-q peaks are not a controlled q->0 velocity; weaker charge poles are retained separately.'}
    (OUT/'triplet_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    shutil.copy2(__file__, OUT/'triplet_extractor_source_snapshot.py')
    print('Saved triplet comparison:', OUT, flush=True)


if __name__ == '__main__':
    main()
