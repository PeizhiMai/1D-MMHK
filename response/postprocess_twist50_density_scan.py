#!/usr/bin/env python3
"""Response panels and joint spin/charge velocity plots for the L10/50 scan.

Peak velocities are finite-q/eta-dependent estimators, not q->0 slopes.
Thermodynamic BA references are evaluated at each case's actual mean density.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
from math import ceil

import numpy as np
from matplotlib.lines import Line2D

from twist_averaged_velocity_scan import digest, tag
from mmhk_density_workflow import number_minimizers
from extract_zero_twist_velocities import write_csv
from fast_lorentzian_spectrum import spectrum, peaks, GRID

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'essler_bethe_ansatz'))
from hubbard_velocities import solve

COLORS = dict(charge='#c43c64', spin='#234f91')
STYLES = dict(charge='-', spin='--')


def peak_record(E, W, eta, q):
    found = peaks(E, W, eta, cutoff=5.)
    first, strong = found[0], max(found, key=lambda p: p['height'])
    significant = W > 1e-6*W.sum()
    low = int(np.flatnonzero(significant)[0])
    return dict(first_omega=first['omega'], dominant_omega=strong['omega'],
                first_velocity=first['omega']/q, dominant_velocity=strong['omega']/q,
                lowest_significant_pole=float(E[low]), lowest_pole_residue=float(W[low]),
                lowest_pole_fraction=float(W[low]/W.sum()), integrated_weight=float(W.sum()),
                first_moment=float(E@W), lowest_pole_velocity=float(E[low]/q))


def save_figure(fig, path):
    fig.savefig(path.with_suffix('.png'), dpi=180)
    fig.savefig(path.with_suffix('.svg'))


def response_panels(plt, rows, cache, config, output, eta, full_max):
    low_omega = np.linspace(0, 2, 2001)
    ymax = 1.12*max(spectrum(low_omega, *poles, eta).max() for poles in cache.values())
    for n in config['physical_densities']:
        fig, axes = plt.subplots(4, 4, figsize=(10.8, 7.6), sharex=True, sharey=True)
        fig.subplots_adjust(left=.07, right=.985, bottom=.11, top=.83, hspace=.26, wspace=.15)
        for U, ax in zip(config['U_list'], axes.flat):
            selected = next((r for r in rows if r['U_over_t'] == U and r['target_density'] == n), None)
            if selected is None:
                ax.text(.5, .5, 'Not calculated', transform=ax.transAxes, ha='center', color='#777777')
            else:
                for ch in COLORS:
                    E, W = cache[(U, n, ch)]
                    s = spectrum(low_omega, E, W, eta)
                    ax.plot(low_omega, s, color=COLORS[ch], ls=STYLES[ch], lw=1.5)
                    x = selected[f'{ch}_dominant_omega']
                    if x <= 2:
                        ax.plot(x, spectrum(np.array([x]), E, W, eta)[0], 'o', color=COLORS[ch], ms=3)
                if abs(selected['actual_density']-n) > 1e-10:
                    ax.text(.96, .72, rf'$\bar n={selected["actual_density"]:.3f}$', transform=ax.transAxes, ha='right', fontsize=8)
            ax.text(.96, .94, rf'$U/t={U:g}$', transform=ax.transAxes, ha='right', va='top', fontsize=10)
            ax.set(xlim=(0, 2), ylim=(0, ymax), xticks=[0, 1, 2])
            ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
        for ax in list(axes.flat)[13:]:
            ax.set_visible(False)
        for ax in axes[2, 1:]:
            ax.tick_params(labelbottom=True)
        fig.supylabel(r'$S_{c,s}(q_1,\omega)$', x=.012, fontsize=11)
        fig.supxlabel(r'$\omega/t$', y=.045, fontsize=11)
        fig.suptitle(rf'10-MMHK: target $n={n:g}$, 50 twists, $q_1=\pi/5$, $\eta/t={eta:g}$', fontsize=13, y=.98)
        handles = [Line2D([], [], color=COLORS[ch], ls=STYLES[ch], label=label, lw=1.8)
                   for ch, label in [('charge', 'Charge'), ('spin', r'Spin ($2S^z$)')]]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .94), ncol=2, frameon=False)
        fig.text(.5, .015, 'Absolute intensities; identical axes across all densities and U. Dots: dominant low-energy maxima.', ha='center', fontsize=8)
        save_figure(fig, output/f'response_panel_n{tag(n)}'); plt.close(fig)
        # Supplemental full computed frequency range, never normalized to unit peak.
        for page, group in enumerate([config['U_list'][:5], config['U_list'][5:10], config['U_list'][10:]], 1):
            fig, axes = plt.subplots(len(group), 2, figsize=(9.4, 1.65*len(group)+1.65), squeeze=False)
            fig.subplots_adjust(left=.08, right=.985, bottom=.12, top=.83, hspace=.32, wspace=.17)
            full = np.linspace(0, full_max, round(full_max*100)+1)
            for i, U in enumerate(group):
                row = next((r for r in rows if r['U_over_t'] == U and r['target_density'] == n), None)
                for j, ax in enumerate(axes[i]):
                    if row:
                        for ch in COLORS:
                            E, W = cache[(U, n, ch)]
                            omega = low_omega if j == 0 else full
                            ax.plot(omega, spectrum(omega, E, W, eta), color=COLORS[ch], ls=STYLES[ch], lw=1.3)
                    ax.text(.97, .88, rf'$U/t={U:g}$', transform=ax.transAxes, ha='right', fontsize=9)
                    if row and abs(row['actual_density']-n) > 1e-10:
                        ax.text(.97, .70, rf'$\bar n={row["actual_density"]:.3f}$', transform=ax.transAxes, ha='right', fontsize=8)
                    ax.set(xlim=(0, 2 if j == 0 else full_max), ylabel=r'$S(q_1,\omega)$')
                    if j == 0:
                        ax.set_ylim(0, ymax)
                    else:
                        ax.set_yscale('log'); ax.set_ylim(1e-10, ymax)
                    ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
                    if i == len(group)-1:
                        ax.set_xlabel(r'$\omega/t$')
            fig.suptitle(rf'$n_{{\rm target}}={n:g}$: low-energy and full-range spectra, 50 twists', fontsize=12, y=.98)
            fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .94), ncol=2, frameon=False)
            fig.text(.5, .02, r'Left: absolute intensity, common scale. Right: log intensity, entire stored pole range plus tails. $eta/t=0.2$.', ha='center', fontsize=8)
            save_figure(fig, output/f'full_spectra_n{tag(n)}_page{page}'); plt.close(fig)


def velocity_panel(plt, rows, config, output, eta):
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.15), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.99, bottom=.25, top=.76, wspace=.18)
    ymax = max(2.55, 1.08*max(r[f'{ch}_{quantity}'] for r in rows for ch in COLORS
                            for quantity in ['dominant_velocity', 'first_velocity', 'BA_velocity_actual_n']))
    for ax, n in zip(axes, config['physical_densities']):
        chosen = sorted([r for r in rows if r['target_density'] == n], key=lambda r: r['U_over_t'])
        U = [r['U_over_t'] for r in chosen]
        for ch, marker in [('charge', 'o'), ('spin', 's')]:
            ax.plot(U, [r[f'{ch}_BA_velocity_actual_n'] for r in chosen], color=COLORS[ch], ls=STYLES[ch], lw=1.8)
            ax.plot(U, [r[f'{ch}_dominant_velocity'] for r in chosen], marker, color=COLORS[ch], ms=4)
            different = [r for r in chosen if abs(r[f'{ch}_first_omega']-r[f'{ch}_dominant_omega']) > 1e-5]
            if different:
                ax.plot([r['U_over_t'] for r in different], [r[f'{ch}_first_velocity'] for r in different],
                        '^', color=COLORS[ch], ms=5, mfc='white', mew=1.)
        actual = [r['actual_density'] for r in chosen]
        if any(abs(x-n) > 1e-10 for x in actual):
            ax.text(.04, .06, rf'$\bar n\in[{min(actual):.3f},{max(actual):.3f}]$', transform=ax.transAxes, fontsize=8)
        ax.set(title=rf'$n_{{\rm target}}={n:g}$', xlabel=r'$U/t$', xlim=(-.25, 12.4), ylim=(0, ymax), xticks=[0, 3, 6, 9, 12])
        ax.grid(alpha=.17); ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel(r'$v/(t a_0/\hbar)$')
    fig.suptitle('Spin and charge velocities: 10-MMHK, 50 twists', fontsize=13, y=.98)
    handles = [Line2D([], [], color=COLORS[ch], ls=STYLES[ch], label=f'BA {ch}', lw=1.8) for ch in COLORS]
    handles += [Line2D([], [], color=COLORS[ch], marker=marker, ls='none', label=f'Peak/q, {ch}', ms=5)
                for ch, marker in [('charge', 'o'), ('spin', 's')]]
    if any(abs(r[f'{ch}_first_omega']-r[f'{ch}_dominant_omega']) > 1e-5 for r in rows for ch in COLORS):
        handles += [Line2D([], [], color='#555555', marker='^', mfc='white', ls='none', label='First peak when different', ms=5)]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .92), ncol=3, fontsize=9, frameon=False)
    fig.text(.5, .06, rf'$\widehat v=\omega_{{\rm peak}}(q_1)/q_1$, $q_1=\pi/5$, $\eta/t={eta:g}$: finite-q/eta-dependent estimates, not a $q\to0$ extrapolation.'
             '\nBA lines are thermodynamic slopes evaluated at the actual mean density for each U.', ha='center', fontsize=8.5)
    save_figure(fig, output/'spin_charge_velocity_vs_U'); plt.close(fig)


def density_panel(plt, density_rows, output):
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.65), sharex=True)
    fig.subplots_adjust(left=.07, right=.99, top=.82, bottom=.22, wspace=.28)
    for ax, n in zip(axes, [.7, .8, .9]):
        chosen = [r for r in density_rows if r['target_n'] == n]
        U = [r['U'] for r in chosen]
        actual = [r['actual_n'] for r in chosen]
        ax.plot(U, actual, 'o-', color='#234f91', lw=1.3, ms=4)
        ax.axhline(n, color='#777777', ls='--', lw=1)
        ax.set(title=rf'Target $n={n:g}$', xlabel=r'$U/t$', xlim=(-.25, 12.25),
               ylim=(n-.009, n+.009), xticks=[0, 3, 6, 9, 12], ylabel=r'$\bar n=N_e/L$')
        ax.grid(alpha=.17); ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Density attained by one common chemical potential, 50 twists', fontsize=12, y=.97)
    fig.text(.5, .055, r'No separate per-twist density constraint. U=0 uses the full grand ground projector if the Fermi level is degenerate.', ha='center', fontsize=8.5)
    save_figure(fig, output/'actual_density_vs_U'); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9')
    p.add_argument('--eta', type=float, default=.2)
    p.add_argument('--allow-partial', action='store_true')
    p.add_argument('--verify-raw', action='store_true')
    a = p.parse_args(); config = json.loads((a.input/'config.json').read_text())
    output = a.input/('analysis_partial' if a.allow_partial else 'analysis'); output.mkdir(exist_ok=True)
    rows, cache, sensitivity, checks, populations, hashes, density_rows = [], {}, [], [], [], {}, []
    missing = []; tolerance = json.loads((a.input/'density_policy.json').read_text()).get('nearest_stable_tolerance', 0.) if (a.input/'density_policy.json').exists() else 0.
    for U in config['U_list']:
        energy_file = a.input/f'U{U:g}'/'ground_energy_table.npz'
        if not energy_file.exists():
            E0 = None
        else:
            with np.load(energy_file) as table:
                E0 = table['E0_by_N'].copy()
            hashes[str(energy_file.resolve())] = digest(energy_file)
        for n in config['physical_densities']:
            directory = a.input/f'U{U:g}'/f'n{tag(n)}'
            if not (directory/'COMPLETE').exists():
                missing.append(dict(U=U, density=n)); continue
            m = json.loads((directory/'metadata.json').read_text())
            actual = m['actual_density_Ne_over_L']; q = m['q_min']; mu = m['mu_selection']['mu']
            assert m['L'] == 10 and m['twist_count'] == 50 and abs(q-np.pi/5) < 1e-12
            assert abs(actual-n) <= tolerance+1e-10
            if U > 0:
                counts = [b['N'] for b in m['twists']]
                assert abs(np.mean(counts)/10-actual) < 1e-12
                mins = number_minimizers(E0, mu)
                assert all([int(N)] == matches for N, matches in zip(counts, mins))
                for b in m['twists']:
                    assert abs(sum(c['probability'] for c in b['components'])-1) < 1e-12
                    assert b['ground_dimension'] == sum(c['ground_dimension']*c['spin_flip_multiplicity'] for c in b['components'])
                max_weight = max(c['weight_error'] for c in m['source_moment_checks'])
                max_moment = max(c['first_moment_error'] for c in m['source_moment_checks'])
                residual = max(g['residual'] for b in m['twists'] for c in b['components'] for g in c['ground_states'])
                assert max_weight < 1e-9 and max_moment < 1e-8 and residual < 1e-7
            else:
                max_weight = max_moment = residual = 0.
            row = dict(U_over_t=U, target_density=n, actual_density=actual, L=10, twist_count=50, mu=mu, eta=a.eta, q_min=q)
            ba_actual, ba_nominal = solve(U, actual, 400), solve(U, n, 400)
            data_file = directory/'averaged_qmin_poles.npz'; hashes[str(data_file.resolve())] = digest(data_file)
            with np.load(data_file) as d:
                assert abs(float(d['actual_density'])-actual) < 1e-12
                for ch in COLORS:
                    ee, ww = d[f'{ch}_energies'], d[f'{ch}_weights']
                    assert np.all(ee > 0) and np.all(ww > 0)
                    cache[(U, n, ch)] = ee.copy(), ww.copy()
                    rec = peak_record(ee, ww, a.eta, q)
                    row.update({f'{ch}_{key}': value for key, value in rec.items()})
                    row[f'{ch}_BA_velocity_actual_n'] = ba_actual[f'v_{ch}']
                    row[f'{ch}_BA_velocity_nominal_n'] = ba_nominal[f'v_{ch}']
                    row[f'{ch}_relative_error_actual_BA'] = rec['dominant_velocity']/ba_actual[f'v_{ch}']-1
                    for eta in [.05, .1, .2, .3]:
                        rec2 = peak_record(ee, ww, eta, q)
                        sensitivity.append(dict(U=U, target_density=n, actual_density=actual, channel=ch, eta=eta, **rec2))
                fsum_error = abs(d['charge_energies']@d['charge_weights']-d['spin_energies']@d['spin_weights'])
                assert fsum_error < 1e-8
                if U == 0:
                    np.testing.assert_array_equal(d['charge_energies'], d['spin_energies'])
                    np.testing.assert_array_equal(d['charge_weights'], d['spin_weights'])
                    assert d['charge_energies'].max() <= 4*np.sin(np.pi/10)+1e-10
            checks.append(dict(U=U, target_density=n, actual_density=actual, maximum_weight_error=max_weight,
                               maximum_first_moment_error=max_moment, maximum_ground_residual=residual, fsum_channel_difference=float(fsum_error)))
            rows.append(row)
            density_rows.append(dict(U=U, target_n=n, actual_n=actual, difference=actual-n, mu=mu,
                                     mu_left=m['mu_selection'].get('mu_left'), mu_right=m['mu_selection'].get('mu_right'),
                                     policy=m['mu_selection']['policy']))
            for b in m['twists']:
                populations.append(dict(U=U, target_density=n, actual_density=actual, twist_index=b['twist_index'],
                                        theta=b['theta_boundary'], mean_N=b['N'], ground_dimension=b['ground_dimension'],
                                        N_support=json.dumps(b.get('N_support', {str(b['N']): 1.})),
                                        Sz_nonnegative=';'.join(str((c['Nup']-c['Ndn'])/2) for c in b['components'])))
            if a.verify_raw:
                for path, expected in m.get('raw_fraction_sha256', {}).items():
                    assert digest(path) == expected
            print(f'U={U:2g} target={n:g} actual={actual:.3f}: '
                  f'vc={row["charge_dominant_velocity"]:.6f}/BA={ba_actual["v_charge"]:.6f}, '
                  f'vs={row["spin_dominant_velocity"]:.6f}/BA={ba_actual["v_spin"]:.6f}', flush=True)
    if missing and not a.allow_partial:
        raise RuntimeError(f'Missing response cases: {missing}')
    write_csv(output/'velocity_comparison.csv', rows); write_csv(output/'broadening_sensitivity.csv', sensitivity)
    write_csv(output/'actual_density.csv', density_rows); write_csv(output/'selected_ground_states.csv', populations)
    (output/'validation.json').write_text(json.dumps(dict(passed=not missing, missing=missing, case_checks=checks,
                                                         raw_hashes_verified=a.verify_raw, student_solver_executed=False), indent=2)+'\n')
    (output/'metadata.json').write_text(json.dumps(dict(source_sha256=digest(__file__), input_sha256=hashes,
        BA_source_sha256=digest(HERE/'essler_bethe_ansatz/hubbard_velocities.py'), eta=a.eta,
        estimator='dominant peak of averaged positive Lehmann spectrum/q_min; finite q, not sound slope',
        q_min=np.pi/5, density_policy_tolerance=tolerance, references='thermodynamic BA at actual mean density',
        peak_policy='0<omega<5, exact-refined peaks; height>=0.1%, prominence>=1%; first and lowest significant pole also reported',
        spectral_grid_method='linear residue deposition and FFT Lorentzian convolution; all poles retained; exact derivative-refined peak frequencies',
        spectral_energy_grid_step=GRID, spectral_grid_absolute_error_bound='sum(W)*h^2/(4*pi*eta^3); double if off-grid frequency interpolation needed',
        fast_spectrum_sha256=digest(HERE/'fast_lorentzian_spectrum.py')), indent=2)+'\n')
    os.environ.setdefault('MPLCONFIGDIR', str(HERE/'.mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10})
    full_max = max(60, 5*ceil(max(E.max() for E, W in cache.values())/5)+5)
    response_panels(plt, rows, cache, config, output, a.eta, full_max)
    velocity_panel(plt, rows, config, output, a.eta)
    density_panel(plt, density_rows, output)
    low_omega = np.linspace(0, 2, 2001); full_omega = np.linspace(0, full_max, round(full_max*100)+1)
    for U, n, ch in cache:
        E, W = cache[(U, n, ch)]
        # Preserve full-frequency arrays as well as the compact plotted view.
        np.savez_compressed(output/f'U{U:g}_n{tag(n)}_{ch}_spectra.npz', omega_low=low_omega,
                            S_low=spectrum(low_omega, E, W, a.eta), omega_full=full_omega,
                            S_full=spectrum(full_omega, E, W, a.eta), eta=a.eta)
    print('OUTPUT:', output, flush=True)


if __name__ == '__main__':
    main()
