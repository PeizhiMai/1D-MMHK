#!/usr/bin/env python3
"""Replot the completed 50-twist scan with directly calculated fine-grid BA.

This script NEVER invokes a Hubbard solver. It verifies and reads existing
response spectra/poles and peak estimates, preserving the original analysis.
BA lines have fixed nominal density; actual-density corrections are markers.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'essler_bethe_ansatz'))
from hubbard_velocities import solve
from fast_lorentzian_spectrum import spectrum

COLORS = {'charge': '#c43c64', 'spin': '#234f91'}
STYLES = {'charge': '-', 'spin': '--'}
MARKERS = {'charge': 'o', 'spin': 's'}
NOTE = HERE.parent / 'note/spin_charge_response_definitions.tex'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def tag(n):
    return f'{n:g}'.replace('.', 'p')


def save(fig, output, name):
    fig.savefig(output / f'{name}.png', dpi=180, facecolor='white')
    fig.savefig(output / f'{name}.svg', facecolor='white')


def velocity_plot(plt, rows, ba, densities, output, step, eta=.2):
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.5), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.99, bottom=.25, top=.74, wspace=.18)
    correction_count = 0
    for ax, n in zip(axes, densities):
        chosen = sorted([r for r in rows if r['target_density'] == n], key=lambda r: r['U_over_t'])
        reference = [r for r in ba if r['density'] == n]
        uu = [r['U_over_t'] for r in reference]
        changed = [r for r in chosen if abs(r['actual_density'] - n) > 1e-10]
        correction_count += len(changed)
        for ch in COLORS:
            ax.plot(uu, [r[f'v_{ch}'] for r in reference], color=COLORS[ch],
                    ls=STYLES[ch], lw=1.9, zorder=2)
            ax.plot([r['U_over_t'] for r in chosen],
                    [r[f'{ch}_dominant_velocity'] for r in chosen],
                    marker=MARKERS[ch], ls='none', color=COLORS[ch], ms=4.5, zorder=4)
            # Do not connect changes of n(U) into a purported fixed-n BA line.
            for r in changed:
                ax.plot([r['U_over_t']] * 2,
                        [r[f'{ch}_BA_velocity_nominal_n'], r[f'{ch}_BA_velocity_actual_n']],
                        color=COLORS[ch], lw=.8, alpha=.65, zorder=3)
            if changed:
                ax.plot([r['U_over_t'] for r in changed],
                        [r[f'{ch}_BA_velocity_actual_n'] for r in changed],
                        marker='D', ls='none', color=COLORS[ch], mfc='none',
                        mew=1.2, ms=7, zorder=5)
        if changed:
            ax.text(.045, .065, '$\\bar n=0.896,\\ 0.902,\\ 0.898$\n'
                    '$U/t=7,\\ 9,\\ 10$, respectively',
                    transform=ax.transAxes, fontsize=8)
        ax.set(title=rf'$n_{{\rm target}}={n:g}$', xlabel=r'$U/t$',
               xlim=(-.25, 12.4), ylim=(0, 2.55), xticks=[0, 3, 6, 9, 12])
        ax.grid(alpha=.17)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel(r'$v/(t a_0/\hbar)$')
    fig.suptitle('Spin and charge velocities: 10-MMHK, 50 twists', y=.98, fontsize=13)
    handles = [Line2D([], [], color=COLORS[ch], ls=STYLES[ch], lw=1.9,
                      label=f'BA {ch}: fixed n') for ch in COLORS]
    handles += [Line2D([], [], color=COLORS[ch], marker=MARKERS[ch], ls='none',
                       ms=5, label=f'Peak/q, {ch}') for ch in COLORS]
    if correction_count:
        handles.append(Line2D([], [], color='#555555', marker='D', ls='none',
                              mfc='none', ms=7, label='BA at actual density'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .925),
               ncol=3, fontsize=9, frameon=False)
    fig.text(.5, .07, rf'BA: directly solved at $\Delta(U/t)={step:g}$; lines use fixed nominal density, diamonds use actual density.'
             '\n' + rf'$\widehat v=\omega_{{\rm peak}}(q_1)/q_1$, $q_1=\pi/5$, $\eta/t={eta:g}$: finite-q/eta estimates, not a $q\to0$ extrapolation.',
             ha='center', fontsize=8.5)
    save(fig, output, 'spin_charge_velocity_vs_U')
    plt.close(fig)
    return correction_count


def response_plots(plt, rows, spectra, densities, couplings, output, full=False, eta=.2):
    from matplotlib.lines import Line2D

    mode = 'full' if full else 'low'
    # One global absolute-intensity scale across all U and all densities.
    ymax = 1.12 * max(d[f'S_{mode}'].max() for d in spectra.values())
    if not full:
        ymax = .05 * np.ceil(ymax / .05)
    record = []
    for n in densities:
        fig, axes = plt.subplots(4, 4, figsize=(10.8, 8.0), sharex=True, sharey=True)
        fig.subplots_adjust(left=.073, right=.985, bottom=.115, top=.855, hspace=.29, wspace=.17)
        for U, ax in zip(couplings, axes.flat):
            row = next(r for r in rows if r['target_density'] == n and r['U_over_t'] == U)
            for ch in COLORS:
                data = spectra[(U, n, ch)]
                omega, values = data[f'omega_{mode}'], data[f'S_{mode}']
                ax.plot(omega, values, color=COLORS[ch], ls=STYLES[ch], lw=1.5)
                if not full:
                    peak = row[f'{ch}_dominant_omega']
                    assert omega[0] < peak < omega[-1]
                    ax.plot(peak, np.interp(peak, omega, values), 'o', color=COLORS[ch], ms=3)
            ax.text(.96, .94, rf'$U/t={U:g}$', transform=ax.transAxes,
                    ha='right', va='top', fontsize=10)
            if abs(row['actual_density'] - n) > 1e-10:
                ax.text(.96, .74, rf'$\bar n={row["actual_density"]:.3f}$',
                        transform=ax.transAxes, ha='right', va='top', fontsize=8.5)
            if full:
                ax.set(yscale='log', ylim=(1e-6, ymax), xlim=(0, omega[-1]),
                       xticks=[0, 20, 40, 60], yticks=[1e-6, 1e-4, 1e-2])
            else:
                ax.set(xlim=(0, 2), ylim=(0, ymax), xticks=[0, 1, 2])
            ax.tick_params(labelsize=8.5)
            ax.spines[['top', 'right']].set_visible(False)
        for ax in list(axes.flat)[13:]:
            ax.set_visible(False)
        # U9--U11 have no live axes below them; keep their frequency ticks.
        for ax in axes[2, 1:]:
            ax.tick_params(labelbottom=True)
        handles = [Line2D([], [], color=COLORS[ch], ls=STYLES[ch], lw=1.8, label=label)
                   for ch, label in [('charge', 'Charge'), ('spin', r'Spin ($2S^z$)')]]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .94),
                   ncol=2, fontsize=10, frameon=False)
        fig.suptitle(rf'10-MMHK: target $n={n:g}$, 50 twists, $q_1=\pi/5$, $\eta/t={eta:g}$',
                     fontsize=13, y=.98)
        fig.supylabel(r'$S_{c,s}(q_1,\omega)$', x=.012, fontsize=12)
        fig.supxlabel(r'$\omega/t$', y=.052, fontsize=12)
        description = ('Full stored frequency range plus tails; log intensity; identical axes for all densities and U.'
                       if full else 'Low-energy view; absolute intensities; identical axes for all densities and U. Dots: dominant maxima.')
        fig.text(.5, .014, description, ha='center', fontsize=8.5)
        name = f'response_panel_n{tag(n)}' + ('_full_range' if full else '')
        save(fig, output, name)
        plt.close(fig)
        record.append({'target_density': n, 'panel_count': len(couplings),
                       'U_list': couplings, 'omega_range': [0, float(omega[-1])],
                       'y_scale': 'log' if full else 'linear', 'y_limits': [1e-6 if full else 0, ymax],
                       'channels': list(COLORS), 'file_stem': name})
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=HERE / 'results/twist50_common_mu_U0to12_n0p7_0p8_0p9')
    p.add_argument('--output', type=Path)
    p.add_argument('--ba-step', type=float, default=.2)
    a = p.parse_args()
    if a.ba_step <= 0 or abs(12 / a.ba_step - round(12 / a.ba_step)) > 1e-8:
        p.error('BA step must be positive and divide the U/t=0..12 interval')
    root = a.input.resolve()
    previous = root / 'analysis'
    output = (a.output or root / f'analysis_ba_step{tag(a.ba_step)}').resolve()
    if output == previous:
        p.error('Use a separate output directory to preserve the original analysis')
    output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    note_before = digest(NOTE)
    config = json.loads((root / 'config.json').read_text())
    old_meta = json.loads((previous / 'metadata.json').read_text())
    assert config['L'] == 10 and config['twist_count'] == 50
    assert config['U_list'] == list(range(13)) and config['physical_densities'] == [.7, .8, .9]
    assert abs(old_meta['eta'] - .2) < 1e-12
    assert digest(HERE / 'essler_bethe_ansatz/hubbard_velocities.py') == old_meta['BA_source_sha256']
    with (previous / 'velocity_comparison.csv').open() as f:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]
    assert len(rows) == 39
    inputs = {str(previous / 'velocity_comparison.csv'): digest(previous / 'velocity_comparison.csv'),
              str(previous / 'metadata.json'): digest(previous / 'metadata.json'),
              str(root / 'config.json'): digest(root / 'config.json')}
    spectra, cache_error = {}, 0.
    for row in rows:
        U, n = int(row['U_over_t']), row['target_density']
        assert row['L'] == 10 and row['twist_count'] == 50 and abs(row['q_min'] - np.pi / 5) < 1e-12
        assert row['eta'] == .2
        poles = root / f'U{U}' / f'n{tag(n)}' / 'averaged_qmin_poles.npz'
        inputs[str(poles)] = digest(poles)
        assert inputs[str(poles)] == old_meta['input_sha256'][str(poles)]
        with np.load(poles) as pole_data:
            assert abs(float(pole_data['actual_density']) - row['actual_density']) < 1e-12
            for ch in COLORS:
                f = previous / f'U{U}_n{tag(n)}_{ch}_spectra.npz'
                inputs[str(f)] = digest(f)
                with np.load(f) as d:
                    data = {key: d[key].copy() for key in d.files}
                assert float(data['eta']) == .2
                for mode in ['low', 'full']:
                    xx, yy = data[f'omega_{mode}'], data[f'S_{mode}']
                    assert xx.shape == yy.shape and np.all(np.diff(xx) > 0)
                    assert np.isfinite(yy).all() and np.min(yy) >= 0
                    recomputed = spectrum(xx, pole_data[f'{ch}_energies'], pole_data[f'{ch}_weights'], .2)
                    error = float(np.max(np.abs(yy - recomputed)))
                    cache_error = max(cache_error, error)
                    np.testing.assert_allclose(yy, recomputed, rtol=1e-12, atol=1e-14)
                spectra[(U, n, ch)] = data
    assert len(spectra) == 78
    print('Verified all 78 cached channel spectra against the unchanged original poles.', flush=True)

    couplings = np.round(np.linspace(0, 12, round(12 / a.ba_step) + 1), 10)
    ba, mesh = [], []
    for n in config['physical_densities']:
        for U in couplings:
            result, refined = solve(U, n, 400), solve(U, n, 600)
            errors = {ch: abs(result[f'v_{ch}'] - refined[f'v_{ch}']) for ch in COLORS}
            assert max(errors.values()) < 2e-7
            mesh.append({'density': n, 'U_over_t': float(U), 'velocity_absolute_errors': errors,
                         'nodes': result['nodes'], 'refined_nodes': refined['nodes']})
            ba.append(result)
            if abs(U - round(U)) < 1e-10:
                old_row = next(r for r in rows if r['target_density'] == n and r['U_over_t'] == round(U))
                for ch in COLORS:
                    assert abs(result[f'v_{ch}'] - old_row[f'{ch}_BA_velocity_nominal_n']) < 1e-10
        print(f'n={n:g}: {len(couplings)} direct BA points, U/t=0..12, step={a.ba_step:g}', flush=True)
    write_csv(output / 'BA_velocities_fine_grid.csv', ba)
    shutil.copy2(previous / 'velocity_comparison.csv', output / 'velocity_comparison.csv')

    os.environ.setdefault('MPLCONFIGDIR', str(HERE / '.mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
    correction_count = velocity_plot(plt, rows, ba, config['physical_densities'], output, a.ba_step)
    panels = response_plots(plt, rows, spectra, config['physical_densities'], config['U_list'], output)
    panels += response_plots(plt, rows, spectra, config['physical_densities'], config['U_list'], output, full=True)
    assert digest(NOTE) == note_before
    assert all(digest(path) == sha for path, sha in inputs.items())
    checks = dict(passed=True, response_cases=39, cached_channel_spectra=78,
                  maximum_cached_spectrum_recomputation_error=cache_error,
                  BA_point_count=len(ba), BA_points_per_density=len(couplings),
                  BA_all_point_mesh_checks=mesh,
                  maximum_BA_mesh_velocity_error=max(max(r['velocity_absolute_errors'].values()) for r in mesh),
                  BA_integer_U_matches_original_nominal_reference=True,
                  actual_density_correction_cases=correction_count,
                  response_panels=panels, input_files_unchanged=True,
                  open_latex_note_unchanged=True, Hubbard_solver_executed=False,
                  student_solver_executed=False)
    (output / 'validation.json').write_text(json.dumps(checks, indent=2) + '\n')
    metadata = dict(source=str(Path(__file__).resolve()), source_sha256=digest(__file__),
                    BA_source=str(HERE / 'essler_bethe_ansatz/hubbard_velocities.py'),
                    BA_source_sha256=old_meta['BA_source_sha256'], input_sha256=inputs,
                    L=10, twist_count=50, eta=.2, q_min=np.pi / 5,
                    BA_U_step=a.ba_step, BA_grid='direct integral-equation solutions; not interpolated coarse velocities',
                    BA_line_density='fixed nominal n=Ne/L',
                    BA_diamonds='actual density only where different from target',
                    response_source='unchanged cached twist/ground-averaged positive inelastic Lehmann spectra',
                    spin_operator='2Sz; normalized the same way as charge',
                    velocity_estimator='dominant peak/q_min; finite-q and eta-dependent',
                    BA_book_equations=['5.102', '5.103', '5.104', '5.105', '6.44', '6.45', '6.46', '6.47'],
                    open_latex_note_sha256=note_before, elapsed_seconds=time.monotonic() - start)
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    (output / 'README.txt').write_text(
        'Replot of the existing L=10, 50-twist U/t=0..12 spin/charge response campaign.\n\n'
        f'BA: {len(couplings)} directly calculated U points per density, Delta(U/t)={a.ba_step:g}.\n'
        'The U grid is five times finer than the original step-1 grid when step=0.2.\n'
        'Lines: fixed nominal n=Ne/L=0.7,0.8,0.9. Not splines of the old 13 points.\n'
        'Open diamonds: BA at actual density when different (U7:0.896, U9:0.902, U10:0.898).\n'
        'No density interpolation and no arbitrary ensemble mixing introduced.\n'
        'BA mesh convergence checked at every point, using 400/600 base quadrature nodes.\n'
        'Source: Essler et al., The One-Dimensional Hubbard Model, Eqs.5.102--5.105 and6.44--6.47.\n'
        'Thermodynamic T=B=0; u=U/(4t); t=a0=hbar=1.\n\n'
        'Response: one 13-panel figure for each density; spin and charge overlaid in every panel.\n'
        'U/t=0 through12, eta/t=0.2, q_min=pi/5. Spin operator is longitudinal2Sz.\n'
        'Identical absolute-intensity scales across densities and U; no per-panel normalization.\n'
        'response_panel_n*: low-energy omega/t=0..2 view with dominant-peak dots.\n'
        'response_panel_n*_full_range: full saved omega/t=0..65, log-intensity supplemental view.\n'
        'The log view displays intensities down to1e-6; smaller stored tails remain in source NPZs.\n'
        'Peak/q estimates retain their original finite-q/eta interpretation and values.\n\n'
        'Original analysis, native solver, raw poles, spectra and open LaTeX note were not changed.\n'
        'No Hubbard or student solver was run. All cached spectra were checked against original poles.\n'
        'Files: PNG/SVG figures; BA_velocities_fine_grid.csv; unchanged velocity_comparison.csv;\n'
        'validation.json and metadata.json. Source spectrum NPZs remain in ../analysis.\n')
    print(f'OUTPUT: {output}\nBA mesh max error: {checks["maximum_BA_mesh_velocity_error"]:.3e}\n'
          f'Elapsed: {metadata["elapsed_seconds"]:.1f}s', flush=True)


if __name__ == '__main__':
    main()
