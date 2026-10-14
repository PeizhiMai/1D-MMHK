#!/usr/bin/env python3
"""Extract peaks AFTER common-mu twist/ground-trace averaging.

The observable is the positive-frequency density structure factor, as in the
previous zero-twist comparison. These finite-q estimates are not q->0 fits.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
from scipy.optimize import brentq
from scipy.signal import find_peaks

from extract_zero_twist_velocities import write_csv

HERE = Path(__file__).resolve().parent


def spectrum(omega, E, W, eta):
    # Chunking avoids allocating a huge (all twist poles) x (frequency) array.
    values = np.zeros(len(omega))
    for start in range(0, len(E), 256):
        e, w = E[start:start+256], W[start:start+256]
        values += np.sum(w[:, None]*eta/np.pi/((omega[None, :]-e[:, None])**2+eta**2), axis=0)
    return values


def peaks(E, W, eta, cutoff=5.):
    omega = np.linspace(0, cutoff, round(cutoff/.0005)+1)
    s = spectrum(omega, E, W, eta)
    indices, _ = find_peaks(s, height=.001*s.max(), prominence=.01*s.max())
    result = []
    def slope(x):
        z = x-E
        return float(np.sum(-2*W*eta*z/(np.pi*(z*z+eta*eta)**2)))
    for i in indices:
        lo, hi = omega[i-1], omega[i+1]
        if not slope(lo) > 0 > slope(hi):
            raise RuntimeError('A peak could not be bracketed at the advertised grid resolution')
        x = brentq(slope, lo, hi, xtol=1e-12)
        height = float(np.sum(W*eta/np.pi/((x-E)**2+eta**2)))
        result.append({'omega': float(x), 'height': height})
    if not result:
        raise RuntimeError('No resolved positive-frequency local maximum')
    return result


def plots(rows, BA, cache, output, eta, twist_count):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    densities = sorted({r['target_density'] for r in rows})
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(len(densities), 2, figsize=(9.2, 3.0*len(densities)+1.8),
                             squeeze=False, sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.985, top=.78 if len(densities)==1 else .88,
                        bottom=.30 if len(densities)==1 else .22, hspace=.26, wspace=.12)
    for i, n in enumerate(densities):
        selected = [r for r in rows if r['target_density'] == n]
        for j, ch, color in zip([0, 1], ['charge', 'spin'], ['#c43c64', '#234f91']):
            ax = axes[i, j]; ba = [r for r in BA[n] if r['U_over_t'] <= max(r['U_over_t'] for r in rows)]
            ax.plot([r['U_over_t'] for r in ba], [r[f'v_{ch}'] for r in ba], color=color,
                    lw=2, label='Thermodynamic Bethe ansatz')
            ax.plot([r['U_over_t'] for r in selected], [r[f'{ch}_dominant_velocity'] for r in selected],
                    'o', color='#222222', ms=5, label='Twist-averaged dominant peak')
            if any(abs(r[f'{ch}_first_omega']-r[f'{ch}_dominant_omega']) > 1e-5 for r in selected):
                ax.plot([r['U_over_t'] for r in selected], [r[f'{ch}_first_velocity'] for r in selected],
                        'D', mfc='none', mec='#c47616', ms=5, label='First resolved peak')
            ax.set(title=rf'{ch.capitalize()}, $n={n:g}$', xlim=(0, 12.5), ylim=(0, 2.35))
            ax.grid(alpha=.18); ax.spines[['top', 'right']].set_visible(False)
            ax.legend(fontsize=7.5, frameon=False, loc='upper right')
            if j == 0:
                ax.set_ylabel(r'$\widehat v/(ta_0/\hbar)$')
            if i == len(densities)-1:
                ax.set_xlabel(r'$U/t$')
    fig.suptitle(rf'$L=10$ Hubbard ring: {twist_count} twists, common $\mu$, actual ground-state trace', fontsize=12, y=.97)
    fig.text(.5, .045, rf'$\widehat v=\omega_{{\rm peak}}(q_1)/q_1$, $q_1=\pi/5$, $\eta/t={eta:g}$.'
             '\nSpectra averaged before peak extraction; one common '+r'$\mu$'+' and each twist\'s actual number sector and full spin ground trace.'
             '\nFinite-q estimates, not a '+r'$q\to0$'+' fit or finite-size extrapolation.', ha='center', fontsize=8.5)
    fig.savefig(output/'twist20_velocities_vs_bethe_ansatz.png', dpi=180)
    fig.savefig(output/'twist20_velocities_vs_bethe_ansatz.svg')
    plt.close(fig)
    couplings = sorted({r['U_over_t'] for r in rows})
    fig, axes = plt.subplots(len(couplings), len(densities), figsize=(9.1, 2.0*len(couplings)),
                             squeeze=False, layout='constrained')
    omega = np.linspace(0, 2, 4001)
    for i, U in enumerate(couplings):
        for j, n in enumerate(densities):
            ax = axes[i, j]; row = next(r for r in rows if r['U_over_t'] == U and r['target_density'] == n)
            for ch, color, style in [('charge', '#c43c64', '-'), ('spin', '#234f91', '--')]:
                E, W = cache[(U, n, ch)]; s = spectrum(omega, E, W, eta)
                ax.plot(omega, s/s.max(), color=color, ls=style, lw=1.4, label=ch.capitalize())
                ax.axvline(row[f'{ch}_dominant_omega'], color=color, ls=':', lw=.8)
            ax.set(title=rf'$U/t={U:g}$, $n={n:g}$', xlim=(0, 2), ylim=(0, 1.08))
            if i == 0:
                ax.legend(fontsize=8, frameon=False)
            if j == 0:
                ax.set_ylabel('Normalized intensity')
            if i == len(couplings)-1:
                ax.set_xlabel(r'$\omega/t$')
            ax.grid(alpha=.15)
    fig.suptitle(rf'{twist_count}-twist physical density spectra at $q_1=\pi/5$, $\eta/t={eta:g}$'
                 '\nEach channel separately normalized; dotted lines mark dominant peaks.', fontsize=11)
    fig.savefig(output/'twist20_qmin_peak_audit.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(len(couplings), len(densities), figsize=(9.1, 1.7*len(couplings)),
                             squeeze=False, layout='constrained')
    for i, U in enumerate(couplings):
        for j, n in enumerate(densities):
            row = next(r for r in rows if r['U_over_t'] == U and r['target_density'] == n)
            metadata = json.loads((Path(row['directory'])/'metadata.json').read_text())
            theta = [b['theta_boundary'] for b in metadata['twists']]
            N = [b['N'] for b in metadata['twists']]
            order = np.argsort(theta)
            axes[i, j].plot(np.asarray(theta)[order], np.asarray(N)[order], 'o-', color='#234f91', ms=3, lw=1)
            axes[i, j].axhline(10*n, color='#888888', ls=':', lw=1)
            axes[i, j].set(title=rf'$U/t={U:g}$, $n={n:g}$, $\mu/t={row["mu"]:.4f}$',
                           xticks=[-np.pi, 0, np.pi], xticklabels=[r'$-\pi$', '0', r'$\pi$'], ylim=(5.5, 10.5))
            if j == 0:
                axes[i, j].set_ylabel(r'$\langle N\rangle_\theta$')
            if i == len(couplings)-1:
                axes[i, j].set_xlabel(r'Boundary twist $\theta$')
    fig.suptitle('Number sectors selected by one common chemical potential per density and U', fontsize=11)
    fig.savefig(output/'twist20_selected_particle_numbers.png', dpi=160)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=HERE/'results/twist20_common_mu_velocity_scan')
    p.add_argument('--eta', type=float, default=.2)
    p.add_argument('--densities', nargs='+', type=float, help='Process a requested subset without touching other densities')
    args = p.parse_args()
    if args.eta <= 0:
        p.error('Require eta>0')
    output = args.input/'analysis'; output.mkdir(exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(HERE/'.mplconfig'))
    config = json.loads((args.input/'config.json').read_text())
    active_request = args.input/'active_request.json'
    if args.densities is None and active_request.exists():
        args.densities = json.loads(active_request.read_text())['physical_densities']
    if args.densities is not None:
        if not set(args.densities).issubset(config['physical_densities']):
            p.error('Requested densities are not part of this scan')
        config['physical_densities'] = args.densities
    rows, sensitivity, cache, hashes, population_rows, BA = [], [], {}, {}, [], {}
    for n in config['physical_densities']:
        BA[n] = [{key: float(value) for key, value in row.items()} for row in
                 csv.DictReader((HERE/f'essler_bethe_ansatz/results_n{str(n).replace(".", "p")}/velocities.csv').open())]
        assert all(abs(b['density']-n) < 1e-12 for b in BA[n])
    for U in config['U_list']:
        for n in config['physical_densities']:
            directory = args.input/f'U{U:g}'/f'n{str(n).replace(".", "p")}'
            if not (directory/'COMPLETE').exists():
                raise RuntimeError(f'Incomplete twist calculation: {directory}')
            meta = json.loads((directory/'metadata.json').read_text())
            ba = next(r for r in BA[n] if r['U_over_t'] == U)
            row = {'U_over_t': U, 'target_density': n, 'actual_density': meta['actual_density_Ne_over_L'],
                   'mu': meta['mu_selection']['mu'], 'q_min': meta['q_min'], 'eta': args.eta,
                   'twist_count': meta['twist_count'], 'directory': str(directory.resolve())}
            assert abs(row['actual_density']-n) < 1e-10
            for block in meta['twists']:
                population_rows.append({'U_over_t': U, 'target_density': n, 'twist_index': block['twist_index'],
                                        'theta_boundary': block['theta_boundary'], 'N': block['N'],
                                        'number_support': json.dumps(block.get('N_support',{str(block['N']):1.})),
                                        'ground_dimension': block['ground_dimension'],
                                        'magnetizations_in_ground_trace': ';'.join(str((c['Nup']-c['Ndn'])/2) for c in block['components'])})
            path = directory/'averaged_qmin_poles.npz'
            hashes[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path, allow_pickle=False) as data:
                save = {'omega': np.linspace(0, 5, 5001), 'eta': args.eta}
                for ch in ['charge', 'spin']:
                    E, W = data[f'{ch}_energies'], data[f'{ch}_weights']; cache[(U, n, ch)] = E.copy(), W.copy()
                    found = peaks(E, W, args.eta); first, strong = found[0], max(found, key=lambda p: p['height'])
                    row.update({f'{ch}_first_omega': first['omega'], f'{ch}_dominant_omega': strong['omega'],
                                f'{ch}_first_velocity': first['omega']/row['q_min'],
                                f'{ch}_dominant_velocity': strong['omega']/row['q_min'],
                                f'{ch}_BA_velocity': ba[f'v_{ch}'],
                                f'{ch}_dominant_relative_error': strong['omega']/row['q_min']/ba[f'v_{ch}']-1})
                    save[f'S_{ch}'] = spectrum(save['omega'], E, W, args.eta)
                    for eta in [.05, .1, .2, .3]:
                        found = peaks(E, W, eta); strong = max(found, key=lambda p: p['height'])
                        sensitivity.append({'U_over_t': U, 'density': n, 'channel': ch, 'eta': eta,
                                            'first_omega': found[0]['omega'], 'dominant_omega': strong['omega'],
                                            'first_velocity': found[0]['omega']/row['q_min'],
                                            'dominant_velocity': strong['omega']/row['q_min']})
                np.savez_compressed(directory/'averaged_qmin_spectra.npz', **save)
            rows.append(row)
            print(f'U={U:g}, n={n:g}: charge={row["charge_dominant_velocity"]:.6f} (BA {row["charge_BA_velocity"]:.6f}), '
                  f'spin={row["spin_dominant_velocity"]:.6f} (BA {row["spin_BA_velocity"]:.6f})', flush=True)
    write_csv(output/'velocity_comparison.csv', rows)
    write_csv(output/'broadening_sensitivity.csv', sensitivity)
    write_csv(output/'selected_ground_states.csv', population_rows)
    plots(rows, BA, cache, output, args.eta, config['twist_count'])
    metadata = {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'input_pole_sha256': hashes, 'eta': args.eta, 'search_window': [0, 5],
                'estimator': 'peak of averaged positive-frequency physical density spectrum divided by q_min',
                'peak_policy': 'first and strongest local maxima reported separately; height>=0.1%, prominence>=1%',
                'caveat': 'Finite-q estimates, not a fitted long-wavelength velocity or finite-size extrapolation'}
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    shutil.copy2(__file__, output/'postprocessor_source_snapshot.py')
    print('Analysis output:', output, flush=True)


if __name__ == '__main__':
    main()
