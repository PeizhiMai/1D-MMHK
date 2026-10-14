#!/usr/bin/env python3
"""Exact T=0 density staircase from the existing twist-resolved E0(N) tables.

No response or Hamiltonian solver is run. Steps are lower-convex-envelope
chemical potentials, not an interpolated root or a coarse sampled mu grid.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mmhk_density_workflow import choose_mu, number_minimizers
from extract_zero_twist_velocities import write_csv

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=HERE/'results/twist20_common_mu_velocity_scan')
    args = p.parse_args()
    output = args.input/'analysis'; output.mkdir(exist_ok=True)
    config = json.loads((args.input/'config.json').read_text()); L = config['L']
    all_rows, stairs, hashes, jumps = [], {}, {}, []
    for U in config['U_list']:
        file = args.input/f'U{U:g}'/'ground_energy_table.npz'
        hashes[str(file.resolve())] = hashlib.sha256(file.read_bytes()).hexdigest()
        with np.load(file, allow_pickle=False) as data:
            E, weights, k = data['E0_by_N'], data['weights'], data['k']
        selection = choose_mu(E, weights, .4, L)
        ps = selection['plateaus']; stairs[U] = ps
        for row in ps:
            all_rows.append({'U_over_t': U, 'mu_left': row['mu_left'], 'mu_right': row['mu_right'],
                             'physical_density_Ne_over_L': 2*row['filling'],
                             'N_per_twist': ';'.join(map(str, row['numbers']))})
        for lower, upper in zip(ps[:-1], ps[1:]):
            nlo, nhi = 2*lower['filling'], 2*upper['filling']
            if nlo+1e-10 < .9 < nhi-1e-10:
                critical_mu = .5*(lower['mu_right']+upper['mu_left'])
                assert abs(lower['mu_right']-upper['mu_left']) < 1e-8
                minimizers = number_minimizers(E, critical_mu)
                switchers = []
                for i, (a, b) in enumerate(zip(lower['numbers'], upper['numbers'])):
                    if a != b:
                        assert a in minimizers[i] and b in minimizers[i]
                        switchers.append({'twist_index': i, 'theta_boundary_over_pi': float(-L*k[i]/np.pi),
                                          'N_below': a, 'N_above': b, 'degenerate_numbers_at_crossing': minimizers[i]})
                left = float(weights@np.argmin(E-(critical_mu-1e-7)*np.arange(2*L+1), axis=1)/L)
                right = float(weights@np.argmin(E-(critical_mu+1e-7)*np.arange(2*L+1), axis=1)/L)
                np.testing.assert_allclose([left, right], [nlo, nhi], atol=1e-10)
                jumps.append({'U_over_t': U, 'mu_critical': critical_mu, 'density_below': nlo,
                              'density_above': nhi, 'jump_size': nhi-nlo, 'switching_twists': switchers,
                              'interpretation': 'Finite-cluster/grid zero-temperature number step; not evidence of thermodynamic phase separation'})
    write_csv(output/'density_mu_plateaus.csv', all_rows)
    (output/'density_jumps_across_n0p9.json').write_text(json.dumps(jumps, indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.25))
    fig.subplots_adjust(left=.075, right=.98, top=.85, bottom=.25, wspace=.23)
    colors = ['#444444', '#234f91', '#c43c64', '#2c8066', '#b17416']
    for U, color in zip(config['U_list'], colors):
        ps = stairs[U]
        edges = [ps[0]['mu_left']]+[p['mu_right'] for p in ps]
        values = [2*p['filling'] for p in ps]
        axes[0].stairs(values, edges, baseline=None, color=color, lw=1.5, label=rf'$U/t={U:g}$')
    for ax in axes:
        ax.axhline(.9, color='#777777', ls='--', lw=1)
        ax.grid(alpha=.16); ax.spines[['top', 'right']].set_visible(False)
        ax.set(xlabel=r'$\mu/t$', ylabel=r'$n=N_e/L$')
    axes[0].axhline(.8, color='#777777', ls=':', lw=1)
    axes[0].set(title='Twenty-twist density staircase', xlim=(-2.1, 3.5), ylim=(0, 1.02))
    axes[0].legend(fontsize=8, loc='lower right', frameon=False)
    jump = next(j for j in jumps if j['U_over_t'] == 7)
    U = jump['U_over_t']; critical = jump['mu_critical']; ps = stairs[U]
    edges = [ps[0]['mu_left']]+[p['mu_right'] for p in ps]
    values = [2*p['filling'] for p in ps]
    axes[1].stairs(values, edges, baseline=None, color='#c43c64', lw=2)
    axes[1].axvline(critical, color='#555555', ls=':', lw=1)
    axes[1].plot([critical, critical], [jump['density_below'], jump['density_above']], 'o',
                 mfc='white', mec='#c43c64', ms=5)
    axes[1].set(title=r'$U/t=7$: jump across $n=0.9$', xlim=(critical-.015, critical+.009), ylim=(.878, .925))
    axes[1].annotate(r'$0.89\ \longrightarrow\ 0.91$', xy=(critical, .902),
                     xytext=(critical-.013, .917), fontsize=10,
                     arrowprops={'arrowstyle': '->', 'color': '#444444'})
    axes[1].text(.03, .045, rf'$\mu_c/t={critical:.9f}$'+'\n'+r'Twists $\theta=\pm0.9\pi$: $N=8\to10$',
                 transform=axes[1].transAxes, fontsize=8.5, va='bottom')
    fig.suptitle(r'$n(\mu)=\frac{1}{20L}\sum_\theta N_0(\theta;\mu)$, $L=10$, $T=0$', fontsize=13, y=.97)
    fig.text(.5, .06, r'No stable $\mu$ interval has $n=0.9$ at $U/t=7$; at $\mu_c$ the density depends on the degenerate-sector choice or mixture.'
             '\nThese are finite-cluster/twist-grid number steps, not a thermodynamic phase-separation diagnosis.',
             ha='center', fontsize=8.5)
    fig.savefig(output/'twist20_density_vs_mu.png', dpi=180)
    fig.savefig(output/'twist20_density_vs_mu.svg')
    plt.close(fig)
    metadata = {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'ground_energy_table_sha256': hashes, 'L': L, 'twist_count': config['twist_count'],
                'physical_density_convention': 'Ne/L, not Ne/(2L)',
                'method': 'exact lower convex envelopes of E0(theta,N)-mu*N; full spin ground energy',
                'additional_solver_runs': 0, 'jump_at_crossing_is_multivalued': True}
    (output/'density_mu_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(jumps, indent=2), flush=True)
    print('Density plot:', output/'twist20_density_vs_mu.png', flush=True)


if __name__ == '__main__':
    main()
