#!/usr/bin/env python3
"""Inspect density attainability from convex-envelope E0(N,theta) tables."""
import json
from pathlib import Path

import numpy as np

from mmhk_density_workflow import choose_mu, number_minimizers
from twist_averaged_velocity_scan import digest
from extract_zero_twist_velocities import write_csv

HERE = Path(__file__).resolve().parent
ROOT = HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'


def main():
    output = ROOT/'density_diagnostics'; output.mkdir(exist_ok=True)
    rows, jumps, hashes, stairs = [], [], {}, {}
    for U in range(13):
        file = ROOT/f'U{U}'/'ground_energy_table.npz'
        if not file.exists():
            continue
        hashes[str(file.resolve())] = digest(file)
        with np.load(file) as d:
            E, w, k = d['E0_by_N'], d['weights'], d['k']
        ps = choose_mu(E, w, .4, 10)['plateaus']; stairs[U] = ps
        for n in [.7, .8, .9]:
            s = choose_mu(E, w, n/2, 10)
            lo = max((p for p in ps if 2*p['filling'] <= n+1e-10), key=lambda p: p['filling'])
            hi = min((p for p in ps if 2*p['filling'] >= n-1e-10), key=lambda p: p['filling'])
            free = ROOT/'U0'/f'free_selection_n{str(n).replace(".", "p")}.json'
            projector = json.loads(free.read_text()) if U == 0 else None
            rows.append(dict(U=U, target_n=n, exact_stable_plateau=s['exact_target_attained'],
                             full_free_projector_target_attained=bool(projector and projector['exact_target_attained']),
                             nearest_stable_n=2*s['filling'], delta_n=2*s['filling']-n,
                             density_below=2*lo['filling'], density_above=2*hi['filling'],
                             nearest_mu=s['mu'], mu_left=s['mu_left'], mu_right=s['mu_right'],
                             plateau_width=s['mu_right']-s['mu_left']))
            if 2*lo['filling']+1e-10 < n < 2*hi['filling']-1e-10:
                assert abs(lo['mu_right']-hi['mu_left']) < 1e-8
                critical = (lo['mu_right']+hi['mu_left'])/2
                minima = number_minimizers(E, critical)
                changed = []
                for i, (a, b) in enumerate(zip(lo['numbers'], hi['numbers'])):
                    if a == b:
                        continue
                    assert a in minima[i] and b in minima[i]
                    grand = E[i]-critical*np.arange(21)
                    changed.append(dict(twist_index=i, theta_over_pi=float(-10*k[i]/np.pi),
                                        N_below=a, N_above=b, numbers_degenerate=minima[i],
                                        grand_energy_N9_above_min=float(grand[9]-grand.min())))
                jumps.append(dict(U=U, target_n=n, density_below=2*lo['filling'], density_above=2*hi['filling'],
                                  mu_critical=critical, switching_twists=changed))
    write_csv(output/'density_attainability.csv', rows)
    (output/'density_jumps.json').write_text(json.dumps(jumps, indent=2)+'\n')
    (output/'metadata.json').write_text(json.dumps(dict(source_sha256=digest(__file__), energy_table_sha256=hashes,
        method='exact lower convex envelope, stable common-mu plateaus', physical_density='Ne/L',
        student_solver_executed=False, free_projector_is_not_a_stable_plateau=True), indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    U = 7; ps = stairs[U]
    edges = [ps[0]['mu_left']]+[p['mu_right'] for p in ps]
    values = [2*p['filling'] for p in ps]
    fig, axes = plt.subplots(1, 3, figsize=(10.7, 3.7))
    fig.subplots_adjust(left=.075, right=.985, top=.79, bottom=.23, wspace=.30)
    for ax, n in zip(axes, [.7, .8, .9]):
        entry = next(r for r in rows if r['U'] == U and r['target_n'] == n)
        s = next(p for p in ps if abs(p['mu_left']-entry['mu_left']) < 1e-12)
        width = max(.005, s['mu_right']-s['mu_left'])
        ax.stairs(values, edges, baseline=None, color='#234f91', lw=1.7)
        ax.axhline(n, color='#777777', ls='--', lw=1)
        ax.plot(entry['nearest_mu'], entry['nearest_stable_n'], 'o', color='#234f91', ms=5)
        ax.set(title=rf'Target $n={n:g}$', xlabel=r'$\mu/t$', ylabel=r'$\bar n$',
               ylim=(n-.018, n+.018), xlim=(entry['mu_left']-width*.35, entry['mu_right']+width*.35))
        ax.grid(alpha=.15); ax.spines[['top', 'right']].set_visible(False)
        if not entry['exact_stable_plateau']:
            jump = next(j for j in jumps if j['U'] == U and j['target_n'] == n)
            critical = jump['mu_critical']
            ax.set_xlim(critical-.005, critical+.003)
            ax.plot([critical]*2, [jump['density_below'], jump['density_above']], 'o', mfc='white', color='#c43c64')
            ax.text(.05, .87, f'{jump["density_below"]:.3f} → {jump["density_above"]:.3f}', transform=ax.transAxes, fontsize=10)
    fig.suptitle(r'$U/t=7$: 50-twist density staircases, $L=10$', fontsize=13, y=.96)
    fig.text(.5, .055, r'$n=0.7$ and $0.8$ are stable plateaus. No stable interval gives exactly $n=0.9$; its adjacent plateaus are 0.896 and 0.904.',
             ha='center', fontsize=8.5)
    fig.savefig(output/'U7_density_vs_mu.png', dpi=180); fig.savefig(output/'U7_density_vs_mu.svg'); plt.close(fig)
    print(json.dumps(jumps, indent=2)); print(output)


if __name__ == '__main__':
    main()
