#!/usr/bin/env python3
"""Configurable q_min common-mu 1D MMHK response. Plans by default.

Reuse the validated project qbasis solver; never execute student code. Supports
even L>=4, one zero twist or an even periodic grid >=4, and physical 0<n<1.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import itertools
import json
from math import comb, ceil
import os
from pathlib import Path
import shutil
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
import numpy as np

PROJECT = Path('/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d')


def free_case(root, L, k, density, config, merge_poles, write_csv, tag, energies_only=False):
    """Exact full grand ground projector, including degenerate Fermi shells."""
    nk = len(k); K = 2*np.pi*np.arange(L)/L
    eps = -2*np.cos(K[None, :]+k[:, None])
    levels = np.sort(np.repeat(eps, 2, axis=1), axis=1)
    E0 = np.column_stack([np.zeros(nk), np.cumsum(levels, axis=1)])
    np.savez_compressed(root/'ground_energy_table.npz', k=k, weights=np.ones(nk)/nk, E0_by_N=E0)
    target = density*L*nk
    if abs(target-round(target)) > 1e-9 or not 0 < round(target) < 2*L*nk:
        if energies_only:
            from mmhk_density_workflow import choose_mu
            selection = choose_mu(E0, np.ones(nk)/nk, density/2, L)
            selection.update(target_density_Ne_over_L=density,
                             actual_density_Ne_over_L=2*selection['filling'])
            (root/f'free_selection_n{tag(density)}.json').write_text(json.dumps(selection, indent=2)+'\n')
            return selection
        raise RuntimeError('U=0 target is not attainable by a stable plateau or full grand projector on this grid; choose an explicit ensemble or another grid')
    total = round(target); flat = np.sort(levels.ravel())
    mu = float((flat[total-1]+flat[total])/2)
    f = np.where(eps < mu-1e-10, 1., np.where(eps > mu+1e-10, 0., .5))
    mean_N = 2*f.sum(axis=1); actual = float(mean_N.mean()/L)
    selection = {'mu': mu, 'target_density_Ne_over_L': density,
                 'actual_density_Ne_over_L': actual, 'filling': actual/2,
                 'exact_target_attained': abs(actual-density) < 1e-10,
                 'density_below_crossing': float(2*(eps < mu-1e-10).sum()/(L*nk)),
                 'density_above_crossing': float(2*(eps <= mu+1e-10).sum()/(L*nk)),
                 'policy': 'stable free plateau or complete zero-temperature grand ground projector'}
    (root/f'free_selection_n{tag(density)}.json').write_text(json.dumps(selection, indent=2)+'\n')
    if energies_only:
        return selection
    if not selection['exact_target_attained']:
        raise RuntimeError('U=0 full grand projector does not give the requested density; no arbitrary shell occupation was chosen')
    blocks, E, W = [], [], []; elastic = 0.
    for i in range(nk):
        filled = [(s, j) for s in [0, 1] for j in range(L) if f[i, j] == 1]
        shell = [(s, j) for s in [0, 1] for j in range(L) if f[i, j] == .5]
        groups, histogram = {}, {}
        for bits in itertools.product([0, 1], repeat=len(shell)):
            occupied = filled+[mode for bit, mode in zip(bits, shell) if bit]
            Nu = sum(s == 0 for s, _ in occupied); Nd = len(occupied)-Nu
            Q = (-sum(j for _, j in occupied)) % L
            groups.setdefault((Nu, Nd), []).append(Q)
            histogram[Nu+Nd] = histogram.get(Nu+Nd, 0)+1
        dimension = 2**len(shell); components = []
        for (Nu, Nd), momenta in sorted(groups.items()):
            if Nu < Nd:
                continue
            mult = 1 if Nu == Nd else 2
            components.append({'N': Nu+Nd, 'Nup': Nu, 'Ndn': Nd,
                               'ground_dimension': len(momenta), 'spin_flip_multiplicity': mult,
                               'probability': mult*len(momenta)/dimension,
                               'ground_states': [{'Q': Q, 'id': j} for j, Q in enumerate(momenta)]})
        assert sum(c['ground_dimension']*c['spin_flip_multiplicity'] for c in components) == dimension
        blocks.append({'twist_index': i, 'k': float(k[i]), 'phase_per_bond': float(-k[i]),
                       'theta_boundary': float(-L*k[i]), 'N': float(mean_N[i]),
                       'N_support': {str(N): count/dimension for N, count in histogram.items()},
                       'ground_dimension': dimension, 'components': components})
        for dest in range(L):
            src = (dest+1) % L; probability = f[i, src]*(1-f[i, dest])
            delta = float(eps[i, dest]-eps[i, src])
            if probability > 0 and delta > 1e-8:
                E.append(delta); W.append(2*probability/(L*nk))
            elif probability > 0 and abs(delta) <= 1e-8:
                elastic += 2*probability/(L*nk)
            # delta=0 is elastic and stored separately from inelastic poles.
    E, W = merge_poles(np.asarray(E), np.asarray(W))
    assert not len(E) or E.max() <= 4*np.sin(np.pi/L)+1e-10
    out = root/f'n{tag(density)}'; out.mkdir(exist_ok=True)
    np.savez_compressed(out/'averaged_qmin_poles.npz', k=k, theta_boundary=-L*k,
                        N_per_twist=mean_N, twist_weights=np.ones(nk)/nk, q_min=2*np.pi/L,
                        U=0., mu=mu, actual_density=actual,
                        charge_energies=E, charge_weights=W, spin_energies=E, spin_weights=W,
                        elastic_charge=elastic, elastic_spin=elastic)
    write_csv(out/'averaged_qmin_poles.csv', [{'channel': ch, 'omega': float(e), 'weight': float(w)}
                                           for ch in ['charge', 'spin'] for e, w in zip(E, W)])
    meta = {'L': L, 'U': 0., 'target_density_Ne_over_L': density, 'actual_density_Ne_over_L': actual,
            'native_filling_Ne_over_2L': actual/2, 'twist_count': nk, 'mu_selection': selection,
            'q_min': 2*np.pi/L, 'twists': blocks, 'response_method': 'exact free particle-hole sum',
            'ensemble': 'complete grand ground projector; all number/spin states at the Fermi level',
            'observable': 'coherent physical charge and longitudinal 2Sz densities',
            'frequency_origin': 'neutral Em-Eg, no mu shift', 'student_solver_executed': False,
            'elastic_weight_charge_and_spin': elastic,
            'solver_source_sha256': config['solver_source_sha256'],
            'solver_executable_sha256': config['solver_executable_sha256']}
    (out/'metadata.json').write_text(json.dumps(meta, indent=2)+'\n')
    np.savez_compressed(out/'exact_free_reference.npz', epsilon=eps, f=f, k=k, mu=mu, E0_by_N=E0)
    (out/'COMPLETE').write_text('success\n')
    print(f'COMPLETE U=0 L={L} n={actual:g} twists={nk}: charge=spin', flush=True)
    return selection


def analyze(a, config, spectrum, peaks, write_csv, tag, digest):
    from bethe_reference import load_ba
    solve = load_ba(a.project)
    out = a.output/f'skill_analysis_eta{tag(a.eta)}'; out.mkdir(exist_ok=True)
    analysis_config = {key: getattr(a, key) for key in ['eta', 'omega_max', 'points', 'peak_cutoff', 'nodes']}
    hashes = {}
    for U in config['U_list']:
        for n in config['physical_densities']:
            directory = a.output/f'U{U:g}'/f'n{tag(n)}'
            if not (directory/'COMPLETE').exists():
                raise RuntimeError(f'Incomplete case: {directory}')
            for name in ['metadata.json', 'averaged_qmin_poles.npz']:
                hashes[str((directory/name).resolve())] = digest(directory/name)
    analysis_config.update(input_sha256=hashes, script_sha256=digest(__file__),
                           peak_helper_sha256=digest(a.project/'response/postprocess_twist_velocity_scan.py'),
                           BA_source_sha256=digest(a.project/'response/essler_bethe_ansatz/hubbard_velocities.py'))
    manifest = out/'config.json'
    if manifest.exists() and json.loads(manifest.read_text()) != analysis_config:
        raise RuntimeError('Different analysis settings or input data; choose a fresh analysis destination')
    manifest.write_text(json.dumps(analysis_config, indent=2)+'\n')
    omega = np.linspace(0, a.omega_max, a.points); rows, sensitivity, cache = [], [], {}
    populations = []
    for U in config['U_list']:
        for n in config['physical_densities']:
            directory = a.output/f'U{U:g}'/f'n{tag(n)}'
            m = json.loads((directory/'metadata.json').read_text())
            if m['L'] != config['L'] or abs(m['actual_density_Ne_over_L']-n) > 1e-10:
                raise RuntimeError('Wrong mixing size or actual density in analysis input')
            ba = solve(U, n, a.nodes)
            row = dict(U_over_t=U, L=m['L'], target_density=n, actual_density=m['actual_density_Ne_over_L'],
                       mu=m['mu_selection']['mu'], q_min=m['q_min'], eta=a.eta,
                       twist_count=m['twist_count'])
            with np.load(directory/'averaged_qmin_poles.npz', allow_pickle=False) as d:
                arrays = dict(omega=omega, eta=a.eta)
                for ch in ['charge', 'spin']:
                    E, W = d[f'{ch}_energies'], d[f'{ch}_weights']
                    arrays[f'S_{ch}'] = spectrum(omega, E, W, a.eta)
                    cache[(U, n, ch)] = arrays[f'S_{ch}']
                    significant = W > 1e-6*W.sum()
                    lowest = float(E[significant].min()) if np.any(significant) else None
                    row.update({f'{ch}_BA_velocity': ba[f'v_{ch}'], f'{ch}_lowest_significant_pole': lowest,
                                f'{ch}_integrated_inelastic_weight': float(W.sum())})
                    for eta in sorted(set([.05, .1, .2, .3, a.eta])):
                        try:
                            found = peaks(E, W, eta, a.peak_cutoff)
                        except RuntimeError as error:
                            if str(error) != 'No resolved positive-frequency local maximum':
                                raise
                            found = []
                        first = found[0]['omega'] if found else None
                        dominant = max(found, key=lambda r: r['height'])['omega'] if found else None
                        sensitivity.append(dict(U=U, density=n, channel=ch, eta=eta,
                                                first_omega=first, dominant_omega=dominant))
                        if eta == a.eta:
                            row.update({f'{ch}_first_omega': first, f'{ch}_dominant_omega': dominant,
                                        f'{ch}_first_velocity': None if first is None else first/m['q_min'],
                                        f'{ch}_dominant_velocity': None if dominant is None else dominant/m['q_min']})
                np.savez_compressed(out/f'U{U:g}_n{tag(n)}_spectra.npz', **arrays)
            for block in m['twists']:
                populations.append(dict(U=U, density=n, theta=block['theta_boundary'], mean_N=block['N'],
                                        N_support=json.dumps(block.get('N_support', {str(block['N']): 1.})),
                                        ground_dimension=block['ground_dimension']))
            rows.append(row)
    write_csv(out/'velocity_comparison.csv', rows); write_csv(out/'broadening_sensitivity.csv', sensitivity)
    write_csv(out/'selected_ground_states.csv', populations)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for n in config['physical_densities']:
        selected = [r for r in rows if r['target_density'] == n]; count = len(selected)
        fig, axes = plt.subplots(ceil(count/4), min(4, count), squeeze=False,
                                 figsize=(min(4, count)*2.6, ceil(count/4)*2.1+.6), layout='constrained')
        for ax, row in zip(axes.flat, selected):
            for ch, color, style in [('charge', '#c43c64', '-'), ('spin', '#234f91', '--')]:
                ax.plot(omega, cache[(row['U_over_t'], n, ch)], color=color, ls=style, label=ch)
            ax.set(title=f'U/t={row["U_over_t"]:g}', xlabel=r'$\omega/t$', ylabel=r'$S(q_1,\omega)$', xlim=(0, a.omega_max))
        axes.flat[0].legend(frameon=False, fontsize=8)
        for ax in list(axes.flat)[count:]:
            ax.set_visible(False)
        fig.suptitle(f'L={config["L"]}, n={n:g}, {config["twist_count"]} twists, eta/t={a.eta:g}', fontsize=11)
        fig.savefig(out/f'n{tag(n)}_spectra.png', dpi=160); plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.1), layout='constrained')
        for ax, ch, color in zip(axes, ['charge', 'spin'], ['#c43c64', '#234f91']):
            ax.plot([r['U_over_t'] for r in selected], [r[f'{ch}_BA_velocity'] for r in selected],
                    color=color, label='Thermodynamic BA')
            ax.plot([r['U_over_t'] for r in selected], [r[f'{ch}_dominant_velocity'] for r in selected],
                    'o', color='black', label='Finite-q dominant peak/q')
            ax.set(title=ch, xlabel='U/t', ylabel=r'$v/(ta_0/\hbar)$'); ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f'L={config["L"]}, n={n:g}: finite-q/eta-dependent estimates, not q→0 fits', fontsize=10)
        fig.savefig(out/f'n{tag(n)}_velocity_comparison.png', dpi=160); plt.close(fig)
    (out/'metadata.json').write_text(json.dumps({
        'normalization': 'absolute intensity; charge and longitudinal 2Sz; source 1/sqrt(L)',
        'peak_policy': '0<omega<cutoff, height>=0.1% maximum, prominence>=1%; first and strongest separately',
        'lowest_pole_policy': 'residue > 1e-6 of total inelastic weight',
        'velocity_caveat': 'Finite-q, eta-dependent peak estimates, not thermodynamic slopes',
        'student_solver_executed': False
    }, indent=2)+'\n')
    print('Analysis:', out.resolve())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, default=PROJECT)
    p.add_argument('--L', type=int, default=10)
    p.add_argument('--densities', type=float, nargs='+', default=[.8])
    p.add_argument('--couplings', type=float, nargs='+', default=list(range(13)))
    p.add_argument('--twists', type=int, default=20, help='1 for zero twist, or even >=4')
    p.add_argument('--depth', type=int, default=300)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--eta', type=float, default=.2)
    p.add_argument('--omega-max', type=float, default=5)
    p.add_argument('--points', type=int, default=5001)
    p.add_argument('--peak-cutoff', type=float, default=5)
    p.add_argument('--nodes', type=int, default=400, help='BA quadrature nodes')
    p.add_argument('--output', type=Path, required=True)
    modes = p.add_mutually_exclusive_group()
    modes.add_argument('--execute', action='store_true')
    modes.add_argument('--analyze-only', action='store_true')
    p.add_argument('--energies-only', action='store_true')
    a = p.parse_args(); a.project = a.project.resolve(); a.output = a.output.resolve()
    if (a.L < 4 or a.L % 2 or (a.twists != 1 and (a.twists < 4 or a.twists % 2))
            or min(a.workers, a.threads) < 1 or a.depth < 3 or a.nodes < 50 or a.points < 3
            or min(a.eta, a.omega_max, a.peak_cutoff) <= 0
            or not np.isfinite([*a.densities, *a.couplings, a.eta, a.omega_max, a.peak_cutoff]).all()
            or any(not 0 < n < 1 for n in a.densities) or any(U < 0 for U in a.couplings)):
        p.error('Require even L>=4, metallic 0<n<1, U>=0, and valid resources/grid')
    if a.energies_only and a.analyze_only:
        p.error('Cannot analyze an energies-only request')
    a.couplings = sorted(set(a.couplings)); a.densities = sorted(set(a.densities))
    density_names = [f'{n:.8f}'.rstrip('0').rstrip('.') for n in a.densities]
    if len(set(density_names)) != len(density_names) or len({f'{U:g}' for U in a.couplings}) != len(a.couplings):
        p.error('Parameter values collide in output names; choose clearly separated values')
    plan = dict(L=a.L, physical_densities=a.densities, U_list=a.couplings, twist_count=a.twists,
                q_min=2*np.pi/a.L, workers=a.workers, OMP_threads_per_worker=a.threads,
                maximum_canonical_dimension=comb(a.L, a.L//2)**2,
                estimated_minimum_energy_sector_solves_per_U=(a.L+1)*(1 if a.twists == 1 else a.twists//2+1),
                ensemble='actual complete ground trace; common mu; reject unattainable target',
                output=str(a.output), solver='fast qbasis; no student solver')
    if not a.execute and not a.analyze_only:
        print(json.dumps(plan, indent=2)); print('Plan only. Add --execute to calculate.'); return
    sys.path.insert(0, str(a.project/'response'))
    from twist_averaged_velocity_scan import Runs, energy_table, ground_components, save_density, digest, tag
    from mmhk_density_workflow import k_grid, choose_mu
    from extract_zero_twist_velocities import merge_poles, write_csv
    from postprocess_twist_velocity_scan import spectrum, peaks
    if a.analyze_only:
        config = json.loads((a.output/'config.json').read_text())
        analyze(a, config, spectrum, peaks, write_csv, tag, digest); return
    source = a.project/'create_code/chain_Hubbard_tbc_response.cc'
    exe = a.project/'response/bin/chain_Hubbard_tbc_response.x'
    if not exe.exists():
        p.error('Build response solver first with the project response/build_response.sh')
    helpers = ['twist_averaged_velocity_scan.py', 'mmhk_density_workflow.py', 'postprocess_response.py',
               'extract_zero_twist_velocities.py']
    config = dict(L=a.L, U_list=a.couplings, physical_densities=a.densities, twist_count=a.twists,
                  grid='periodic', depth=a.depth, q_indices=[1], solver_source_sha256=digest(source),
                  solver_executable_sha256=digest(exe), runner_sha256=digest(__file__),
                  helper_sha256={name: digest(a.project/'response'/name) for name in helpers})
    a.output.mkdir(parents=True, exist_ok=True); manifest = a.output/'config.json'
    if manifest.exists() and json.loads(manifest.read_text()) != config:
        raise RuntimeError('Different campaign settings or source hashes; use a fresh output root')
    if not manifest.exists():
        shutil.copy2(source, a.output/'solver_source_snapshot.cc')
        shutil.copy2(exe, a.output/'solver_executable_snapshot.x')
        shutil.copy2(__file__, a.output/'runner_source_snapshot.py')
        for name in helpers:
            shutil.copy2(a.project/'response'/name, a.output/f'{name[:-3]}_source_snapshot.py')
        manifest.write_text(json.dumps(config, indent=2)+'\n')
    frozen = a.output/'solver_executable_snapshot.x'
    if digest(frozen) != config['solver_executable_sha256']:
        raise RuntimeError('Frozen executable hash mismatch')
    k, weights = k_grid(a.L, a.twists, 'periodic'); summary = []
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for U in a.couplings:
            root = a.output/f'U{U:g}'; root.mkdir(exist_ok=True)
            if U == 0:
                for density in a.densities:
                    selected = free_case(root, a.L, k, density, config, merge_poles, write_csv, tag, a.energies_only)
                    summary.append(dict(U=U, density=density, mu=selected['mu'],
                                        actual_density=selected['actual_density_Ne_over_L'],
                                        exact_target=selected['exact_target_attained']))
                write_csv(a.output/'density_mu_summary.csv', summary)
                continue
            runs = Runs(root, U, a.L, a.depth, a.threads, frozen,
                        config['solver_source_sha256'], config['solver_executable_sha256'])
            if a.twists == 1:
                E = np.array([[runs.run(0, 0., (N+1)//2, N//2)[1]['E0'] for N in range(2*a.L+1)]])
                np.savez_compressed(root/'ground_energy_table.npz', k=k, weights=weights, E0_by_N=E)
            else:
                E = energy_table(runs, k, pool)
            selections = {}
            for density in a.densities:
                selected = choose_mu(E, weights, density/2, a.L)
                selected.update(target_density_Ne_over_L=density, actual_density_Ne_over_L=2*selected['filling'])
                selections[density] = selected
                summary.append(dict(U=U, density=density, mu=selected['mu'], actual_density=2*selected['filling'],
                                    exact_target=selected['exact_target_attained']))
            (root/'density_selections.json').write_text(json.dumps(selections, indent=2)+'\n')
            write_csv(a.output/'density_mu_summary.csv', summary)
            if a.energies_only:
                continue
            if not all(s['exact_target_attained'] for s in selections.values()):
                raise RuntimeError(f'U={U:g}: target lies in a density step. Inspect density_selections.json plateaus; no replacement density or crossing mixture used')
            tasks = sorted({(i, N) for s in selections.values() for i, N in enumerate(s['numbers'])})
            futures = {pool.submit(ground_components, runs, i, k[i], N, E[i, N]): (i, N) for i, N in tasks}
            components = {futures[f]: f.result() for f in as_completed(futures)}
            for density, selected in selections.items():
                save_density(runs, density, selected, k, components)
    if not a.energies_only:
        analyze(a, config, spectrum, peaks, write_csv, tag, digest)


if __name__ == '__main__':
    main()
