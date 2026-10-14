#!/usr/bin/env python3
"""Common-mu T=0 physical Hubbard density responses on a uniform twist grid.

Physical n=Ne/L. Energy scans use the minimum-|Sz| sector, which contains all
SU(2) multiplets; responses trace ALL magnetizations in the actual ground
subspace, not a fixed Sz state. Neutral excitation frequencies need no mu shift.
This never imports or runs the student's solver.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time

import numpy as np

from mmhk_density_workflow import choose_mu, number_minimizers, k_grid, DEGENERACY_TOL
from postprocess_response import read_fraction, poles_from_fraction
from extract_zero_twist_velocities import merge_poles, write_csv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tag(x):
    return f'{x:.8f}'.rstrip('0').rstrip('.').replace('.', 'p')


class Runs:
    def __init__(self, root, U, L, depth, threads, exe, source_hash, exe_hash):
        self.root, self.U, self.L = Path(root), U, L
        self.depth, self.exe = depth, Path(exe)
        self.env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS='1')
        self.source_hash, self.exe_hash = source_hash, exe_hash
        self.lock = threading.Lock(); self.new_calls = 0; self.seconds = 0.

    def run(self, ik, k, Nu, Nd, ground_only=True):
        kind = 'energy_scan' if ground_only else 'qmin_responses'
        path = self.root/kind/f'twist{ik:03d}'/f'Nup{Nu}_Ndn{Nd}'
        command = [str(self.exe), '--L', str(self.L), '--U', str(self.U),
                   '--nup', str(Nu), '--ndn', str(Nd), '--phase', repr(float(-k)),
                   '--depth', str(self.depth), '--observable', 'density', '--output', str(path.resolve())]
        command += ['--ground-only'] if ground_only else ['--q', '1']
        identity = {'command': command, 'source_sha256': self.source_hash,
                    'executable_sha256': self.exe_hash}
        lineage = path/'lineage.json'
        if (path/'COMPLETE').exists():
            if not lineage.exists() or json.loads(lineage.read_text()) != identity:
                raise RuntimeError(f'Incompatible cached run: {path}')
        else:
            if (path/'metadata.json').exists():
                raise RuntimeError(f'Interrupted run requires inspection: {path}')
            path.mkdir(parents=True, exist_ok=True)
            lineage.write_text(json.dumps(identity, indent=2)+'\n')
            start = time.monotonic()
            with (path/'run.log').open('w') as log:
                subprocess.run(command, env=self.env, check=True, stdout=log, stderr=subprocess.STDOUT)
            with self.lock:
                self.new_calls += 1; self.seconds += time.monotonic()-start
        m = json.loads((path/'metadata.json').read_text())
        if (m['L'], m['U'], m['Nup'], m['Ndn']) != (self.L, self.U, Nu, Nd):
            raise RuntimeError('Canonical metadata mismatch')
        if abs(m['phase_per_bond']+k) > 1e-12 or not m['ground_scan_all_Q']:
            raise RuntimeError('Wrong twist or restricted ground scan')
        return path, m


def energy_table(runs, k, pool):
    """Exact spectral symmetries, not a truncation of the number scan.

    E(N,theta)=E(N,-theta) by complex conjugation. On an even nearest-neighbor
    ring, E(2L-N,theta)=E(N,theta)+U*(L-N) by particle-hole plus conjugation.
    Responses are still calculated directly at EACH of the twenty twists.
    """
    nk, L = len(k), runs.L
    partner = (-np.arange(nk)) % nk
    representatives = sorted({min(i, int(partner[i])) for i in range(nk)})
    low = np.empty((nk, L+1))
    futures = {pool.submit(runs.run, i, k[i], (N+1)//2, N//2): (i, N)
               for i in representatives for N in range(L+1)}
    for count, future in enumerate(as_completed(futures), 1):
        i, N = futures[future]
        _, m = future.result(); low[i, N] = m['E0']; low[partner[i], N] = m['E0']
        if count % 20 == 0 or count == len(futures):
            print(f'U={runs.U:g}: ground energies {count}/{len(futures)}', flush=True)
    E = np.empty((nk, 2*L+1)); E[:, :L+1] = low
    for N in range(L):
        E[:, 2*L-N] = low[:, N]+runs.U*(L-N)
    # Direct nonzero-twist checks of BOTH shortcuts at every U.
    i = nk//2-2; j = int(partner[i]); N = L-2
    checks = []
    for index, number in [(j, L-1), (i, 2*L-N)]:
        _, m = runs.run(index, k[index], (number+1)//2, number//2)
        error = abs(m['E0']-E[index, number]); assert error < 1e-8
        checks.append({'twist_index': index, 'N': number, 'energy_error': error})
    np.savez_compressed(runs.root/'ground_energy_table.npz', E0_by_N=E, k=k,
                        weights=np.ones(nk)/nk, conjugate_partner=partner,
                        explicitly_scanned_twist_indices=representatives)
    (runs.root/'energy_symmetry_checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    return E


def ground_components(runs, i, k, N, E_min):
    active = []
    for Nu in range((N+1)//2, min(runs.L, N)+1):
        Nd = N-Nu
        if not 0 <= Nd <= runs.L:
            continue
        _, m = runs.run(i, k, Nu, Nd)
        if m['E0'] > E_min+DEGENERACY_TOL:
            break  # Higher-|Sz| spectra are subsets of the lower-|Sz| spectra.
        if abs(m['E0']-E_min) > DEGENERACY_TOL:
            raise RuntimeError('Minimum-|Sz| energy scan missed a lower ground state')
        path, m = runs.run(i, k, Nu, Nd, False)
        if not m.get('ground_subspace_complete') or m['q_indices'] != [1]:
            raise RuntimeError('Uncertified ground subspace or wrong momentum transfer')
        active.append({'N': N, 'Nup': Nu, 'Ndn': Nd,
                       'spin_flip_multiplicity': 1 if Nu == Nd else 2,
                       'ground_dimension': len(m['ground_states']),
                       'directory': str(path.resolve()), 'metadata': m})
    if not active:
        raise RuntimeError('Empty ground subspace')
    return active


def fractions(component, channel, hashes, checks):
    path = Path(component['directory']); m = component['metadata']
    seen = set(); energies, weights = [], []
    for file in sorted(path.glob(f'{channel}_gsQ*_q1.dat')):
        d = read_fraction(file)
        identity = int(d['ground_Q']), int(d.get('ground_id', 0))
        if identity in seen:
            raise RuntimeError('Duplicate ground-state source')
        seen.add(identity)
        assert abs(d['ensemble_weight']-1/len(m['ground_states'])) < 1e-12
        e, w = poles_from_fraction(d)
        assert not np.any((e < -1e-8) & (w > 1e-8))
        weight_error = abs(w.sum()-d['source_weight'])
        moment_error = abs(e@w-d['first_moment'])
        assert weight_error < 1e-9 and moment_error < 1e-8
        checks.append({'file': str(file.resolve()), 'weight_error': float(weight_error),
                       'first_moment_error': float(moment_error)})
        good = (e > 1e-8) & (w > 0)
        energies.extend(e[good]); weights.extend(w[good]*d['ensemble_weight'])
        hashes[str(file.resolve())] = digest(file)
    assert seen == {(g['Q'], g['id']) for g in m['ground_states']}
    return merge_poles(np.asarray(energies), np.asarray(weights))


def save_density(runs, density, selection, k, components):
    output = runs.root/f'n{tag(density)}'; output.mkdir(exist_ok=True)
    positive = {ch: [[], []] for ch in ['charge', 'spin']}
    hashes, checks, blocks, numbers = {}, [], [], []
    nk = len(k)
    for i in range(nk):
        active = components[(i, selection['numbers'][i])]
        dimension = sum(c['spin_flip_multiplicity']*c['ground_dimension'] for c in active)
        block = {'twist_index': i, 'k': float(k[i]), 'phase_per_bond': float(-k[i]),
                 'theta_boundary': float(-runs.L*k[i]), 'N': selection['numbers'][i],
                 'ground_dimension': dimension, 'components': []}
        numbers.append(block['N'])
        for c in active:
            probability = c['spin_flip_multiplicity']*c['ground_dimension']/dimension
            block['components'].append({key: value for key, value in c.items() if key != 'metadata'}
                                       | {'probability': probability, 'ground_states': c['metadata']['ground_states']})
            for ch in positive:
                e, w = fractions(c, ch, hashes, checks)
                positive[ch][0].extend(e); positive[ch][1].extend(w*probability/nk)
        blocks.append(block)
    actual = float(np.mean(numbers)/runs.L)
    assert abs(actual-density) < 1e-10
    arrays = {'k': k, 'theta_boundary': -runs.L*k, 'N_per_twist': np.asarray(numbers),
              'twist_weights': np.ones(nk)/nk, 'q_min': 2*np.pi/runs.L,
              'U': runs.U, 'mu': selection['mu'], 'actual_density': actual}
    pole_rows = []
    for ch, (es, ws) in positive.items():
        e, w = merge_poles(np.asarray(es), np.asarray(ws))
        np.testing.assert_allclose(w.sum(), np.sum(ws), atol=1e-11, rtol=1e-11)
        arrays[f'{ch}_energies'], arrays[f'{ch}_weights'] = e, w
        for energy, weight in zip(e, w):
            pole_rows.append({'channel': ch, 'omega': float(energy), 'weight': float(weight)})
    np.savez_compressed(output/'averaged_qmin_poles.npz', **arrays)
    write_csv(output/'averaged_qmin_poles.csv', pole_rows)
    metadata = {'L': runs.L, 'U': runs.U, 'target_density_Ne_over_L': density,
                'actual_density_Ne_over_L': actual, 'native_filling_Ne_over_2L': actual/2,
                'twist_count': nk, 'grid': f'{nk} distinct uniform boundary twists, periodic grid includes zero',
                'mu_selection': selection, 'ensemble': 'complete zero-temperature ground trace of H(theta)-mu*N',
                'q_min': 2*np.pi/runs.L, 'observable': 'coherent physical charge and longitudinal 2Sz densities',
                'frequency_origin': 'neutral Em-Eg; chemical potential cancels, no frequency shift',
                'normalization': 'source 1/sqrt(L); equal twist weights sum to one',
                'twists': blocks, 'source_moment_checks': checks, 'raw_fraction_sha256': hashes,
                'solver_source_sha256': runs.source_hash, 'solver_executable_sha256': runs.exe_hash,
                'student_solver_executed': False}
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (output/'COMPLETE').write_text('success\n')
    print(f'COMPLETE U={runs.U:g}, n={density:g}: mu={selection["mu"]:.9f}, actual n={actual:.9f}, N(theta)={numbers}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--couplings', nargs='+', type=float, default=[4, 6, 7, 8, 12])
    p.add_argument('--densities', nargs='+', type=float, default=[.8])
    p.add_argument('--twists', type=int, default=20)
    p.add_argument('--depth', type=int, default=300)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--energies-only', action='store_true')
    p.add_argument('--output', type=Path, default=HERE/'results/twist20_common_mu_velocity_scan')
    args = p.parse_args(); L = 10
    if args.twists < 4 or args.twists % 2 or args.workers < 1 or args.threads < 1 or args.depth < 3:
        p.error('Require an even twist count and positive resources')
    if any(not 0 < n < 1 for n in args.densities) or any(U <= 0 for U in args.couplings):
        p.error('Require 0<n<1 and U>0')
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    exe = HERE/'bin/chain_Hubbard_tbc_response.x'; source = ROOT/'create_code/chain_Hubbard_tbc_response.cc'
    config = {'L': L, 'U_list': args.couplings, 'physical_densities': args.densities,
              'twist_count': args.twists, 'grid': 'periodic', 'depth': args.depth,
              'q_indices': [1], 'solver_source_sha256': digest(source), 'solver_executable_sha256': digest(exe)}
    file = output/'config.json'
    if file.exists():
        previous = json.loads(file.read_text())
        if (any(previous.get(key) != value for key, value in config.items() if key != 'physical_densities')
                or not set(config['physical_densities']).issubset(previous['physical_densities'])):
            raise RuntimeError('Existing scan configuration differs; choose a new output directory')
    if not file.exists():
        shutil.copy2(source, output/'solver_source_snapshot.cc')
        shutil.copy2(exe, output/'solver_executable_snapshot.x')
        file.write_text(json.dumps(config, indent=2)+'\n')
    (output/'active_request.json').write_text(json.dumps({'physical_densities': args.densities,
                                                       'energies_only': args.energies_only}, indent=2)+'\n')
    frozen_exe = output/'solver_executable_snapshot.x'
    assert digest(frozen_exe) == config['solver_executable_sha256']
    shutil.copy2(__file__, output/'runner_source_snapshot.py')
    for helper in ['mmhk_density_workflow.py', 'postprocess_response.py', 'extract_zero_twist_velocities.py']:
        shutil.copy2(HERE/helper, output/f'{helper[:-3]}_source_snapshot.py')
    k, weights = k_grid(L, args.twists, 'periodic')
    start = time.monotonic(); all_selections = []; stats = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for U in args.couplings:
            root = output/f'U{U:g}'; root.mkdir(exist_ok=True)
            runs = Runs(root, U, L, args.depth, args.threads, frozen_exe,
                        config['solver_source_sha256'], config['solver_executable_sha256'])
            E = energy_table(runs, k, pool)
            selections = {}
            for density in args.densities:
                selected = choose_mu(E, weights, density/2, L)
                selected['target_density_Ne_over_L'] = density
                selected['actual_density_Ne_over_L'] = 2*selected['filling']
                selections[density] = selected
                print(f'U={U:g}, target n={density:g}: mu={selected["mu"]:.9f}, achievable n={2*selected["filling"]:.9f}', flush=True)
                all_selections.append({'U': U, 'density': density, 'mu': selected['mu'],
                                       'actual_density': 2*selected['filling'], 'exact_target': selected['exact_target_attained']})
            (root/'density_selections.json').write_text(json.dumps(selections, indent=2)+'\n')
            if not args.energies_only:
                if not all(s['exact_target_attained'] for s in selections.values()):
                    raise RuntimeError('A target lies in a zero-temperature density jump: choose an explicit crossing mixture or report nearest density before responses')
                tasks = sorted({(i, N) for s in selections.values() for i, N in enumerate(s['numbers'])})
                futures = {pool.submit(ground_components, runs, i, k[i], N, E[i, N]): (i, N) for i, N in tasks}
                components = {}
                for done, future in enumerate(as_completed(futures), 1):
                    key = futures[future]; components[key] = future.result()
                    print(f'U={U:g}: actual-ground response {done}/{len(tasks)} (twist={key[0]}, N={key[1]})', flush=True)
                for density, selection in selections.items():
                    save_density(runs, density, selection, k, components)
            stats.append({'U': U, 'new_solver_calls': runs.new_calls, 'summed_solver_seconds': runs.seconds})
            (output/'run_statistics.json').write_text(json.dumps({'elapsed_seconds': time.monotonic()-start,
                                                                'workers': args.workers, 'OMP_threads_per_worker': args.threads,
                                                                'couplings': stats}, indent=2)+'\n')
            write_csv(output/'density_mu_summary.csv', all_selections)
    print('Scan output:', output, 'elapsed seconds:', time.monotonic()-start, flush=True)


if __name__ == '__main__':
    main()
