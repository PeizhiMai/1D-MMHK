#!/usr/bin/env python3
"""L10, 50 distinct twists: inspect density plateaus, then accepted responses.

All energies from fast qbasis, all spin ground components traced. Prior runs
are copied only after matching model, source/executable hashes, phase and depth.
Student code is not executed. No approximate density is accepted implicitly.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

from twist_averaged_velocity_scan import Runs, energy_table, ground_components, save_density, digest, tag
from mmhk_density_workflow import k_grid, choose_mu
from extract_zero_twist_velocities import merge_poles, write_csv

HERE = Path(__file__).resolve().parent
SKILL = Path('/Users/cosdis/.codex/skills/mmhk-1d-spin-charge/scripts')
sys.path.insert(0, str(SKILL))
from run_response import free_case

L, NK, DEPTH = 10, 50, 300
TARGETS = [.7, .8, .9]
OUTPUT = HERE/'results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
OLD = [HERE/'results/twist20_common_mu_U0to12_n0p8', HERE/'results/twist20_common_mu_velocity_scan']


class ReusingRuns(Runs):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.reused_calls = 0

    def run(self, ik, k, Nu, Nd, ground_only=True):
        kind = 'energy_scan' if ground_only else 'qmin_responses'
        dest = self.root/kind/f'twist{ik:03d}'/f'Nup{Nu}_Ndn{Nd}'
        if not (dest/'COMPLETE').exists():
            # Grids overlap at 10 physical phases, not at equal array indices.
            old_grid, _ = k_grid(L, 20, 'periodic')
            matches = np.flatnonzero(abs(old_grid-k) < 2e-14)
            for old in OLD if len(matches) else []:
                src = old/f'U{self.U:g}'/kind/f'twist{int(matches[0]):03d}'/f'Nup{Nu}_Ndn{Nd}'
                if not (src/'COMPLETE').exists() or not (src/'lineage.json').exists():
                    continue
                lineage = json.loads((src/'lineage.json').read_text())
                m = json.loads((src/'metadata.json').read_text())
                if (lineage['source_sha256'] != self.source_hash or lineage['executable_sha256'] != self.exe_hash
                        or (m['L'], m['U'], m['t'], m['Nup'], m['Ndn'], m['depth'], m['ground_only'])
                        != (L, self.U, 1, Nu, Nd, DEPTH, ground_only)
                        or abs(m['phase_per_bond']+k) > 2e-14 or m['observable'] != 'density'
                        or not m['ground_scan_all_Q'] or (not ground_only and (m['q_indices'] != [1] or not m['ground_subspace_complete']))):
                    continue
                if dest.exists() and any(dest.iterdir()):
                    raise RuntimeError('Refusing to replace a partially written destination')
                shutil.copytree(src, dest, dirs_exist_ok=True)
                command = [str(self.exe), '--L', str(self.L), '--U', str(self.U), '--nup', str(Nu), '--ndn', str(Nd),
                           '--phase', repr(float(-k)), '--depth', str(self.depth), '--observable', 'density',
                           '--output', str(dest.resolve())]+(['--ground-only'] if ground_only else ['--q', '1'])
                (dest/'lineage.json').write_text(json.dumps({'command': command, 'source_sha256': self.source_hash,
                                                           'executable_sha256': self.exe_hash}, indent=2)+'\n')
                (dest/'reuse_receipt.json').write_text(json.dumps({'original_directory': str(src.resolve()),
                    'original_lineage': lineage, 'copied_file_sha256': {p.name: digest(p) for p in src.iterdir() if p.is_file()},
                    'phase_difference': m['phase_per_bond']+k}, indent=2)+'\n')
                with self.lock:
                    self.reused_calls += 1
                break
        return super().run(ik, k, Nu, Nd, ground_only)


def density_report(U, E, weights):
    selections, rows = {}, []
    for n in TARGETS:
        selected = choose_mu(E, weights, n/2, L)
        selected.update(target_density_Ne_over_L=n, actual_density_Ne_over_L=2*selected['filling'])
        below = [s for s in selected['plateaus'] if 2*s['filling'] <= n+1e-10]
        above = [s for s in selected['plateaus'] if 2*s['filling'] >= n-1e-10]
        lo = max(below, key=lambda s: s['filling']); hi = min(above, key=lambda s: s['filling'])
        rows.append(dict(U=U, target_n=n, nearest_stable_n=2*selected['filling'],
                         difference=2*selected['filling']-n, exact_stable_target=selected['exact_target_attained'],
                         density_below=2*lo['filling'], density_above=2*hi['filling'], mu=selected['mu'],
                         mu_left=selected['mu_left'], mu_right=selected['mu_right'],
                         selected_N_histogram=json.dumps({str(N): selected['numbers'].count(N) for N in sorted(set(selected['numbers']))})))
        selections[str(n)] = selected
    return selections, rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=['energies', 'responses'], required=True)
    p.add_argument('--workers', type=int, default=6)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output', type=Path, default=OUTPUT)
    p.add_argument('--couplings', type=int, nargs='+', default=list(range(13)))
    a = p.parse_args(); out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    source = HERE.parent/'create_code/chain_Hubbard_tbc_response.cc'; exe = HERE/'bin/chain_Hubbard_tbc_response.x'
    config = dict(L=L, U_list=list(range(13)), physical_densities=TARGETS, twist_count=NK, grid='periodic',
                  depth=DEPTH, q_indices=[1], solver_source_sha256=digest(source), solver_executable_sha256=digest(exe),
                  runner_sha256=digest(__file__), free_helper_sha256=digest(SKILL/'run_response.py'),
                  ground_helper_sha256=digest(HERE/'twist_averaged_velocity_scan.py'))
    file = out/'config.json'
    if file.exists() and json.loads(file.read_text()) != config:
        raise RuntimeError('Existing campaign physics/source settings differ')
    if not file.exists():
        file.write_text(json.dumps(config, indent=2)+'\n')
        for origin, name in [(source, 'solver_source_snapshot.cc'), (exe, 'solver_executable_snapshot.x'),
                             (Path(__file__), 'runner_source_snapshot.py'), (SKILL/'run_response.py', 'free_helper_source_snapshot.py')]:
            shutil.copy2(origin, out/name)
    frozen = out/'solver_executable_snapshot.x'
    assert digest(frozen) == config['solver_executable_sha256']
    k, weights = k_grid(L, NK, 'periodic'); reports, statistics = [], []
    approval = out/'density_policy.json'
    policy = json.loads(approval.read_text()) if approval.exists() else {'nearest_stable_tolerance': 0.}
    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for U in a.couplings:
            root = out/f'U{U}'; root.mkdir(exist_ok=True)
            if U == 0:
                for n in TARGETS:
                    selected = free_case(root, L, k, n, config, merge_poles, write_csv, tag, energies_only=a.stage == 'energies')
                    print('FREE DENSITY', n, selected, flush=True)
                continue
            runs = ReusingRuns(root, U, L, DEPTH, a.threads, frozen, config['solver_source_sha256'], config['solver_executable_sha256'])
            if (root/'ground_energy_table.npz').exists():
                with np.load(root/'ground_energy_table.npz') as d:
                    np.testing.assert_allclose(d['k'], k, rtol=0, atol=1e-14); E = d['E0_by_N'].copy()
            else:
                E = energy_table(runs, k, pool)
            selections, rows = density_report(U, E, weights); reports.extend(rows)
            (root/'density_selections.json').write_text(json.dumps(selections, indent=2)+'\n')
            write_csv(out/f'density_plateaus_{a.stage}.csv', reports)
            for r in rows:
                print(f'DENSITY U={U}: target {r["target_n"]:.3f}, stable {r["nearest_stable_n"]:.3f}, '
                      f'bracket [{r["density_below"]:.3f},{r["density_above"]:.3f}], exact={r["exact_stable_target"]}', flush=True)
            if a.stage == 'responses':
                accepted = {float(n): s for n, s in selections.items()
                            if s['exact_target_attained'] or abs(2*s['filling']-float(n)) <= policy.get('nearest_stable_tolerance', 0.)+1e-12}
                accepted = {n: s for n, s in accepted.items() if not (root/f'n{tag(n)}'/'COMPLETE').exists()}
                tasks = sorted({(i, N) for s in accepted.values() for i, N in enumerate(s['numbers'])})
                futures = {pool.submit(ground_components, runs, i, k[i], N, E[i, N]): (i, N) for i, N in tasks}
                components = {}
                for count, f in enumerate(as_completed(futures), 1):
                    components[futures[f]] = f.result()
                    if count % 20 == 0 or count == len(futures):
                        print(f'U={U}: response grounds {count}/{len(futures)}', flush=True)
                for n, s in accepted.items():
                    actual = float(2*s['filling'])
                    save_density(runs, actual, s, k, components)
                    produced = root/f'n{tag(actual)}'; desired = root/f'n{tag(n)}'
                    if produced != desired:
                        if desired.exists():
                            raise RuntimeError('Approximate density output name already exists')
                        produced.rename(desired)
                        m = json.loads((desired/'metadata.json').read_text())
                        m.update(target_density_Ne_over_L=n, density_approximation_policy=policy,
                                 actual_density_Ne_over_L=actual, native_filling_Ne_over_2L=actual/2)
                        (desired/'metadata.json').write_text(json.dumps(m, indent=2)+'\n')
            statistics.append(dict(U=U, new_calls=runs.new_calls, reused_calls=runs.reused_calls, solver_seconds=runs.seconds))
            (out/f'run_statistics_{a.stage}.json').write_text(json.dumps(dict(elapsed_seconds=time.monotonic()-start,
                                                        workers=a.workers, threads=a.threads, per_U=statistics), indent=2)+'\n')
    print('STAGE COMPLETE', a.stage, out, flush=True)


if __name__ == '__main__':
    main()
