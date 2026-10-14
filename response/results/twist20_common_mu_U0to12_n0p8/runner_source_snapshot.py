#!/usr/bin/env python3
"""Extend n=0.8, twenty-twist physical density responses to U/t=0,...,12.

Positive-U cases use the fast qbasis solver and a common chemical potential.
U=0 is evaluated exactly with the full grand-canonical Fermi-level projector,
not a selected fixed-number or fixed-spin Slater determinant. Existing cases
are copied with provenance, leaving the original campaign untouched.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import itertools
import json
import os
from pathlib import Path
import shutil
import time

import numpy as np

from twist_averaged_velocity_scan import Runs, energy_table, ground_components, save_density, digest
from mmhk_density_workflow import k_grid, choose_mu
from extract_zero_twist_velocities import merge_poles, write_csv

HERE = Path(__file__).resolve().parent
OLD = HERE/'results/twist20_common_mu_velocity_scan'
OUT = HERE/'results/twist20_common_mu_U0to12_n0p8'
L, NK, DEPTH = 10, 20, 300
DENSITY = .8


def free_case(root, k, config):
    """Exact finite-grid free-fermion E0(N), ground projector and rho spectra."""
    K = 2*np.pi*np.arange(L)/L
    epsilon = -2*np.cos(K[None, :]+k[:, None])  # phase=-k, epsilon(K)=-2*cos(K-phase)
    levels = np.sort(np.repeat(epsilon, 2, axis=1), axis=1)
    E0 = np.column_stack([np.zeros(NK), np.cumsum(levels, axis=1)])
    np.savez_compressed(root/'ground_energy_table.npz', k=k, weights=np.ones(NK)/NK, E0_by_N=E0)
    global_levels = np.sort(levels.ravel())
    Ntotal = round(DENSITY*L*NK)
    mu = float(.5*(global_levels[Ntotal-1]+global_levels[Ntotal]))
    tol = 1e-10
    f = np.where(epsilon < mu-tol, 1., np.where(epsilon > mu+tol, 0., .5))
    meanN = 2*f.sum(axis=1)
    assert abs(meanN.mean()/L-DENSITY) < 1e-12
    blocks, gap, residue = [], [], []
    for i in range(NK):
        # Explicitly count all degenerate free ground configurations.
        filled = [(s, j) for s in [0, 1] for j in range(L) if f[i, j] == 1]
        shell = [(s, j) for s in [0, 1] for j in range(L) if f[i, j] == .5]
        groups = {}
        number_histogram = {}
        for bits in itertools.product([0, 1], repeat=len(shell)):
            occupied = filled+[mode for bit, mode in zip(bits, shell) if bit]
            Nu = sum(s == 0 for s, _ in occupied); Nd = len(occupied)-Nu
            # qbasis labels translation by Q=-sum(K indices) for the stated
            # c_K=sum_j exp(+iKj)c_j convention (rho(+q) then targets Q+q).
            Q = (-sum(j for _, j in occupied)) % L
            groups.setdefault((Nu, Nd), []).append(Q)
            number_histogram[Nu+Nd] = number_histogram.get(Nu+Nd, 0)+1
        dimension = 2**len(shell)
        components = []
        for (Nu, Nd), momenta in sorted(groups.items()):
            if Nu < Nd:
                continue
            mult = 1 if Nu == Nd else 2
            components.append({'N': Nu+Nd, 'Nup': Nu, 'Ndn': Nd,
                               'ground_dimension': len(momenta), 'spin_flip_multiplicity': mult,
                               'probability': mult*len(momenta)/dimension,
                               'ground_states': [{'Q': Q, 'id': j} for j, Q in enumerate(momenta)],
                               'method': 'exact Slater ground-projector trace'})
        assert sum(c['ground_dimension']*c['spin_flip_multiplicity'] for c in components) == dimension
        assert abs(sum(N*count/dimension for N, count in number_histogram.items())-meanN[i]) < 1e-12
        blocks.append({'twist_index': i, 'k': float(k[i]), 'phase_per_bond': float(-k[i]),
                       'theta_boundary': float(-L*k[i]), 'N': float(meanN[i]),
                       'N_support': {str(N): count/dimension for N, count in number_histogram.items()},
                       'number_variance': float(sum(count*(N-meanN[i])**2/dimension for N, count in number_histogram.items())),
                       'ground_dimension': dimension, 'components': components})
        # rho(+q)=sum_K c_K^dagger c_(K+q)/sqrt(L), with c_K=sum_j exp(+iKj)c_j/sqrt(L).
        for dest in range(L):
            src = (dest+1) % L
            probability = f[i, src]*(1-f[i, dest])
            if probability > 0:
                delta = float(epsilon[i, dest]-epsilon[i, src])
                assert delta > 1e-8
                gap.append(delta); residue.append(2*probability/L/NK)
    energies, weights = merge_poles(np.asarray(gap), np.asarray(residue))
    assert abs(weights.sum()-.2) < 1e-12
    assert energies.max() <= 4*np.sin(np.pi/L)+1e-12
    output = root/'n0p8'; output.mkdir(exist_ok=True)
    np.savez_compressed(output/'averaged_qmin_poles.npz', k=k, theta_boundary=-L*k,
                        N_per_twist=meanN, twist_weights=np.ones(NK)/NK, q_min=2*np.pi/L,
                        U=0., mu=mu, actual_density=DENSITY,
                        charge_energies=energies, charge_weights=weights,
                        spin_energies=energies, spin_weights=weights)
    rows = [{'channel': ch, 'omega': float(e), 'weight': float(w)}
            for ch in ['charge', 'spin'] for e, w in zip(energies, weights)]
    write_csv(output/'averaged_qmin_poles.csv', rows)
    selection = {'mu': mu, 'target_density_Ne_over_L': DENSITY, 'actual_density_Ne_over_L': DENSITY,
                 'filling': .4, 'exact_target_attained': True,
                 'policy': 'complete zero-temperature grand ground projector; half occupied exactly degenerate Fermi modes',
                 'finite_grid_density_below_crossing': float(2*(epsilon < mu-tol).sum()/(L*NK)),
                 'finite_grid_density_above_crossing': float(2*(epsilon <= mu+tol).sum()/(L*NK))}
    metadata = {'L': L, 'U': 0., 'target_density_Ne_over_L': DENSITY, 'actual_density_Ne_over_L': DENSITY,
                'native_filling_Ne_over_2L': .4, 'twist_count': NK, 'mu_selection': selection,
                'q_min': 2*np.pi/L, 'twists': blocks,
                'response_method': 'exact free-fermion particle-hole sum for all twenty twists',
                'ensemble': 'full grand-canonical ground projector; all number and spin states at the Fermi level',
                'observable': 'coherent physical charge and longitudinal 2Sz densities',
                'frequency_origin': 'neutral excitation energy, no chemical-potential shift',
                'source_moment_checks': [], 'raw_fraction_sha256': {},
                'solver_source_sha256': config['solver_source_sha256'],
                'solver_executable_sha256': config['solver_executable_sha256'],
                'student_solver_executed': False}
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    (root/'density_selections.json').write_text(json.dumps({'0.8': selection}, indent=2)+'\n')
    np.savez_compressed(output/'exact_free_reference.npz', k=k, epsilon=epsilon, f=f, mu=mu,
                        E0_by_N=E0, mean_N=meanN, energies=energies, weights=weights)
    (output/'COMPLETE').write_text('success\n')
    print(f'COMPLETE U=0 exact free response: mu={mu:.12f}, n={meanN.mean()/L:.12f}, charge=spin', flush=True)


def reuse_case(U, config):
    source = OLD/f'U{U:g}'; dest = OUT/f'U{U:g}'; original = source/'n0p8'
    if not (original/'COMPLETE').exists():
        return False
    m = json.loads((original/'metadata.json').read_text())
    if (m['L'] != L or m['twist_count'] != NK or m['actual_density_Ne_over_L'] != DENSITY
            or m['solver_source_sha256'] != config['solver_source_sha256']
            or m['solver_executable_sha256'] != config['solver_executable_sha256']):
        raise RuntimeError('Old case is not physics/provenance compatible')
    old_config = json.loads((OLD/'config.json').read_text())
    if old_config['depth'] != DEPTH or old_config['grid'] != 'periodic' or old_config['q_indices'] != [1]:
        raise RuntimeError('Old solver settings are incompatible')
    dest.mkdir(exist_ok=True); (dest/'n0p8').mkdir(exist_ok=True)
    receipt = {}
    paths = ['ground_energy_table.npz', 'energy_symmetry_checks.json',
             'n0p8/averaged_qmin_poles.npz', 'n0p8/averaged_qmin_poles.csv', 'n0p8/metadata.json', 'n0p8/COMPLETE']
    for name in paths:
        old, new = source/name, dest/name
        if new.exists() and digest(new) != digest(old):
            raise RuntimeError('Refusing to overwrite a different reused output')
        shutil.copy2(old, new); receipt[str(old.resolve())] = digest(old)
    (dest/'reuse_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(f'REUSED U={U:g}: completed common-mu ground traces, no new solver runs', flush=True)
    return True


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = HERE.parent/'create_code/chain_Hubbard_tbc_response.cc'
    exe = HERE/'bin/chain_Hubbard_tbc_response.x'
    config = {'L': L, 'U_list': list(range(13)), 'physical_densities': [DENSITY],
              'twist_count': NK, 'grid': 'periodic', 'depth': DEPTH, 'q_indices': [1],
              'solver_source_sha256': digest(source), 'solver_executable_sha256': digest(exe)}
    file = OUT/'config.json'
    if file.exists() and json.loads(file.read_text()) != config:
        raise RuntimeError('Existing U sweep has different configuration')
    if not file.exists():
        shutil.copy2(source, OUT/'solver_source_snapshot.cc')
        shutil.copy2(exe, OUT/'solver_executable_snapshot.x')
        file.write_text(json.dumps(config, indent=2)+'\n')
    frozen = OUT/'solver_executable_snapshot.x'; assert digest(frozen) == config['solver_executable_sha256']
    (OUT/'active_request.json').write_text(json.dumps({'physical_densities': [DENSITY], 'energies_only': False}, indent=2)+'\n')
    shutil.copy2(__file__, OUT/'runner_source_snapshot.py')
    k, weights = k_grid(L, NK, 'periodic')
    root = OUT/'U0'; root.mkdir(exist_ok=True)
    free_case(root, k, config)
    summaries, statistics = [], []
    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:
        for U in range(1, 13):
            if reuse_case(U, config):
                m = json.loads((OUT/f'U{U}'/'n0p8/metadata.json').read_text())
                summaries.append({'U': U, 'mu': m['mu_selection']['mu'], 'actual_density': .8, 'reused': True})
                continue
            root = OUT/f'U{U}'; root.mkdir(exist_ok=True)
            runs = Runs(root, U, L, DEPTH, 2, frozen, config['solver_source_sha256'], config['solver_executable_sha256'])
            E = energy_table(runs, k, pool)
            selection = choose_mu(E, weights, .4, L)
            selection.update(target_density_Ne_over_L=DENSITY, actual_density_Ne_over_L=2*selection['filling'])
            (root/'density_selections.json').write_text(json.dumps({'0.8': selection}, indent=2)+'\n')
            print(f'U={U}: mu={selection["mu"]:.12f}, selected n={2*selection["filling"]:.12f}', flush=True)
            if not selection['exact_target_attained']:
                raise RuntimeError('Positive-U target lies in a density jump; do not silently substitute a different density')
            futures = {pool.submit(ground_components, runs, i, k[i], N, E[i, N]): (i, N)
                       for i, N in enumerate(selection['numbers'])}
            components = {}
            for count, future in enumerate(as_completed(futures), 1):
                key = futures[future]; components[key] = future.result()
                print(f'U={U}: actual-ground response {count}/{NK}', flush=True)
            save_density(runs, DENSITY, selection, k, components)
            summaries.append({'U': U, 'mu': selection['mu'], 'actual_density': .8, 'reused': False})
            statistics.append({'U': U, 'new_solver_calls': runs.new_calls, 'summed_solver_seconds': runs.seconds})
            write_csv(OUT/'density_mu_summary.csv', summaries)
            (OUT/'run_statistics.json').write_text(json.dumps({'elapsed_seconds': time.monotonic()-start,
                                                               'new_cases': statistics}, indent=2)+'\n')
    summaries.insert(0, {'U': 0, 'mu': float(np.load(OUT/'U0/n0p8/averaged_qmin_poles.npz')['mu']),
                         'actual_density': .8, 'reused': False})
    write_csv(OUT/'density_mu_summary.csv', summaries)
    print('ALL 13 COUPLINGS COMPLETE:', OUT, flush=True)


if __name__ == '__main__':
    main()
