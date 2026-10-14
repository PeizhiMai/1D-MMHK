#!/usr/bin/env python3
"""Validate common-mu ground traces, averaging, depth and a nonzero-twist ED."""
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from twist_averaged_velocity_scan import Runs, ground_components, fractions, digest
from postprocess_response import read_fraction, poles_from_fraction
from postprocess_twist_velocity_scan import spectrum, peaks
from mmhk_density_workflow import number_minimizers

BASE = HERE/'results/twist20_common_mu_velocity_scan'
OUT = BASE/'analysis'


def component_spectrum(components, omega, eta):
    dimension = sum(c['spin_flip_multiplicity']*c['ground_dimension'] for c in components)
    result = {}
    for ch in ['charge', 'spin']:
        result[ch] = np.zeros(len(omega))
        for c in components:
            e, w = fractions(c, ch, {}, [])
            result[ch] += c['spin_flip_multiplicity']*c['ground_dimension']/dimension*spectrum(omega, e, w, eta)
    return result


def main():
    config = json.loads((BASE/'config.json').read_text())
    omega = np.linspace(0, 2, 1001); eta = .2
    report = {'passed': True, 'density_and_ground_trace_checks': [], 'average_checks': [], 'depth_checks': []}
    for U in config['U_list']:
        root = BASE/f'U{U:g}'; directory = root/'n0p8'
        m = json.loads((directory/'metadata.json').read_text())
        with np.load(root/'ground_energy_table.npz') as table:
            E, k = table['E0_by_N'], table['k']
        minimizers = number_minimizers(E, m['mu_selection']['mu'])
        assert minimizers == [[b['N']] for b in m['twists']]
        assert abs(np.mean([b['N'] for b in m['twists']])/10-.8) < 1e-12
        direct = {ch: np.zeros(len(omega)) for ch in ['charge', 'spin']}
        for block in m['twists']:
            assert abs(sum(c['probability'] for c in block['components'])-1) < 1e-12
            assert block['ground_dimension'] == sum(c['spin_flip_multiplicity']*c['ground_dimension'] for c in block['components'])
            for c in block['components']:
                source = json.loads((Path(c['directory'])/'metadata.json').read_text())
                assert source['ground_subspace_complete']
                assert abs(source['E0']-E[block['twist_index'], block['N']]) < 1e-8
                assert all(g['residual'] < 1e-8 for g in source['ground_states'])
                for ch in direct:
                    for path in Path(c['directory']).glob(f'{ch}_gsQ*_q1.dat'):
                        assert digest(path) == m['raw_fraction_sha256'][str(path.resolve())]
                        d = read_fraction(path); e, w = poles_from_fraction(d)
                        good = e > 1e-8
                        direct[ch] += c['probability']/20*d['ensemble_weight']*spectrum(omega, e[good], w[good], eta)
        with np.load(directory/'averaged_qmin_poles.npz') as poles:
            for ch in direct:
                merged = spectrum(omega, poles[f'{ch}_energies'], poles[f'{ch}_weights'], eta)
                error = float(np.max(abs(merged-direct[ch])))
                assert error < 1e-8
                report['average_checks'].append({'U': U, 'channel': ch, 'unmerged_vs_saved_average_max_error': error})
        report['density_and_ground_trace_checks'].append({'U': U, 'mu': m['mu_selection']['mu'],
                       'actual_density': .8, 'number_sectors_selected': sorted(set(b['N'] for b in m['twists'])),
                       'ground_dimensions': [b['ground_dimension'] for b in m['twists']]})
    # Complete spin traces at zero, pi/2 and pi boundary twists; not just M=0.
    root = BASE/'U8'; m = json.loads((root/'n0p8/metadata.json').read_text())
    exe = BASE/'solver_executable_snapshot.x'
    deep_runs = Runs(root/'validation/depth600', 8., 10, 600, 2, exe,
                     config['solver_source_sha256'], config['solver_executable_sha256'])
    deep_components = {}
    with np.load(root/'ground_energy_table.npz') as table:
        energies, k = table['E0_by_N'], table['k']
    for i in [0, 5, 10]:
        block = m['twists'][i]
        deep = ground_components(deep_runs, i, k[i], block['N'], energies[i, block['N']])
        deep_components[i] = deep
        shallow = [c | {'metadata': json.loads((Path(c['directory'])/'metadata.json').read_text())}
                   for c in block['components']]
        for eta in [.05, .2]:
            a = component_spectrum(shallow, omega, eta); b = component_spectrum(deep, omega, eta)
            for ch in a:
                error = float(np.max(abs(a[ch]-b[ch]))); relative = error/float(b[ch].max())
                assert relative < 1e-7
                report['depth_checks'].append({'twist_index': i, 'native_theta': float(-10*k[i]),
                                               'channel': ch, 'eta': eta, 'low_energy_max_error': error,
                                               'low_energy_relative_error': relative})
    # Independently represent a nonzero-twist H in boundary-bond gauge.
    L, Nu, Nd, U, theta = 10, 4, 4, 8., -10*k[5]
    masks = lambda N: [sum(1 << j for j in inds) for inds in itertools.combinations(range(L), N)]
    states = [u | (d << L) for u in masks(Nu) for d in masks(Nd)]
    where = {s: i for i, s in enumerate(states)}
    rows, cols, values = [], [], []
    for column, state in enumerate(states):
        rows.append(column); cols.append(column)
        values.append(U*sum(((state >> j) & 1)*((state >> (L+j)) & 1) for j in range(L)))
        for spin in [0, 1]:
            bonds = [(j+1+spin*L, j+spin*L, -1.) for j in range(L-1)]
            # Native theta=L*phase: row(last,0) has exp(+i theta).
            bonds.append((L-1+spin*L, spin*L, -np.exp(1j*theta)))
            for a, b, hopping in bonds:
                for dest, src, coefficient in [(a, b, hopping), (b, a, np.conj(hopping))]:
                    if not (state >> src) & 1 or (state >> dest) & 1:
                        continue
                    sign = (-1)**((state & ((1 << src)-1)).bit_count())
                    new = state ^ (1 << src)
                    sign *= (-1)**((new & ((1 << dest)-1)).bit_count()); new |= 1 << dest
                    rows.append(where[new]); cols.append(column); values.append(coefficient*sign)
    H = coo_matrix((values, (rows, cols)), shape=(len(states), len(states))).tocsr()
    assert (H-H.conj().T).nnz == 0
    E, V = eigsh(H, k=2, ncv=30, which='SA', tol=1e-11, maxiter=4000,
                 v0=np.random.default_rng(42).normal(size=len(states)))
    order = np.argsort(E); E, V = E[order], V[:, order]; g = V[:, 0]
    assert E[1]-E[0] > 1e-8
    np.testing.assert_allclose(E[0], energies[5, 8], atol=1e-8, rtol=0)
    raised_states = [u | (d << L) for u in masks(5) for d in masks(3)]
    raised_where = {s: i for i, s in enumerate(raised_states)}; raised = np.zeros(len(raised_states), complex)
    for column, state in enumerate(states):
        for j in range(L):
            dest, src = j, L+j
            if not (state >> src) & 1 or (state >> dest) & 1:
                continue
            sign = (-1)**((state & ((1 << src)-1)).bit_count()); new = state ^ (1 << src)
            sign *= (-1)**((new & ((1 << dest)-1)).bit_count()); new |= 1 << dest
            raised[raised_where[new]] += sign*g[column]
    S2 = float(np.vdot(raised, raised).real); assert S2 < 1e-7
    expected = component_spectrum(deep_components[5], omega, .2)
    independent = {'dimension': len(states), 'theta_native': float(theta), 'E0': float(E[0]),
                   'S_squared': S2, 'ground_residual': float(np.linalg.norm(H@g-E[0]*g)), 'spectra': {}}
    for spin, ch in [(False, 'charge'), (True, 'spin')]:
        density = np.array([sum(np.exp(2j*np.pi*j/L)*(((state >> j) & 1)
                           + (-1 if spin else 1)*((state >> (L+j)) & 1))
                           for j in range(L))/np.sqrt(L) for state in states])
        source = density*g; weight = np.vdot(source, source).real
        v = source/np.sqrt(weight); previous = np.zeros_like(v); beta = 0.; a, b = [], []
        for iteration in range(600):
            residual = H@v-beta*previous
            alpha = np.vdot(v, residual).real; residual -= alpha*v
            next_beta = np.linalg.norm(residual); a.append(alpha)
            if iteration:
                b.append(beta)
            if next_beta < 1e-12:
                break
            previous, v = v, residual/next_beta; beta = next_beta
        ep, z = eigh_tridiagonal(np.asarray(a), np.asarray(b), lapack_driver='stev')
        actual = spectrum(omega, ep-E[0], weight*z[0]**2, .2)
        error = float(np.max(abs(actual-expected[ch]))); relative = error/float(expected[ch].max())
        assert relative < 1e-6, (ch, error, relative)
        independent['spectra'][ch] = {'maximum_absolute_error': error, 'relative_to_maximum': relative}
    report['independent_nonzero_twist_ED'] = independent
    assert digest(ROOT_NOTE := HERE.parent/'note/spin_charge_response_definitions.tex') == '11860f6cf8fc410c07f595093b528fcd2737398f0fbdca023228c053d07f100f'
    report['open_note_unchanged'] = True; report['student_solver_executed'] = False
    (OUT/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS density selection, full ground traces, direct spectrum averaging, depth and independent nonzero-twist ED.', flush=True)
    print(json.dumps(independent, indent=2), flush=True)


if __name__ == '__main__':
    main()
