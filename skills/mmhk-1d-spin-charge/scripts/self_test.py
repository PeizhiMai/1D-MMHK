#!/usr/bin/env python3
"""Behavioral regression tests; optional tiny fast-qbasis integration tests.

Writes only inside a temporary directory. Does not run student code.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

import bethe_reference as ba
import run_response as campaign


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, default=campaign.PROJECT)
    p.add_argument('--native-smoke', action='store_true')
    a = p.parse_args(); sys.path.insert(0, str(a.project/'response'))
    from twist_averaged_velocity_scan import tag
    from mmhk_density_workflow import k_grid, choose_mu
    from extract_zero_twist_velocities import merge_poles, write_csv
    solve = ba.load_ba(a.project)
    from hubbard_velocities import kernel_R
    for n in [.3, .6, .8, .9]:
        for U in [0., 1., 8., 1000.]:
            row = ba.reference(solve, U, n, 200)
            if U == 0:
                assert row['K_charge'] == 1 and row['v_charge'] == row['v_spin']
                continue
            _, state = solve(U, n, 200, return_state=True)
            k, w, Q, u = (state[key] for key in ['k', 'w', 'Q', 'u'])
            R = kernel_R(np.sin(k)[None, :]-np.sin(k)[:, None], u)
            xi = np.linalg.solve(np.eye(len(k))-R*(w*np.cos(k))[None, :], np.ones(len(k)))
            xiQ = 1+(kernel_R(np.sin(k)-np.sin(Q), u)*w*np.cos(k))@xi
            assert abs(row['K_charge']-xiQ**2/2) < 2e-6
            if U == 1000:
                assert abs(row['K_charge']-.5) < .003
    for n in [.1, .5, .8, .95]:
        for q in [.04, .4, 2.5, np.pi]:
            E, W = ba.free_poles(n, q)
            assert np.all(E >= 0)
            assert abs(W.sum()-min(n, q/np.pi)) < 1e-12
            assert abs(E@W-8*np.sin(np.pi*n/2)*np.sin(q/2)**2/np.pi) < 1e-11
    # A T=0 density jump must not be turned into a continuously adjustable plateau.
    selection = choose_mu(np.array([[0., 2., 1., 6., 10.]]), np.ones(1), .25, 2)
    assert not selection['exact_target_attained']
    scripts = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix='mmhk_skill_selftest_') as temporary:
        root = Path(temporary)
        dry = root/'dry'
        cmd = [sys.executable, str(scripts/'run_response.py'), '--L', '12', '--densities', '.6', '.8',
               '--twists', '40', '--couplings', '0', '4', '8', '--output', str(dry)]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        assert not dry.exists()
        plan = json.loads(result.stdout.split('\nPlan only.')[0])
        assert plan['L'] == 12 and plan['twist_count'] == 40 and abs(plan['q_min']-np.pi/6) < 1e-12
        config = {'solver_source_sha256': 'analytic-test', 'solver_executable_sha256': 'analytic-test'}
        for L, n, nk in [(4, .25, 4), (4, .5, 4), (6, .5, 4), (12, .5, 4), (6, 1/3, 1), (10, .8, 20)]:
            case = root/f'L{L}_n{n:g}_twist{nk}'; case.mkdir()
            k, _ = k_grid(L, nk, 'periodic')
            campaign.free_case(case, L, k, n, config, merge_poles, write_csv, tag)
            d = np.load(case/f'n{tag(n)}'/'averaged_qmin_poles.npz')
            np.testing.assert_array_equal(d['charge_weights'], d['spin_weights'])
            np.testing.assert_array_equal(d['charge_energies'], d['spin_energies'])
            assert abs(float(d['actual_density'])-n) < 1e-10
            assert not len(d['charge_energies']) or d['charge_energies'].max() <= 4*np.sin(np.pi/L)+1e-10
            if L == 4 and n == .5:
                assert float(d['elastic_charge']) > 0
            if L == 10:
                m = json.loads((case/'n0p8/metadata.json').read_text())
                z = next(b for b in m['twists'] if abs(b['theta_boundary']) < 1e-12)
                assert z['ground_dimension'] == 16
                assert set(z['N_support']) == set(map(str, range(6, 11)))
        if a.native_smoke:
            for L, densities in [(4, ['.5']), (6, ['.5', '0.6666666666666666'])]:
                output = root/f'native_L{L}'
                cmd = [sys.executable, str(scripts/'run_response.py'), '--project', str(a.project),
                       '--L', str(L), '--densities', *densities, '--couplings', '4', '--twists', '4',
                       '--depth', '80', '--workers', '2', '--threads', '1', '--nodes', '200',
                       '--points', '1001', '--omega-max', '3', '--execute', '--output', str(output)]
                subprocess.run(cmd, check=True)
                for density in map(float, densities):
                    m = json.loads((output/'U4'/f'n{tag(density)}'/'metadata.json').read_text())
                    assert m['L'] == L and abs(m['actual_density_Ne_over_L']-density) < 1e-10
                    assert max(c['first_moment_error'] for c in m['source_moment_checks']) < 1e-8
                    for b in m['twists']:
                        assert abs(sum(c['probability'] for c in b['components'])-1) < 1e-12
                    if L == 6 and density == .5:
                        assert len({b['N'] for b in m['twists']}) > 1
    print('PASS: configuration, free projector/continuum, BA Kc, density-step guard'
          + (', and L4/L6 native qbasis smoke tests' if a.native_smoke else ''))


if __name__ == '__main__':
    main()
