#!/usr/bin/env python3
"""BA velocities, small-q linear-LL curves, and exact U=0 lattice response.

Interacting curves are NOT exact finite-q Bethe-ansatz structure factors.
This helper uses the maintained zero-field BA solver in the specified project.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.special import roots_legendre

PROJECT = Path(__file__).resolve().parents[3]


def load_ba(project):
    path = Path(project)/'response/essler_bethe_ansatz'
    sys.path.insert(0, str(path))
    from hubbard_velocities import solve
    return solve


def reference(solve, U, n, nodes=400):
    """Kc from BA compressibility; check the finite-density derivative."""
    row = solve(U, n, nodes)
    if U == 0:
        return dict(row, K_charge=1., K_spin=1., K_charge_derivative_error=0.)
    h = min(2e-4, n/8, (1-n)/8)
    if h < 1e-7:
        raise ValueError('Too close to a density endpoint for this numerical Kc derivative')
    derivatives = []
    for step in [h, h/2]:
        plus = solve(U, n+step, nodes)['mu_bare_Hubbard']
        minus = solve(U, n-step, nodes)['mu_bare_Hubbard']
        derivatives.append((plus-minus)/(2*step))
    if min(derivatives) <= 0:
        raise RuntimeError('Unphysical BA compressibility')
    error = abs(derivatives[1]/derivatives[0]-1)
    if error > 2e-4:
        raise RuntimeError('BA compressibility derivative needs refinement')
    Kc = np.pi*row['v_charge']/(2*derivatives[1])
    if not .5-1e-4 <= Kc <= 1+1e-4:
        raise RuntimeError('Unexpected Kc for repulsive metallic zero-field Hubbard')
    return dict(row, K_charge=float(Kc), K_spin=1., K_charge_derivative_error=float(error))


def broaden(omega, energy, weight, eta):
    return weight*eta/(np.pi*((np.asarray(omega)-energy)**2+eta**2))


def free_poles(n, q, nodes=500):
    """Quadrature of the exact thermodynamic tight-binding U=0 continuum.

    S0=2 integral dp/(2*pi) f(p)[1-f(p+q)] delta(w-eps(p+q)+eps(p)).
    Occupation discontinuities are integrated on separate intervals.
    """
    if not 0 < n < 1 or not 0 < q <= np.pi:
        raise ValueError('Require 0<n<1 and 0<q<=pi')
    kf = np.pi*n/2
    cuts = [-kf, kf]
    for sign in [-1, 1]:
        for turn in range(-2, 3):
            value = sign*kf-q+2*np.pi*turn
            if -kf < value < kf:
                cuts.append(value)
    cuts = sorted(set(cuts)); x, w = roots_legendre(nodes)
    E, W = [], []
    for left, right in zip(cuts[:-1], cuts[1:]):
        midpoint = (left+right)/2
        dest = (midpoint+q+np.pi) % (2*np.pi)-np.pi
        if abs(dest) < kf:
            continue
        p = midpoint+(right-left)*x/2
        delta = -2*np.cos(p+q)+2*np.cos(p)
        if np.min(delta) < -1e-10:
            raise RuntimeError('Negative-energy free particle-hole excitation')
        E.extend(delta); W.extend(w*(right-left)/(2*np.pi))
    return np.asarray(E), np.asarray(W)


def continuum(omega, E, W, eta):
    result = np.zeros(len(omega))
    for start in range(0, len(E), 128):
        e, w = E[start:start+128], W[start:start+128]
        result += np.sum(w[:, None]*eta/np.pi /
                         ((np.asarray(omega)[None, :]-e[:, None])**2+eta**2), axis=0)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project', type=Path, default=PROJECT)
    p.add_argument('--density', type=float, required=True)
    p.add_argument('--couplings', type=float, nargs='+', default=list(range(13)))
    p.add_argument('--q', type=float, nargs='+', required=True, help='radians, not q indices')
    p.add_argument('--eta', type=float, default=.2)
    p.add_argument('--nodes', type=int, default=400)
    p.add_argument('--omega-max', type=float, default=5)
    p.add_argument('--points', type=int, default=5001)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if (not 0 < a.density < 1 or a.nodes < 50 or a.points < 3 or a.eta <= 0 or a.omega_max <= 0
            or not np.isfinite([a.density, a.eta, a.omega_max, *a.couplings, *a.q]).all()
            or any(U < 0 for U in a.couplings) or any(not 0 < q <= np.pi for q in a.q)):
        p.error('Invalid metallic zero-field Hubbard reference configuration')
    a.couplings = sorted(set(a.couplings)); a.q = sorted(set(a.q))
    config = {key: value for key, value in vars(a).items() if key != 'output'}
    config['project'] = str(a.project.resolve())
    source = a.project/'response/essler_bethe_ansatz/hubbard_velocities.py'
    config['BA_source_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
    config['helper_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a.output.mkdir(parents=True, exist_ok=True)
    manifest = a.output/'config.json'
    if manifest.exists() and json.loads(manifest.read_text()) != config:
        raise RuntimeError('Different reference settings: choose a fresh output directory')
    manifest.write_text(json.dumps(config, indent=2)+'\n')
    solve = load_ba(a.project)
    rows = [reference(solve, U, a.density, a.nodes) for U in a.couplings]
    with (a.output/'velocities_and_K.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    omega = np.linspace(0, a.omega_max, a.points)
    charge = np.empty((len(rows), len(a.q), len(omega))); spin = np.empty_like(charge)
    for i, row in enumerate(rows):
        for j, q in enumerate(a.q):
            charge[i, j] = broaden(omega, row['v_charge']*q, row['K_charge']*q/np.pi, a.eta)
            spin[i, j] = broaden(omega, row['v_spin']*q, q/np.pi, a.eta)
    arrays = dict(omega=omega, q=a.q, U=a.couplings, density=a.density, eta=a.eta,
                  S_charge_linear_LL=charge, S_spin_2Sz_linear_LL=spin)
    free_checks = []
    if 0 in a.couplings:
        exact = []
        for q in a.q:
            E, W = free_poles(a.density, q)
            exact.append(continuum(omega, E, W, a.eta))
            free_checks.append({'q': q, 'weight': float(W.sum()), 'first_moment': float(E@W),
                                'analytic_weight': min(q/np.pi, a.density),
                                'analytic_first_moment': float(8*np.sin(np.pi*a.density/2)*np.sin(q/2)**2/np.pi)})
            if abs(free_checks[-1]['weight']-free_checks[-1]['analytic_weight']) > 1e-10:
                raise RuntimeError('Free static-weight check failed')
            if abs(free_checks[-1]['first_moment']-free_checks[-1]['analytic_first_moment']) > 1e-10:
                raise RuntimeError('Free f-sum check failed')
        arrays['S_charge_U0_exact'] = np.asarray(exact)
        arrays['S_spin_2Sz_U0_exact'] = np.asarray(exact)
    np.savez_compressed(a.output/'reference_curves.npz', **arrays)
    (a.output/'metadata.json').write_text(json.dumps({
        'regime': 'T=B=0; metallic repulsive nearest-neighbor Hubbard; t=a0=hbar=1',
        'interacting_curves': 'BA-informed linear LL approximation, not exact finite-q response',
        'normalization': '1/sqrt(L) coherent density; spin is 2Sz; delta has unit integral',
        'Kc_method': 'pi*vc/(2*dmu_bare/dn), independent of Galilean invariance',
        'eta': 'external Lorentzian smoothing, not a predicted physical linewidth',
        'free_response_checks': free_checks,
        'warning': 'Small-q approximation; finite-q continuum, edge exponents, marginal spin logs, weak charge spin features omitted'
    }, indent=2)+'\n')
    print(a.output.resolve())


if __name__ == '__main__':
    main()
