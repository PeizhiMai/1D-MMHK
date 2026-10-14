#!/usr/bin/env python3
"""Independent 44,100-state U=8, Ne=8 check and ground-state spin diagnostic."""
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import postprocess_response as post


def masks(L, N):
    return [sum(1 << j for j in inds) for inds in itertools.combinations(range(L), N)]


def main():
    start = time.monotonic()
    L, U, nup, ndn = 10, 8., 4, 4
    states = [u | (d << L) for u in masks(L, nup) for d in masks(L, ndn)]
    where = {s: i for i, s in enumerate(states)}
    row, col, data = [], [], []
    for k, s in enumerate(states):
        row.append(k); col.append(k)
        data.append(U*sum(((s >> j) & 1)*((s >> (j+L)) & 1) for j in range(L)))
        for sigma in [0, 1]:
            for j in range(L):
                a, b = j+L*sigma, (j+1) % L+L*sigma
                for dest, src in [(a, b), (b, a)]:
                    if not (s >> src) & 1 or (s >> dest) & 1:
                        continue
                    sign = (-1)**((s & ((1 << src)-1)).bit_count())
                    v = s ^ (1 << src)
                    sign *= (-1)**((v & ((1 << dest)-1)).bit_count())
                    v |= 1 << dest
                    row.append(where[v]); col.append(k); data.append(-sign)
    H = coo_matrix((data, (row, col)), shape=(len(states), len(states))).tocsr()
    assert (H-H.T).nnz == 0
    print('Independent matrix', len(states), 'states, nnz', H.nnz, flush=True)
    energies, V = eigsh(H, k=4, ncv=40, which='SA', tol=1e-11, maxiter=4000,
                        v0=np.random.default_rng(42).normal(size=len(states)))
    order = np.argsort(energies); energies, V = energies[order], V[:, order]
    ground = np.flatnonzero(energies-energies[0] < 1e-8)
    assert len(ground) == 1, energies
    g = V[:, 0]
    residual = float(np.linalg.norm(H @ g-energies[0]*g))
    assert residual < 1e-8
    # S+ from independently represented fermions: for M=0, ||S+g||^2=S(S+1).
    raised_states = [u | (d << L) for u in masks(L, nup+1) for d in masks(L, ndn-1)]
    raised_where = {s: i for i, s in enumerate(raised_states)}
    raised = np.zeros(len(raised_states))
    for k, s in enumerate(states):
        for j in range(L):
            dest, src = j, L+j
            if not (s >> src) & 1 or (s >> dest) & 1:
                continue
            sign = (-1)**((s & ((1 << src)-1)).bit_count())
            v = s ^ (1 << src)
            sign *= (-1)**((v & ((1 << dest)-1)).bit_count())
            v |= 1 << dest
            raised[raised_where[v]] += sign*g[k]
    s2 = float(raised @ raised)
    S = float((np.sqrt(1+4*s2)-1)/2)
    print('Independent E0', energies[0], 'S(S+1)', s2, 'S', S, flush=True)
    omega = np.linspace(0, 5, 5001); eta = .2
    run = HERE/'results/zero_twist_benchmark_n0p8/L10_U8_n0p4_phi0'
    actual = post.reconstruct(run, omega, eta)
    np.testing.assert_allclose(energies[0], actual['metadata']['E0'], atol=1e-9, rtol=0)
    report = {'L': L, 'Ne': nup+ndn, 'Nup': nup, 'Ndn': ndn, 'U': U,
              'dimension': len(states), 'nnz': int(H.nnz), 'lowest_eigenvalues': energies.tolist(),
              'ground_degeneracy_in_M0': len(ground), 'ground_residual': residual,
              'ground_S_squared': s2, 'ground_total_spin': S, 'channels': {}}
    spectra = {}
    for spin, channel in [(False, 'charge'), (True, 'spin')]:
        density = np.array([sum(np.exp(2j*np.pi*j/L)*(((s >> j) & 1)
                             + (-1 if spin else 1)*((s >> (j+L)) & 1))
                             for j in range(L))/np.sqrt(L) for s in states])
        source = density*g; weight = np.vdot(source, source).real
        v = source/np.sqrt(weight); previous = np.zeros_like(v); beta = 0.; a, b = [], []
        for it in range(600):
            r = H @ v-beta*previous
            alpha = np.vdot(v, r).real; r -= alpha*v
            next_beta = np.linalg.norm(r); a.append(alpha)
            if it:
                b.append(beta)
            if next_beta < 1e-12:
                break
            previous, v = v, r/next_beta; beta = next_beta
        ep, z = eigh_tridiagonal(np.asarray(a), np.asarray(b), lapack_driver='stev')
        spectrum = post.lorentzian(omega, ep-energies[0], weight*z[0]**2, eta)
        spectra[channel] = spectrum
        difference = float(np.max(abs(spectrum-actual[channel][1])))
        np.testing.assert_allclose(spectrum, actual[channel][1], atol=1e-7, rtol=1e-6)
        report['channels'][channel] = {'maximum_absolute_spectrum_error': difference,
                                      'relative_to_maximum': difference/float(spectrum.max()),
                                      'first_peak': post.first_peak(omega, spectrum)}
        print('PASS independent', channel, report['channels'][channel], flush=True)
    report.update(passed=True, elapsed_seconds=time.monotonic()-start)
    out = HERE/'results/zero_twist_velocity_comparison_n0p8'
    out.mkdir(exist_ok=True)
    (out/'independent_U8_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez_compressed(out/'independent_U8_spectra.npz', omega=omega, **spectra, eta=eta)


if __name__ == '__main__':
    main()
