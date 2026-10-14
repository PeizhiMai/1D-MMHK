#!/usr/bin/env python3
"""Fast, error-bounded Lorentzian grids; peak positions refined exactly.

Deposit residues linearly on an energy grid (preserving M0/M1), then convolve
with the full Lorentzian kernel by FFT. No poles/tails are dropped. Grid error
is <= sum(W)*h**2/(4*pi*eta**3) before output interpolation. Peak frequencies
are roots of the exact original-pole derivative, not histogram bin locations.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.signal import fftconvolve, find_peaks

from postprocess_twist_velocity_scan import spectrum as direct_spectrum

GRID = .0005


def spectrum(omega, E, W, eta, h=GRID):
    omega = np.asarray(omega); E = np.asarray(E); W = np.asarray(W)
    if len(E) < 100 or len(omega) < 33 or np.min(omega) < 0:
        return direct_spectrum(omega, E, W, eta)
    if eta <= 0 or h <= 0 or np.any(E < 0) or np.any(W < 0):
        raise ValueError('Invalid Lorentzian grid inputs')
    M = int(np.ceil(E.max()/h))+1
    J = int(np.ceil(omega.max()/h))+1
    u = E/h; i = np.floor(u).astype(int); f = u-i
    masses = np.bincount(i, weights=W*(1-f), minlength=M+1)
    masses += np.bincount(i+1, weights=W*f, minlength=M+1)
    shifts = np.arange(-M, J+1)*h
    kernel = eta/(np.pi*(shifts*shifts+eta*eta))
    conv = fftconvolve(masses, kernel, mode='full')[M:M+J+1]
    return np.interp(omega/h, np.arange(J+1), conv)


def peaks(E, W, eta, cutoff=5.):
    omega = np.linspace(0, cutoff, round(cutoff/GRID)+1)
    values = spectrum(omega, E, W, eta)
    # Generous candidate thresholds, followed by exact height/prominence checks.
    indices, props = find_peaks(values, height=.0009*values.max(), prominence=.009*values.max())
    def slope(x):
        z = x-E
        return float(np.sum(-2*W*eta*z/(np.pi*(z*z+eta*eta)**2)))
    def exact(x):
        return float(np.sum(W*eta/(np.pi*((x-E)**2+eta**2))))
    candidates = []
    for j, i in enumerate(indices):
        lo, hi = omega[max(0, i-2)], omega[min(len(omega)-1, i+2)]
        if not slope(lo) > 0 > slope(hi):
            # Do not manufacture a maximum when a histogram artifact is suspected.
            from postprocess_twist_velocity_scan import peaks as direct_peaks
            return direct_peaks(E, W, eta, cutoff)
        x = brentq(slope, lo, hi, xtol=1e-12)
        bases = []
        for key in ['left_bases', 'right_bases']:
            b = int(props[key][j]); xbase = omega[b]
            if 0 < b < len(omega)-1:
                left, right = omega[b-1], omega[b+1]
                if slope(left) < 0 < slope(right):
                    xbase = brentq(slope, left, right, xtol=1e-12)
            bases.append(exact(xbase))
        candidates.append(dict(omega=float(x), height=exact(x), prominence=exact(x)-max(bases)))
    if not candidates:
        raise RuntimeError('No resolved positive-frequency local maximum')
    maximum = max(p['height'] for p in candidates)
    result = [{'omega': p['omega'], 'height': p['height']} for p in candidates
              if p['height'] >= .001*maximum and p['prominence'] >= .01*maximum]
    if not result:
        raise RuntimeError('No resolved positive-frequency local maximum')
    return result
