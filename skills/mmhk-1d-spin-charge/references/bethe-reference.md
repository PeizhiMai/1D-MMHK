# Bethe ansatz: velocities versus frequency-dependent response

## Exact velocity reference already implemented

Maintained source:
`/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/essler_bethe_ansatz/hubbard_velocities.py`.

It solves the T=B=0 thermodynamic repulsive nearest-neighbor Hubbard integral
equations for physical n=Ne/L, u=U/(4t). Eqs. 5.102–5.105 of Essler et al. fix
root densities and dressed energies; Eqs. 6.44–6.47 give dressed-energy /
dressed-momentum slopes at the Fermi boundaries, not bare rapidity slopes.
The zero-field spin boundary is infinite and requires a stable asymptotic ratio.

Local primary source:
`/Users/cosdis/Downloads/books/Essler, Frahm, Gohmann, Klumper, Korepin - The One-Dimensional Hubbard Model.pdf`.
Printed pp.166–167 and 187–188 (PDF pages184–185 and205–206).

Checks: quadrature refinement, integral-equation residuals, exact U=0
v_c=v_s=2*sin(pi*n/2), half-filled spin Bessel limit, and strong-U asymptotes.
For new n/U, compute the reference instead of interpolating an unrelated
density-specific CSV. At n=1,U>0 the charge sector is gapped; a vanishing
boundary slope is not the velocity of a gapless charge sound peak.

## Why velocities do not determine the full spectrum

At fixed q, S(q,w) needs both excitation energies and density matrix elements.
BA dressed energies/dispersion boundaries do not supply the Lehmann residues.
Do not call a Lorentzian centered at v*q an exact interacting BA line shape.
Edges, peak maxima, weights, nonlinear continuum widths, and high-energy
Hubbard features are different quantities. Finite-L MMHK also need not have
the same thermodynamic line shape as the Hubbard BA reference.

Reference levels, in increasing scope:

1. **U=0 exact full finite-q response:** a tight-binding particle-hole integral
   is straightforward, and charge=spin(2Sz). The helper computes its broadened
   thermodynamic continuum with interval-split Gauss quadrature. This is distinct
   from the finite-grid/ring free response in run_response.py.
2. **BA-informed linear LL response:** valid for the smooth, long-wavelength
   metallic sector at T=B=0. With this skill's unit-area delta convention,

   S_c(q,w) ≈ Kc*abs(q)/pi * delta(w-v_c*abs(q))
   S_2Sz(q,w) ≈ abs(q)/pi * delta(w-v_s*abs(q)).

   This includes a leading absolute intensity, not just the peak position.
   SU(2) gives Ks=1 in the spinful-fermion convention. The helper obtains
   Kc=pi*v_c/[2*d(mu_bare)/dn] from BA compressibility and checks the density
   finite-difference step. A dressed-charge calculation is an alternative.
   **Do not use Kc=vF/vc on a lattice**; that is a Galilean-invariant relation.
   eta broadens the delta externally and is not a predicted physical lifetime.
   Spin marginal logarithms and nonlinear spin-charge features are omitted.
3. **Nonlinear LL / mobile-impurity edge references:** can add continuum
   boundaries, effective masses, and threshold exponents when justified by a
   model-specific BA calculation. They still need amplitudes/matrix elements
   for a full absolute spectrum. Do not transplant a Galilean quantum-wire
   formula directly into the lattice Hubbard model.
4. **Full interacting S(q,w):** requires evaluated form factors and enough
   excited states/continua, or a separately converged dynamical method such as
   dynamical/time-dependent DMRG. An exact half-filled two-spinon contribution
   exists but is neither both full responses nor the n=.8 doped reference.
   Inspect any candidate implementation's model, density, polarization,
   truncation, normalization, sum rules and accessibility before promising it.

At L=10, n=.8, q_min/kF=.5, so q_min is not asymptotically small. A broadened
linear-LL overlay is a leading benchmark, not a controlled full line-shape
prediction at that q. Increasing L changes this ratio; increasing twists alone
does not. A full finite-q reference should be compared with identical source,
density, q, eta and normalization, with size/boundary differences declared.

## Helper command

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/.codex/skills/mmhk-1d-spin-charge/scripts/bethe_reference.py \
  --density 0.8 --couplings 0 1 2 3 4 5 6 7 8 9 10 11 12 \
  --q 0.6283185307179586 --eta 0.2 --nodes 400 --omega-max 5 \
  --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/BA_LL_reference_n0p8_qpi5
```

Outputs: `velocities_and_K.csv`, `reference_curves.npz`, `metadata.json`,
hash-checked `config.json`. Arrays have axes (U,q,omega) for linear LL and
(q,omega) for exact U=0. The latter are included only if U=0 is requested.
Every interacting curve is explicitly labelled **linear_LL**, not exact_BA.
The exact free curve verifies static weight min(q/pi,n) and the f-sum
8*sin(pi*n/2)*sin(q/2)²/pi. Refine free quadrature if using very small eta.

## Primary literature, checked 2026-10-13

- [Pereira and Sela, Spin-charge coupling in quantum wires at zero magnetic field](https://arxiv.org/abs/0911.1391),
  especially [Eqs.36–45 and Sec.IV](https://arxiv.org/html/0911.1391).
  Their time-Fourier S has an extra factor 2*pi relative to our unit-area
  Lehmann convention. Supports the linear-LL delta reference and weaker
  spin-derived charge response beyond it; does not provide a turnkey doped
  Hubbard full-response solver. Their Galilean-specific coefficient relations
  are not automatically applicable to the Hubbard lattice.
- [Essler and Korepin, Formfactors in the half-filled Hubbard model](https://arxiv.org/abs/cond-mat/9808018),
  PRB59,1734 (1999). Two-spinon spin contribution at half-filling; not the
  complete doped spin and charge structure factors.
- [Pereira et al., Charge dynamics in half-filled Hubbard chains with finite on-site interaction](https://arxiv.org/abs/1111.2009),
  PRB85,165132 (2012). Uses dynamical numerical spectra together with exact
  excitation information/edge theory; useful example of why energies alone
  do not fix intensity. Half-filling is not n=.8.

Recheck current research/available implementations when a future task requests
a full interacting reference. Do not turn this bounded source check into a
claim that no such implementation or later result can exist.
