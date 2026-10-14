---
name: mmhk-1d-spin-charge
description: Calculate physical spin/charge dynamical density responses for 1D MMHK or Hubbard TBC rings at configurable density, mixing size, and twists; extract finite-q velocity estimates and thermodynamic Bethe-ansatz references. Use for these neutral two-particle responses, not single-particle Green functions.
---

# 1D MMHK spin/charge response

Use the validated fast qbasis Hubbard-ring solver in the repository's
`response/` directory. The bundled scripts default to this checkout's root;
another checkout can be supplied with `--project`. This skill depends on that
project's maintained solver/helpers. Use the user's Python environment (the
original machine uses `~/.venvs/myenv/bin/python`). Historical absolute paths
in references identify the original calculation, not a required install path.

## Establish the calculation contract

Keep these independent parameters explicit:

- **L**: number of mixed momenta / dual-space orbitals, not the twist count.
  Current native solver supports even L >= 4. Changing L changes the Hilbert
  space and q_min = 2*pi/L; it does not require recompiling this runtime-L solver.
- **n = Ne/L**: physical total density; half-filling is n=1. Legacy input
  x=Ne/(2L)=n/2 uses x=0.5 for half-filling. Never silently interchange them.
- **Twists / reduced-BZ points**: one for a single zero-twist ring, or a grid
  for the reduced-BZ response. Increasing twists does NOT reduce q_min.
- Couplings U/t, supported transfer q_m=2*pi*m/L, eta/t, frequency range,
  Lanczos depth, and CPU resources. Defaults reproduce the n=.8 study, but
  inherit explicit user choices rather than always using L=10 or 20 twists.

The helper supports q_min, t=a0=hbar=1, T=B=0, U>=0, and 0<n<1. For other q,
half-filling, electron doping, finite field/temperature, or a modified mixing
geometry, use [calculation.md](references/calculation.md) rather than pretending
the helper implements them. Do not add interactions to the Hubbard Hamiltonian
and still call the BA reference an exact reference for the modified model.

## Physical invariants

Read [physics.md](references/physics.md) when selecting operators or ensembles.

- Use the **coherent** physical charge density and longitudinal **2Sz** density,
  with 1/sqrt(L) normalization. Student Fourier-mode occupations or incoherent
  orbital-pair intensity sums are different observables.
- Use each twist's actual global ground subspace. Scan all Q and include all
  magnetizations of degenerate ground multiplets; do not force Sz=0 or a singlet.
- For common-mu calculations minimize E0(theta,N)-mu*N across number sectors.
  Do not impose N=round(nL) separately at every twist. Average the full ground
  trace and twists **before** extracting peaks.
- A finite T=0 twist grid has density steps. If no stable mu interval attains
  the target, report adjacent plateaus and ask about an explicit crossing
  ensemble, nearest density, or changing the grid. Never silently fabricate
  the requested density. At U=0 the helper may use the full grand ground
  projector at a Fermi-level crossing when it actually gives the target.
- Neutral frequencies are Em-Eg: **no chemical-potential shift**.
- Distinguish positive Lehmann S(q,w) from the retarded commutator response,
  especially at finite eta. Separate elastic weight from inelastic response.
- Never run/import the student's slow 10-MMHK solver. If testing that specific
  student code is necessary, honor the user's limit of 6-MMHK.

## Run or reuse a campaign

Read [calculation.md](references/calculation.md) for runnable examples, output
formats, and solver validation. `scripts/run_response.py` provides a plan by
default; `--execute` runs it, `--analyze-only` uses existing outputs. It exposes
L, densities, couplings, twists, depth, workers, threads, eta, and omega range.
Only execute calculations within the user's requested scope.

Preserve existing results. Check physics settings and source/executable hashes
before reusing them; use a fresh output root for changed parameters. Estimate
the canonical dimension binom(L,Nup)*binom(L,Ndn) before larger-L runs. Test a
small case and scale CPU concurrency without workers*threads oversubscription.

The maintained U=0..12, L=10, n=.8, 20-twist results are at
`/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/twist20_common_mu_U0to12_n0p8`.
Those are validated examples, not immutable settings for later tasks.

## Velocity and frequency-dependent reference

Read [bethe-reference.md](references/bethe-reference.md) for definitions,
limitations, primary sources, and commands. `scripts/bethe_reference.py`
calculates configurable BA velocities and Kc, and writes **BA-informed linear
Luttinger-liquid** S_c/s(q,w) curves. It also supplies the exact thermodynamic
tight-binding U=0 particle-hole response at finite q. Neither curve is the
exact interacting finite-q Hubbard/MMHK response.

Report first resolved peak, dominant peak, lowest significant pole, residues,
and eta sensitivity separately. omega_peak/q_min is a **finite-q estimator**,
not a proven sound velocity. Do not select peaks just to improve BA agreement.
If actual thermodynamic line shapes are requested, explain the need for
form factors or a converged dynamical method (e.g. dynamical/time-dependent
DMRG). BA dressed energies/velocities alone do not give intensities or widths.

Deliver absolute-intensity spectra, pole tables, selected number/spin grounds,
mu/density evidence, velocity comparisons, checks, and source provenance. A
normalized-shape view is supplementary, not a replacement for absolute weights.
Do not edit the open LaTeX note unless the user asks to revise it.
