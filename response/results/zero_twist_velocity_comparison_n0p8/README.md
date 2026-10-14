# Zero-twist L=10 Hubbard response at physical density n=0.8

Physical density is n=Ne/L=0.8: Ne=8 and the spin-balanced sector is
Nup=Ndn=4. In the older native solver convention the filling is x=Ne/(2L)=0.4.
This is a canonical physical-density response of one zero-twist Hubbard ring,
not the student's reduced-BZ-summed Fourier-mode-occupation response.

The fast qbasis solver was run at U/t=4,6,7,8,12, phase=0, t=1, all ten q values,
depth=300. The student solver was never executed. No LaTeX source was changed.

## Ground-state spin is an important finite-size detail

For all five couplings, the ground state is Q=0 and a spin triplet, S=1:
the unique ground energies in M=Sz=0 and M=1 agree to better than 3e-14,
while the lowest energy in M=2 is strictly higher. Since any S>=2 multiplet
would also appear in M=2, this identifies the ground spin as S=1. An independent
unreduced 44,100-state calculation at U=8 confirms <S^2>=2.

The fixed-Nup=Ndn=4 response is the M=0 triplet member, **not a singlet**.
It is a valid canonical fixed-spin-sector result. A rotationally invariant
zero-field ground ensemble instead gives

```
S_iso(q,w) = [S_M0(q,w) + 2*S_M1(q,w)]/3.
```

M=-1 has the same longitudinal spectrum as M=1 by exact global spin flip.
The spin-scalar charge response is identical in all triplet components; this
is checked numerically. The longitudinal spin response is not identical, so
the state/ensemble must be specified when discussing finite-size peak estimates.
The q=0 inelastic response is zero in either prescription.

## Operational estimates and comparison

```
q1 = 2*pi/L = pi/5,
v_hat = omega_dominant_peak(q1)/q1, eta/t=0.2.
```

The peak is the strongest local maximum in 0<omega/t<5, with 1% prominence
and at least 0.1% of maximum height. Peak positions are refined from the exact
Lorentzian derivative. No definition was selected to force agreement with BA.
Velocities below have units ta0/hbar. BA is thermodynamic, T=0, B=0, n=0.8.

| U/t | Charge estimate | Charge BA | Spin, fixed M=0 | Spin, equal triplet | Spin BA |
| --- | --- | --- | --- | --- | --- |
| 4 | 1.957516 | 2.063589 | 1.501924 | 1.102650 | 1.250097 |
| 6 | 1.795203 | 1.892745 | 1.206098 | 0.874079 | 0.993767 |
| 7 | 1.723944 | 1.813383 | 1.085378 | 0.787373 | 0.893249 |
| 8 | 1.663168 | 1.744857 | 0.982860 | 0.715097 | 0.808356 |
| 12 | 1.501525 | 1.561418 | 0.704105 | 0.520316 | 0.577305 |

The dominant charge estimate is 3.8-5.2% below BA. The fixed-M0 spin estimate
is 20.1-22.0% above BA; the equal-triplet estimate is 9.9-12.0% below BA.
At U=8 the dominant charge peak is omega/t=1.045000, the M0 spin peak is
0.617549, and the equal-triplet spin peak is 0.449309.

## Qualifications: not a thermodynamic velocity determination

- These are finite-q peak/q estimates, not fitted q->0 slopes. q1=pi/5 is
  still coarse, and there has been no finite-size extrapolation.
- The dominant charge peak is more unambiguous than in the Ne=9 data, but
  weaker lower charge poles remain. At U=8 one occurs at omega/t=0.434620.
  Narrower broadening resolves additional first peaks at some U values;
  "first resolved peak at eta=0.2" is not "lowest exact charge excitation".
- Equal-triplet spin maxima move with broadening because several poles
  contribute. The eta=0.02,0.05,0.1,0.2,0.3 audit is saved; this is sensitivity,
  not statistical uncertainty. The charge dominant maxima are much more stable.
- At n=0.8 the folded thermodynamic 4kF magnitude is 0.4*pi=q2, not q1.
  Thus the exact q1=folded-4kF coincidence of n=0.9 is absent, but this alone
  does not establish a linear branch or determine the identity of each pole.
- Distinct dominant spin/charge features with the expected U trends provide
  a better qualitative benchmark here. They do not by themselves establish
  thermodynamic spin-charge velocities or remove the need for finite-size scaling.

## Verification and provenance

- Every qbasis run scans all Q sectors and checks source norms/first moments.
  Zeroth/first moment errors are below 1e-13. Ground-vector residuals are below 1e-12.
- The independent U=8 unreduced ED reproduces M0 charge/spin spectra with
  maximum absolute errors 7.4e-13 and 1.2e-10, respectively.
- U=8 q1 depth=300 versus 600 is checked for both M0 and M1 at eta=0.02,0.05,0.2.
  Dominant peak positions and the omega/t<=2 window agree to numerical precision.
  At narrow eta, weak high-energy spin tails above omega/t=2 differ more and
  are not claimed converged to that same precision; the errors are reported.
- BA density, mesh, kernel, dressed-energy derivative and limiting checks pass.
- Raw-vs-merged pole spectra/peaks, inversion, old n=0.9 regression, raw/source
  hashes and the unchanged open LaTeX note are checked in validation.json.
- Raw runs are in `../zero_twist_benchmark_n0p8` and
  `../zero_twist_benchmark_n0p8_Sz1`; immutable solver/source snapshots are saved.

## Outputs

- `zero_twist_L10_n0p8_ensemble_velocities_vs_bethe_ansatz.png/svg`: main
  comparison, including fixed-M0 spin markers. Charge is ensemble independent.
- `triplet_velocity_comparison.csv`: comparison table, frequencies and deviations.
- `U8_ground_ensemble_and_poles.png`: spectra before/after the triplet average.
- `triplet_broadening_sensitivity.csv`: first/dominant peak choices at five eta.
- `triplet_merged_poles.csv`: retained poles and their fractional residues.
- `velocity_comparison.csv` and `zero_twist_L10_velocities_vs_bethe_ansatz.png/svg`:
  original fixed-M0-only extraction, first/dominant peaks at eta=0.2.
- `pole_threshold_audit.csv`, `momentum_dependence.csv`, `ground_component_peaks.csv`:
  fixed-M0 diagnostics, including q1,q2,q3.
- `independent_U8_validation.json/npz` and `validation.json`: numerical checks.
- `metadata.json`, `triplet_metadata.json`, source snapshots: analysis provenance.

## Reproduce postprocessing (does not rerun a Hamiltonian solver)

```bash
cd /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR="$PWD/response/.mplconfig" \
~/.venvs/myenv/bin/python response/extract_zero_twist_velocities.py --density 0.8

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR="$PWD/response/.mplconfig" \
~/.venvs/myenv/bin/python response/compare_n0p8_spin_ensemble.py
```

The generic density runner is `response/run_density_velocity_benchmark.py`;
use a fresh output path if its provenance differs from an existing run. The
BA solver is `response/essler_bethe_ansatz/hubbard_velocities.py --density 0.8`.
