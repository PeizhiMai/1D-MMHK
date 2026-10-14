# Zero-twist L=10 velocity-estimate comparison at n=0.9

This analysis uses the existing physical charge/spin density spectra only.
**No new student-code, Hubbard ED, or MMHK solver runs were made.** The open
LaTeX response-definition note was not changed.

## System and estimator

- L=10, Ne=9, Nup=5, Ndn=4, t=1, U/t=4,6,7,8,12.
- Physical density n=Ne/L=0.9, equivalent to the earlier filling x=0.45.
- Exactly zero twist; equal average of the degenerate ground momenta Q=2 and 8.
- Physical longitudinal spin source n_up-n_down (2Sz); normalization affects
  intensity, not the extracted frequencies.
- q1=2*pi/L=pi/5=0.62831853.
- Operational estimate `v_hat = omega_peak(q1)/q1`, using omega(0)=0 as the
  gapless-mode reference. No finite-frequency q=0 peak is used: total N and Sz
  are conserved and their inelastic q=0 spectra vanish.

The main plot uses eta/t=0.2. Peaks are local maxima in 0<omega/t<5, with at
least 0.1% of the maximum height and 1% of its prominence. Their positions are
refined by locating a zero of the exact Lorentzian-spectrum derivative.
The **first** significant peak and **strongest** low-energy peak are reported
separately, without selecting a definition based on agreement with Bethe ansatz.

Nearly identical Lanczos Ritz poles are merged within 1e-6, preserving their
total residue and first moment. This avoids treating duplicate ghost/degenerate
poles as distinct branches. A separate pole-threshold audit checks lowest poles
at residue thresholds 1e-4, 1e-3 and 1e-2 of the channel's inelastic total weight.

## Results: velocities in units ta0/hbar

| U/t | Charge, first peak | Charge, strongest peak | Charge, BA | Spin, first/strongest peak | Spin, BA |
| --- | --- | --- | --- | --- | --- |
| 4 | 1.030738 | 1.969934 | 1.749791 | 1.086592 | 1.272087 |
| 6 | 0.779908 | 1.685791 | 1.358547 | 0.837296 | 0.985824 |
| 7 | 0.706943 | 1.569777 | 1.232506 | 0.745569 | 0.877538 |
| 8 | 0.654056 | 1.475601 | 1.138850 | 0.670547 | 0.787995 |
| 12 | 0.542024 | 1.244039 | 0.933448 | 0.475864 | 0.552389 |

The first-peak spin estimate is about 14-15% below thermodynamic Bethe ansatz.
The first charge estimate is about 41-43% below, whereas the strongest-charge
estimate is about 13-33% above. At U/t=8 the charge first and strongest peaks
are at 0.410955 and 0.927147; the first peak has 98.15% of the strongest peak's
height. Thus it is not simply a negligible peak that can be discarded.

The charge split is also visible before averaging: at U/t=8, Q=2 has its
dominant q1 charge peak at 0.402249, while Q=8 has it at 0.935392. The comparison
uses the same symmetric ground-momentum average as the original calculation.

## Why these are not controlled thermodynamic velocities

There is only one hole in this ring, and an odd-electron spin doublet. The
smallest q is not arbitrarily close to zero, and inspecting q2 and q3 does not
give a uniquely trackable linear low-energy branch. The different ground
components and broadening can change which maximum is selected.

There is an additional exact commensurability warning. At n=0.9,

```
kF = pi*n/2 = 0.45*pi,
4*kF = 1.8*pi = -0.2*pi (mod 2*pi),
q1 = 0.2*pi = |folded 4*kF|.
```

The interacting charge density has a 4kF harmonic (Essler et al., Eq. 9.47,
printed pp. 311-312). Consequently the first charge feature at this particular
q cannot automatically be identified with the uniform q->0 charge mode. This
observation is a warning about the estimator; it does not prove the identity
of each finite-ring eigenstate contributing to that peak.

The Bethe-ansatz curves are thermodynamic T=0, B=0 results at n=0.9 from the
previous integral-equation calculation. Those are actual dressed velocities,
not peak estimates. Failure of the L=10 peak/q estimator to reproduce them
does not establish absence of spin-charge separation.

For a controlled velocity determination, use smaller momentum transfers at
the same density (larger systems), finite-size scaling, or appropriately
validated finite-size energy/stiffness methods. Merely replacing the first
peak by the strongest peak is not such a determination.

## Outputs and reproduction

- `zero_twist_L10_velocities_vs_bethe_ansatz.png/svg`: main comparison. ED
  estimates are markers only; lines are the independent Bethe-ansatz curves.
- `qmin_peak_selection_audit.png`: spectra and selected peaks; q1*v_BA is a
  linear small-q guide, not an exact finite-q excitation energy.
- `velocity_comparison.csv`: frequencies, both estimators, BA and deviations.
- `broadening_sensitivity.csv`: eta=0.05,0.1,0.2,0.3, not statistical error bars.
- `pole_threshold_audit.csv`: significant lowest-pole estimates.
- `momentum_dependence.csv`: q1,q2,q3 peak checks.
- `ground_component_peaks.csv`: separate ground-momentum components.
- `metadata.json`: estimator definitions and hashes of all raw inputs.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/.mplconfig \
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/extract_zero_twist_velocities.py
```
