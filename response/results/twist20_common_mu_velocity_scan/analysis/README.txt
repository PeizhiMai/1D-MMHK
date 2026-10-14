20-twist, common-mu Hubbard density responses and n(mu)

Scope: L=10, t=1, U/t=4,6,7,8,12. Physical density n=Ne/L=0.8.
The request for n=0.9 response calculations was withdrawn. The already
completed U=4,6 n=0.9 files are retained unchanged. No additional n=0.9
response runs were made. Its density staircase uses shared energy tables only.

Grid: twenty distinct uniform boundary twists, spacing pi/10, including zero
and one antiperiodic endpoint. Native per-bond phase=theta/L. The native
boundary row(last,0) has exp(+i theta); row(0,last) has exp(-i theta).

State: one common mu(U,n) minimizes E0(theta,N)-mu*N at each twist. Neither
N nor Sz is imposed on the response ground trace. A minimum-|Sz| energy
scan contains every spin multiplet; all magnetizations and degenerate ground
vectors at the minimum are then included with their projector weights.

Result: at n=0.8 this U range has a stable plateau with N=8 on ALL twenty
twists. This is a grand-canonical selection result, not an N=8 constraint.
Most twists have singlet grounds; near zero twist some have triplet grounds,
and the full triplet is traced. Spin is not artificially fixed at Sz=0.

Physical observables: rho_c/s(q)=sum_j exp(i*q*j)(n_up +/- n_dn)/sqrt(L).
Only q1=2*pi/10=pi/5 is measured for this velocity comparison. Spectra are
positive-frequency structure factors, not finite-eta retarded losses.
Neutral Em-Eg frequencies need no chemical-potential shift.

Average all ground components and all twenty twist spectra BEFORE extracting
peaks. Estimator v_hat=omega_peak(q1)/q1. eta/t=0.2, 0<omega/t<5, local
maximum height>=0.1% and prominence>=1%; first/strongest peaks are retained
separately. The strongest peaks agree with thermodynamic BA within about
1.1-3.2% (charge) and 1.4-2.1% (spin) at this eta. These are finite-q,
broadening-dependent estimates, not q->0 fits or finite-size extrapolations.
For example, at smaller eta the dominant U=6 charge peak can switch.

Density plot: n(mu) is a genuine staircase of lower-convex-envelope ground
energies, not an interpolated density root. At U/t=7, mu_c/t=1.430946883264:
n(mu_c-0)=0.89, n(mu_c+0)=0.91. The twists theta=+/-0.9*pi each switch N=8
->10. At the crossing, number sectors are degenerate and the density depends
on the choice/mixture. There is no stable mu INTERVAL with n=0.9. This
finite-cluster/twist-grid step is not evidence of thermodynamic phase separation.

Validation: density selection over every N, complete spin ground traces,
source norms/first moments, unmerged-vs-saved average, and depth=300 vs 600
at U=8 for theta=0,pi/2,pi all pass. Independent 44,100-state boundary-gauge
ED at nonzero twist theta=pi/2 reproduces both low-energy spectra with
relative errors below 3e-13. Ground energy shortcuts are conjugation and
particle-hole symmetries, directly checked at each U. Responses are computed
directly at all twenty twists; they are not generated from mirrored spectra.
No student module was imported or executed. The open LaTeX note is unchanged.

Main artifacts:
  twist20_density_vs_mu.png/svg
  density_mu_plateaus.csv; density_jumps_across_n0p9.json
  twist20_velocities_vs_bethe_ansatz.png/svg; velocity_comparison.csv
  twist20_qmin_peak_audit.png; broadening_sensitivity.csv
  twist20_selected_particle_numbers.png; selected_ground_states.csv
  validation.json; density_mu_metadata.json; metadata.json
Raw run metadata, fractions and frozen solver/source are in the parent tree.

Reproduce postprocessing (no solver runs):
  cd /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    ~/.venvs/myenv/bin/python response/postprocess_twist_velocity_scan.py --densities 0.8
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    ~/.venvs/myenv/bin/python response/plot_twist_density_mu.py
