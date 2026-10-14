L=10 physical Hubbard density responses: n=Ne/L=0.8, 20 twists, U/t=0,...,12.

Same physical operators and q=pi/5 as the previous calculation. Responses
are positive-frequency structure factors, with rho_c/s=sum exp(iqj)(n_up +/- n_dn)/sqrt(L).
Spin is 2Sz. eta/t=0.2 for the primary plots; intensities are absolute.
One common chemical potential is chosen at every U. Full actual ground
multiplets are traced; neither fixed Sz=0 nor a per-twist N=8 constraint
is imposed. All positive-U selections happen to pick N=8 on every twist.

Completed U=4,6,7,8,12 datasets were reused without modifying the originals.
New qbasis cases are U=1,2,3,5,9,10,11. The student's solver was never run.

U=0 is evaluated by the exact free particle-hole formula and checked against
qbasis at theta=0, pi/2, pi. mu/t=-0.618033988749895. At zero twist four
Fermi-level spin orbitals are degenerate. The full 16-state grand ground
projector has N=6..10, with probabilities (1,4,6,4,1)/16, mean N=8 and
variance 1. Other twists have N=8. Thus the average density is exactly 0.8,
without silently replacing it by the neighboring plateaus 0.79 or 0.81.
Charge and longitudinal 2Sz spectra are identical at U=0. The independent
qbasis comparison errors at the three twists are below 6e-14.

13 coupling values and 20 twist contributions at every coupling are complete.
Source norm/first-moment checks, density selection, complete ground traces,
file hashes and reuse invariance pass; see validation.json. The open LaTeX
note was not changed. Source/executable and analysis snapshots are preserved.

Outputs:
  n0p8_spectra_and_peak_audit/U0to12_spin_charge_low_energy.png/svg
  n0p8_spectra_and_peak_audit/all_U_spin_charge_spectra_U0to4.png/svg/npz
  n0p8_spectra_and_peak_audit/all_U_spin_charge_spectra_U5to9.png/svg/npz
  n0p8_spectra_and_peak_audit/all_U_spin_charge_spectra_U10to12.png/svg/npz
  twist20_velocities_vs_bethe_ansatz.png/svg; velocity_comparison.csv
  broadening_sensitivity.csv; selected_ground_states.csv
  n0p8_spectra_and_peak_audit/peak_vs_lowest_pole.csv and related audits

Peak/q values remain finite-q and broadening-dependent proxies, not
extrapolated sound velocities. First and dominant peaks are retained
separately; no assumption from the earlier five-U sweep is imposed on
new couplings. Weak-group classification is skipped for the free case.

Reproduce postprocessing without rerunning solvers:
cd /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 ~/.venvs/myenv/bin/python response/postprocess_twist_velocity_scan.py \
 --input response/results/twist20_common_mu_U0to12_n0p8 --densities 0.8
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 ~/.venvs/myenv/bin/python response/plot_twist_response_spectra.py \
 --input response/results/twist20_common_mu_U0to12_n0p8
