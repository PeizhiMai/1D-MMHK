Completed L=10 1D-MMHK physical density-response campaign

39 cases: U/t=0..12 inclusive; target n=Ne/L=0.7,0.8,0.9; 50 distinct periodic twists.
q1=2*pi/10=pi/5; t=a0=hbar=1; T=B=0. One common chemical potential per U/target.
Each twist uses its actual number sector and complete spin/translation ground projector.
No student solver imported or executed; no open LaTeX note edited.

Actual density:
- n=0.7 and 0.8 attained at every U.
- target n=0.9 attained except U7:0.896, U9:0.902, U10:0.898.
- Nearest stable plateaus within ±0.004 were explicitly approved by the user.
- At U0,n0.8, the full degenerate free grand ground projector gives exactly 0.8;
  there is no stable mu plateau at that exact density (adjacent plateaus 0.796/0.804).
  This is the same natural full free-ground trace convention used previously.
- For positive U no number-crossing mixture was manufactured.
- Positive-U n0.8 happens to select N8 at all twists; it was not imposed.

Observable and plotting:
rho_c/s(q)=sum_j exp(i*q*j)(n_up +/- n_down)/sqrt(L). Spin is longitudinal 2Sz.
Spin intensity can be divided by 4 for Sz; peak energies are unchanged.
Positive-frequency inelastic Lehmann spectrum with Lorentzian eta/t=0.2.
Elastic terms are excluded; neutral Em-Eg frequencies have no chemical-potential shift.
Main response panels use identical absolute-intensity axes for all U/densities.
Full-range log panels and saved arrays include the entire stored pole range plus tails.

Velocities:
vhat=dominant peak frequency/q1 after ground/twist averaging, not averaging per-twist peaks.
First peak, lowest significant pole (>1e-6 total residue), and eta sensitivity also saved.
At eta0.2 the first resolved and dominant peaks coincide in all 39 cases.
BA curves are thermodynamic velocities at each actual mean density, not at the same finite-cluster mu.
Nominal-density BA values are also included in velocity_comparison.csv.
These vhat values are finite-q/eta-dependent estimators, not controlled q->0 sound velocities.
More twists improve averaging but do not change qmin.

Key comparison over U4..12:
- n0.7: charge peak/q about 3.1–3.6% below BA; spin 10.9–14.7% below BA.
- n0.8: charge 1.0–3.2% below BA; spin 1.2–2.1% below BA.
- near n0.9: charge -0.7% to +3.2% versus BA; spin 8.2–12.1% below BA.
Thus n0.8 gives the closest spin/charge velocity benchmark for this finite L, q and eta.
This observation is not a proof of thermodynamic velocity convergence.

Validation/provenance:
- 39 completed cases; unique physical twists, normalized weights and exact number selection checked.
- Source norms/first moments, ground residuals, equal charge/spin f-sum moments checked.
- All raw fraction hashes were rechecked by the final postprocessor.
- U0 charge and spin pole arrays identical; fixed-q free bandwidth bound satisfied.
- BA integer-U references tested at 400/600 nodes for n0.7,0.8,0.9; additional
  n0.7 kernel/dressed-derivative/limiting checks passed.
- Selected 300/600-Lanczos-depth anchors tested for N7,N9,N10, including nonzero twist.
  Primary eta0.2 spectrum changes <1.2e-9 of peak height and dominant peak shifts <1.2e-15.
  At eta0.05, maximum tested spectrum discrepancy is 2.15e-5 of peak height;
  the original strict 1ppm small-eta test therefore does not pass, and this is
  explicitly retained in validation/depth_convergence.json. No global depth proof claimed.
- Plot grids use error-bounded FFT convolution on h=0.0005 with linear residue
  deposition preserving M0/M1; all poles/tails retained. Peak frequencies use
  exact original-pole derivative roots. Direct-sum checks satisfy the analytic
  grid-error bound; tested peak differences <1e-13.
- Prior raw canonical runs are reused only at matching phases, sectors, depth,
  operator and native source/executable hashes. Original data remain unchanged.

Files:
response_panel_n0p7/0p8/0p9.png, .svg: 13 U panels each.
spin_charge_velocity_vs_U.png, .svg: both velocities together, one panel per density.
actual_density_vs_U.png, .svg: mean densities actually attained.
full_spectra_n*_page1/2/3.png, .svg: low-energy and full-frequency views.
velocity_comparison.csv, broadening_sensitivity.csv, actual_density.csv,
selected_ground_states.csv, n0p8_twist_count_comparison.csv.
Per-case low/full spectrum NPZs; original averaged poles in ../U*/n*.
../density_diagnostics: exact mu plateaus, jump metadata and U7 density-vs-mu figure.
../bethe_reference: nominal BA velocities and numerical validation.
../validation: Lanczos-depth and FFT/direct-sum comparisons.

20→50 twist sampling check at n0.8: maximum change in dominant peak/q is
0.4894% for charge and 0.2891% for spin over U0..12.
Reuse audit: 1163 canonical runs reused; 8161 original files rehashed unchanged.
