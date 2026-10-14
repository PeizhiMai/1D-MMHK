Direct single-particle A(k,w) for L10, U/t0..12, targets physical n=Ne/L=.7,.8,.9.
Exactly101 uniform physical momenta k=j*pi/100, j=0..100. No interpolation from the251-point run.
Eta/t=.1; stored omega/t=-6..16 with step.01; Lanczos depth600.
Spectral colormaps display omega/t=-6..10; stored NPZ data are unchanged.
Per-spin A=-Im(Gup+Gdn)/(2*pi) has unit total spectral weight.
Addition poles Em(N+1)-E0-mu; removal poles E0-Em(N-1)-mu.
Reference density and common mu remain those of the original50-twist response campaign.
The101 k grid selects20 folded phases, with11 inversion representatives. It changes
spectral sampling only. Never retune mu using this coarser set of phases.
case_summary.csv reports both reference_n_50_twists and sampled_n_20_folded_phases,
which can differ. Figure n_ref labels refer to the50-twist density constraint.
At off-grid phases, solve actual grand ground states at the same mu. Do not impose N or Sz.
The full number-sector scan is rigorously pruned using the hopping norm bound
||H(phi+d)-H(phi)||<=4Lt|sin(d/2)|. Solve candidate sectors over ALL translation Q,
check the lower bounds of excluded N, then trace every ground spin multiplet/projector.
Shared phases reuse only hash-matched native600-step charged sources.
Physical momentum K-phi=2*pi*q/L+k_RBZ equals the external q=0 twist measurement
after a discrete large-gauge transformation. No averaging of distinct physical k.
U0 uses the exact free retarded Green function, independent of free occupation degeneracy.
Rainbow colour range0..1.5, with an upper saturation arrow. NPZ values are unmodified.
Legacy half-filled notebooks plot -2ImGup=2*pi*A_per_spin under spin symmetry.
spectral_panel_n*:13 U panels per density. U*/n*: individual PNG/SVG, Akw.npz and metadata.
Selected600/900 depth checks and independent L4 full-Fock ED support the implementation.
Selected depth tests do not establish a global all-case numerical error bound.
Spin/charge response plots and the four-slide PPT retain50 twists and eta.1.
The student solver and open LaTeX note remain unchanged.
