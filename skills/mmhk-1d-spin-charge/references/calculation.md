# Calculation and reuse

## Configurable helper

Use the user's virtual environment. The helper defaults to **plan only** and
does not create output directories or launch solvers until `--execute`.

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/.codex/skills/mmhk-1d-spin-charge/scripts/run_response.py \
  --L 12 --densities 0.7 0.8 --couplings 4 8 12 --twists 20 \
  --depth 300 --workers 2 --threads 2 --eta 0.2 --omega-max 5 \
  --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/L12_density_scan
```

Add `--execute` when calculation is requested. Density .7/.8 is illustrative:
either can be unattainable on a particular finite twist grid. No rounding or
silent substitution is performed. Choose smaller test parameters first for a
new mixing size; the largest minimum-|Sz| dimension is binom(L,L/2)² before
translation reduction, and full number scans still reach the largest sectors
even when the target density is low.

- `--L`: even >=4, runtime mixing size; no hard-coded 10-site assumption.
- `--densities`: physical Ne/L, currently 0<n<1; multiple targets share E0 tables.
- `--couplings`: explicit U/t values, default 0..12 inclusive in steps of one.
- `--twists 1`: a single **zero-twist** ring; otherwise even >=4 periodic points.
- `--energies-only --execute`: E0 tables and mu selections, no response Lanczos
  solves/plots. Finite-grid failures can be investigated without the response run.
- `--depth`, `--workers`, `--threads`: independent Lanczos depth and resources;
  avoid CPU oversubscription. Each solver gets OMP_NUM_THREADS=threads and MKL=1.
- `--eta`, `--omega-max`, `--points`, `--peak-cutoff`: analysis settings.
- `--analyze-only`: reads L, couplings, densities, and twists from the existing
  campaign config, ignores their CLI defaults, and never runs ED.
- `--nodes`: thermodynamic BA quadrature, default 400.
- `--project`: alternative checkout root, default the existing 1d project.

Example analysis of the already completed sweep:

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/.codex/skills/mmhk-1d-spin-charge/scripts/run_response.py \
  --analyze-only --eta 0.2 --omega-max 2 \
  --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/twist20_common_mu_U0to12_n0p8
```

This writes a new `skill_analysis_eta0p2` subdirectory, not the historical
`analysis/` files. Changed settings in the same analysis directory are rejected
to avoid overwriting a different figure/data configuration. Use a separate
copy of inputs or explicitly managed new analysis destination when needed.

## Implementation and output

The helper reuses these maintained project functions:

- `twist_averaged_velocity_scan.py`: hash-checked cached canonical Runs,
  symmetry-reduced E0 table, complete ground-component detection, pole averaging.
- `mmhk_density_workflow.py`: periodic k grid, lower-convex-envelope mu selection.
- `postprocess_response.py`: continued fractions to poles and moment checks.
- `postprocess_twist_velocity_scan.py`: Lorentzian sums and refined peak search.
- `essler_bethe_ansatz/hubbard_velocities.py`: BA reference at each actual n/U.

The helper has its own generic free-fermion calculation, avoiding the hard-coded
L=10/n=.8/weight=.2 assumptions in the historical `run_twist_U0to12.py`.
Positive-U targets outside a stable plateau stop after saving all candidate
plateaus in `U*/density_selections.json`. Ask the user before selecting a
different density or adding a number-crossing ensemble. The original scripts'
hard-coded L=10 plot titles and cached density-specific BA CSV lookups are not
used for generalized analysis.

Saved artifacts:

- `config.json`, frozen native executable/source, runner/helper source snapshots.
- Per-U `ground_energy_table.npz`, per-target `n*/metadata.json` and `COMPLETE`.
- Per-target `averaged_qmin_poles.npz/csv`, actual numbers, spin grounds, mu,
  native source moment checks, source/executable and raw-data hashes.
- `skill_analysis_eta*/`: absolute-intensity PNGs and spectra NPZs,
  `velocity_comparison.csv`, `broadening_sensitivity.csv`, selected ground states,
  metadata/config with input hashes. First/dominant peak and lowest significant
  pole are separate. Unresolved peaks are blank, not invented numbers.

For U=0, f=1/2 is used only for a full degenerate grand projector; an arbitrary
continuous f is not used to force density. Free coherent spin/charge spectra
are identical. Poles obey Omega<=4*sin(pi/L) at q_min. Zero-energy terms are
stored separately as elastic_charge/spin, not turned into a low-energy peak.

## Native build and other transfer momenta

Build only if the maintained solver is missing/incompatible:

```bash
bash /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/build_response.sh
```

That script exposes QBASIS_ROOT, QBASIS_INSTALL, MKL_ROOT, TBB_ROOT and CXX.
The known installation is `/Users/cosdis/Downloads/ED/ED_Zhentao/quantum_basis`.
On this Mac the current qbasis/MKL executable uses x86_64 libraries via Rosetta;
do not mix incompatible ARM and Intel artifacts or silently install a new stack.

The native solver accepts `--q INDEX` for a supported q_m or omits it for all
q=0..L-1. It evaluates **one canonical block**, not a common-mu campaign:

```bash
OMP_NUM_THREADS=2 MKL_NUM_THREADS=1 \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/bin/chain_Hubbard_tbc_response.x \
  --L 6 --U 8 --nup 3 --ndn 2 --phase 0 --q 2 --depth 300 \
  --observable density --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/example_L6_q2
```

For a different-q grand calculation, extend the campaign's command, metadata,
fraction glob, completeness check and output schema together; do not just
relabel q_min files. The older `mmhk_density_workflow.py` supports all-q and
runtime L, but is slower and has legacy filling/sign presentation conventions.
Its `k_blocks.theta_boundary`/plot label uses +L*k while the native phase is
-k; correct the reporting to -L*k and convert reported x to physical n before
using it as a deliverable. Its default nearest-plateau policy also needs an
explicit density acceptance check before claiming target n. Do not silently
apply the q_min helper to inter-block q or a different source vertex.

## Validation before a new large campaign

1. Small even L: source norm and first moment match saved Lanczos residues;
   certify the full ground subspace and residuals; q=0 inelastic response vanishes.
2. U=0: charge=spin(2Sz), correct weight and fixed-transfer bandwidth; verify
   a degenerate grand shell and a nonzero twist against native qbasis if extending
   free handling. A generic free total residue is not always 2/L at low density.
3. SU(2): ground energies agree across actual multiplet magnetizations. At a
   single twist S(q,w) need not be q↔-q symmetric; a conjugation-symmetric average is.
4. Check eta and Lanczos-depth convergence, and test an independent small-system
   ED/source matrix. Do not interpret the lowest Ritz pole with negligible
   weight as a stable excitation edge without residue/threshold checks.
5. Use conjugation and particle-hole E0 shortcuts only for this even,
   nearest-neighbor, zero-field Hubbard ring. For changed models scan directly.
6. Reuse canonical runs only when arguments, source/executable hashes and
   ground/source settings agree. COMPLETE alone is not a cache identity.

For half-filling, electron doping or odd L, adapt and validate the native model
and helper deliberately. In particular half-filled U>0 charge is gapped, and
electron-hole mapping for the nearest-neighbor Hubbard BA uses 2-n. Do not
extrapolate the metallic small-q charge reference through n=1.

Skill regression tests (temporary outputs only):

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/.codex/skills/mmhk-1d-spin-charge/scripts/self_test.py
```

Add `--native-smoke` to include short L=4/L=6 tests of the fast qbasis solver.
Tests cover changing L/density/twists, free elastic/projector handling, particle
number variation under common mu, and Kc agreement with independent dressed
charge equations. They do not establish convergence for every future large L.
