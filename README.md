# 1D-MMHK: spin/charge responses and single-particle spectra

qbasis/Lanczos calculations for the one-dimensional repulsive Hubbard ring
with twisted boundary conditions, used as the dual-space block of 1D MMHK.
The repository includes solver sources, calculation/postprocessing scripts,
raw numerical data, validation records, Bethe-ansatz velocity references,
figures, and the four-slide spin/charge presentation.

## Current results

The main completed campaign is in
[`response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9`](response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9).

- **Mixing size:** 10 orbitals; **couplings:** U/t = 0, 1, ..., 12.
- **Physical target densities:** n = Ne/L = 0.7, 0.8, 0.9. Half-filling is n=1.
  Some historical scripts call x=Ne/(2L) the filling, so x=0.5 is half-filling.
- **Neutral responses:** 50 reduced-BZ twists, q_min = 2*pi/10, one common
  chemical potential for each (U, target density). Each twist uses its actual
  grand ground subspace, including degenerate spin/translation sectors;
  neither Sz=0 nor a fixed particle number is imposed at every twist.
- The approved nearest stable density tolerance is 0.004. Near n=0.9,
  the actual 50-twist averages are 0.896 at U/t=7, 0.902 at 9, and 0.898 at 10;
  the other targets in this campaign are attained.
- **Current broadening:** eta/t=0.1. The earlier eta=0.2 analyses are retained.
- **BA reference:** directly solved thermodynamic Hubbard velocities at
  Delta(U/t)=0.2, not an interpolation of the integer-U values.
- **Single-particle A(k,omega):** exactly 101 directly calculated momenta
  k_j=j*pi/100, j=0,...,100; eta/t=0.1; Lanczos depth 600. No momentum
  interpolation. The common chemical potentials remain those of the original
  50-twist density campaign; the 101 momenta use 20 folded phases and are not
  a new 50-twist average. The sampled density can differ from the reference
  density and is recorded separately.
- Spectral colormaps display omega/t from **-6 to 10**. Stored spectral arrays
  retain the full **-6 to 16** range with spacing 0.01.

### Figures, presentation, and data

Inside the main campaign directory:

| Directory | Contents |
| --- | --- |
| `analysis_ba_step0p2_eta0p1/` | Spin/charge velocity comparison; 13-U response panels for each density; full-range panels; peak, density, and BA CSV tables; spectrum NPZ files |
| `analysis_ba_step0p2_eta0p1/presentation/` | `spin_charge_10MMHK_50twists_eta0p1_v2.pptx`: four slides with editable native charts, covering velocities and neutral responses |
| `single_particle_eta0p1_k101_depth600/` | Three 13-U spectral panels; 39 individual PNG/SVG maps and `U*/n*/Akw.npz` arrays; summary and validation records |
| `U*/` | Ground-energy scans, chemical-potential selections, full ground ensembles, and neutral Lanczos sources |
| `density_diagnostics/` | n(mu) plateaus and target-density attainability |
| `validation/` | Neutral-response depth and other checks |

The PPT contains **neutral responses and velocity comparisons**, not the
single-particle spectral maps. The maps are delivered separately.
`single_particle_eta0p1_depth600/canonical/` contains charged sources reused
by the final 101-point campaign and must be preserved. Other pilot and historical
results are retained for provenance; they are **not** substitutes for the final
101-point output. `COMPLETE` refers to the particular directory/run, not to
global validation of every historical result.

## Definitions and limitations

The Hamiltonian is

```text
H(phi) = -t sum_(j,sigma) [exp(i phi) c†_(j,sigma)c_(j+1,sigma) + h.c.]
         + U sum_j n_(j,up)n_(j,down).
```

`phi` is the per-bond phase; the boundary twist is L*phi.
The coherent physical density operators are

```text
rho_c/s(q) = sum_j exp(i q j) [n_(j,up) +/- n_(j,down)] / sqrt(L).
```

Spin intensity uses **2Sz**; divide by four for Sz normalization. Elastic
ground-subspace weight is separated from the inelastic positive Lehmann
structure factor. Neutral frequencies are Em-Eg, with **no mu shift**.
This positive-frequency spectrum is not the full retarded commutator
susceptibility. It is also not the student's Fourier-mode occupation response
or an incoherent sum of orbital-pair intensities.

Single-particle spectra are per spin:
`A=-Im(G_up+G_down)/(2*pi)`, with total frequency-integrated weight one.
Particle/hole frequencies use the selected common mu, not an automatic U/2
shift away from half-filling.

`omega_peak(q_min)/q_min` is a **finite-q, broadening-dependent estimator**,
not a demonstrated q->0 sound velocity. Dominant peaks, first resolved peaks,
and low-weight poles are distinct diagnostics. The BA velocities are exact
thermodynamic reference quantities for the Hubbard model, not predictions
for all finite-cluster peak locations. Dressed BA velocities alone do not
supply the interacting finite-q frequency-dependent response or its weights.

Reference: F. H. L. Essler, H. Frahm, F. Göhmann, A. Klümper and V. E. Korepin,
*The One-Dimensional Hubbard Model*, Cambridge University Press (2005),
Eqs. (5.102)-(5.105), (6.44)-(6.47). See
[`response/essler_bethe_ansatz/README.md`](response/essler_bethe_ansatz/README.md)
for the equations, conventions, and numerical checks. The book is not redistributed.

## Code map

- [`create_code/chain_Hubbard_tbc_response.cc`](create_code/chain_Hubbard_tbc_response.cc):
  runtime-L neutral response/ground solver.
- [`create_code/chain_Hubbard_tbc_spectral.cc`](create_code/chain_Hubbard_tbc_spectral.cc):
  charged-source particle/hole solver, sharing the neutral solver's helpers.
- [`response/mmhk_density_workflow.py`](response/mmhk_density_workflow.py):
  configurable physical response workflow and common-mu selection.
- [`response/run_twist50_density_scan.py`](response/run_twist50_density_scan.py):
  the 50-twist density/neutral campaign.
- [`response/run_single_particle_101_scan.py`](response/run_single_particle_101_scan.py):
  the final direct 101-momentum spectral campaign.
- [`response/replot_eta_response.py`](response/replot_eta_response.py),
  [`response/replot_twist50_response.py`](response/replot_twist50_response.py), and
  [`response/finalize_101_spectra.py`](response/finalize_101_spectra.py):
  postprocessing/replotting, without launching the many-body solver.
- [`response/essler_bethe_ansatz/hubbard_velocities.py`](response/essler_bethe_ansatz/hubbard_velocities.py):
  zero-field BA solver.
- [`skills/mmhk-1d-spin-charge/`](skills/mmhk-1d-spin-charge/):
  reusable workflow, with density, mixing size, twist and response settings.
- [`Gcal/`](Gcal/): historical half-filled Green-function data/notebooks for
  different mixing sizes. The legacy G solver and ground-sector utilities are
  retained under `create_code/`.
- [`note/spin_charge_response_definitions.tex`](note/spin_charge_response_definitions.tex):
  definition note; its source is preserved unchanged.

## Setup and reproduction

Use Python 3.10+ and install `requirements.txt` in your own virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The original machine used `~/.venvs/myenv/bin/python`. Existing plots and arrays
can be inspected without building qbasis. A source build requires compatible
qbasis, GCC/OpenMP, ARPACK ILP64, Intel MKL ILP64 and TBB. qbasis is an external
dependency, not vendored in this repository. The build scripts retain the
original machine's defaults; set the environment variables or adapt library
paths to your installation before compiling. Do not expect committed historical
executable hashes to match a rebuild on another platform.

```bash
export QBASIS_ROOT=/path/to/quantum_basis
export QBASIS_INSTALL=/path/to/qbasis/install
export MKL_ROOT=/path/to/mkl
export TBB_ROOT=/path/to/tbb
export CXX=g++
bash response/build_response.sh
```

`response/build_spectral.sh` records the original Mac/Rosetta build command;
adapt its paths, or use the same link flags as `build_response.sh` with
`chain_Hubbard_tbc_spectral.cc`. Binary snapshots are intentionally excluded.

Plan a configurable neutral calculation **without running it**:

```bash
python skills/mmhk-1d-spin-charge/scripts/run_response.py \
  --project "$PWD" --L 10 --densities 0.7 0.8 0.9 \
  --couplings 0 1 2 3 4 5 6 7 8 9 10 11 12 --twists 50 \
  --eta 0.1 --output response/results/new_campaign
```

Only add `--execute` after reviewing the plan and resource/ensemble settings.
For BA velocities independently of qbasis and without a local copy of the book:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python skills/mmhk-1d-spin-charge/scripts/bethe_reference.py \
  --project "$PWD" --density 0.8 --couplings $(seq 0 0.2 12) \
  --q 0.6283185307179586 --eta 0.1 --nodes 400 --output /tmp/ba_n0p8
```

This writes velocities/Kc and supplementary **BA-informed linear
Luttinger-liquid** response curves (or the exact free lattice response at U=0).
The interacting curves are not exact finite-q Hubbard line shapes. The original
`hubbard_velocities.py` command-line entry point also records a local book-file
hash and therefore expects that file; its imported `solve` function does not.

### Frozen provenance and relocation

Historical manifests intentionally retain original absolute paths and SHA-256
hashes; they have not been rewritten to make an archive appear newly calculated.
Some original replot/runner scripts enforce those paths and executable identities.
After cloning elsewhere, **do not bypass solver-lineage checks to resume an old
campaign**: use a fresh output directory for a new source/build environment.
The repository-wide archive verifier uses relative paths and can be run on a clone:

```bash
python response/verify_repository_archive.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python response/tests/test_density_workflow_helpers.py
```

The archive verifier checks the committed file manifest and current spectral
arrays. Native solver tests additionally require building the solver. Their
saved numerical reports are retained under `response/results/`. Selected
600/900 spectral depth anchors are checks, not a global all-case error bound.

## Archive scope

The initial commit preserves calculation code, numerical `.dat`/`.npz`/`.csv`
outputs and JSON provenance, current/historical plots, the final PPT, notebooks,
and the definition note. It omits verbose qbasis logs, caches, platform-specific
binaries, temporary slide builds, and generated per-index Slurm scripts
(`generate.sh` remains). Downloaded student code, unpublished Overleaf source,
and reference PDFs are not redistributed. Legacy comparison scripts needing
these excluded external inputs are retained as historical diagnostics, not
self-contained reference reproductions. No slow student 10-MMHK calculation
is launched by the archive checks.

See [`response/README.md`](response/README.md),
[`response/PHYSICAL_DENSITY_WORKFLOW.md`](response/PHYSICAL_DENSITY_WORKFLOW.md),
and [`response/STUDENT_CODE_AUDIT.md`](response/STUDENT_CODE_AUDIT.md) for the
development history and observable-definition audit. This archive does not
claim a matched reproduction of the student's Figure 3.
