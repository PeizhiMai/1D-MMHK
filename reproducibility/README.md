# Reproduce on a laptop

The maintained entry point is `reproduce.py`. The historical `response/*.py`
commands and slide builders remain as provenance; several assume the original
machine and are **not** the fresh-clone interface.

## Docker: recommended on macOS, Windows/WSL2, and Linux

Install Docker Desktop or Docker Engine and start it. No Intel MKL, Rosetta,
installed Codex skill, local qbasis checkout, or book PDF is required.
On macOS clone into `~/src`, not privacy-protected Desktop/Documents folders.
On Windows use the commands inside WSL2 with Docker integration enabled.
Use a shell with bash and git:

```bash
mkdir -p ~/src && cd ~/src
git clone https://github.com/PeizhiMai/1D-MMHK.git
cd 1D-MMHK
bash scripts/setup.sh --docker
bash scripts/run_docker.sh doctor
bash scripts/run_docker.sh calculate --config configs/smoke.json --output runs/smoke
bash scripts/run_docker.sh plan --config configs/study.json
bash scripts/run_docker.sh calculate --config configs/study.json --output runs/study
```

The full command starts from **new** ground-energy scans, tunes one common mu
per coupling/density, calculates neutral sources for complete grand ground
subspaces, calculates particle/hole sources at the directly requested momenta,
solves the BA equations, and writes all figures/tables and the presentation.
It does not use any old numerical archive as a calculation shortcut. U=0 is
calculated analytically using the full free grand-ground projector.

### Plot only: no many-body solver installed or launched

```bash
bash scripts/setup.sh --plots
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plot --output runs/archive-plots
```

This smaller image contains Python, Node and public plotting/slide dependencies,
but **no qbasis or native solver**. It reads the committed neutral poles and 39
single-particle arrays, checks relative SHA-256 provenance, recomputes the fine
BA reference (not ED), and generates the figures and seven slides. It also
works with newly calculated portable runs:

```bash
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plot \
  --input runs/study --output runs/study-replot
```

Both modes write under `runs/`, leaving the scientific archive and open TeX
source untouched. `--no-slides` skips PowerPoint export. Fresh portable spectral
runs include continued-fraction plans and support `plot --eta 0.2`; the frozen
101-point archive only exposes eta=0.1 through this interface, so other eta
values are rejected, never merely relabeled.

## Configuration and resources

`configs/study.json` is the requested 10-MMHK study: physical densities
0.9/0.8/0.7, U/t=0..12, 50 periodic reduced-BZ twists, density tolerance 0.004,
q_min=2*pi/10, response/spectral Lanczos depths 300/600, eta=0.1, 101 directly
calculated k in [0,pi], stored omega [-6,16], displayed omega [-6,10], and BA
step 0.2 with 400 quadrature nodes. `configs/smoke.json` is a **toy L=4 test**,
with a deliberately larger 0.13 density tolerance; do not mistake it for the
study's scientific accuracy policy.

Copy a config, edit L (even >=4), densities (0<n<1), twists (1 for zero twist,
or even >=4), U list, depths, broadening, grids or worker/thread counts. No
recompile is needed for different L. Half filling is **n=Ne/L=1**, not 0.5.
The interface currently computes q_min and assumes zero field/temperature.
Unsupported domains fail explicitly. Increasing twists improves density
resolution but does not decrease q_min. Increasing L costs exponentially:
`binom(L,L/2)^2`. The default laptop gate is one million canonical states,
and a conservative LP64 CSR index gate is also enforced.

Provide 8 GB RAM to Docker (16 GB system RAM is a sensible starting point) and
at least 10 GB of free disk (more for new raw sources/logs). Workers=4, threads=1 avoids oversubscription; use workers=2
on small laptops. Runtime is hardware/backend dependent. The clone includes
about 1.2 GB of scientific files, with many small raw-source files. Allow time
for checkout. Larger orbital counts require a separate resource assessment.
See `CLEAN_TEST_REPORT.md` for measured tests, not a claim that every laptop
will finish in the same time.

The 101 physical k points are not 101 density twists. They fold to 20 distinct
RBZ phases, calculated via 11 inversion representatives of the **complete
unpolarized ground trace**. No A(k,omega) interpolation occurs. The common mu
is fixed by the 50-twist density grid, not retuned on the spectral grid.

The neutral response averages complete spin/momentum ground subspaces before
peak extraction. A stable mu plateau nearest to the target is accepted only
within the explicit density tolerance, and its actual density is reported.
No arbitrary number-sector mixture is used to fabricate a target. A true U=0
Fermi-level grand-ground projector is retained. Both channels use coherent
physical densities, spin=2Sz, normalization 1/sqrt(L), and neutral Em-Eg with
no mu shift. Charged frequencies use the common mu and per-spin unit weight.

## Outputs and restart

```
runs/study/config.json                  source/build/settings identity
runs/study/U*/ground_energy_table.npz   complete number scans + symmetry checks
runs/study/U*/n*/metadata.json          mu, actual density, full ground trace
runs/study/U*/n*/averaged_qmin_poles.npz neutral energies/weights
runs/study/spectral/U*/n*/Akw.npz        direct k, omega, per-spin A
runs/study/spectral/U*/n*/spectral_plan.json relative raw fractions for rebroaden
runs/study/plots/                       PNG/SVG panels, individual maps, CSVs
runs/study/plots/spin_charge_spectral.pptx
```

Panels include all 13 U values per density. Full-range neutral panels and
broadening-sensitivity tables supplement the low-energy plots. Velocity plots
show both channels and BA in the same panel for each density. Peak extraction
selects resolved local maxima in 0 < omega/t <= 5; it does not identify a
high-energy Hubbard-band maximum as a sound mode. The MMHK
peak/q values are finite-q, broadening-dependent estimators; BA velocities are
thermodynamic quantities, not exact finite-q peak/line-shape predictions.

Native and postprocessing stages can be separated:

```bash
bash scripts/run_docker.sh calculate --config configs/study.json --output runs/study --stage response
bash scripts/run_docker.sh calculate --config configs/study.json --output runs/study --stage spectral
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plot --input runs/study --output runs/study/plots
```

Rerunning the identical command reuses **only this run's** completed, matching
source/build/settings outputs. Moving a run does not change cache identities.
Do not copy archived native jobs into it. Changed solver/settings need a new
output root. A native job with metadata but no COMPLETE marker is treated as
interrupted and requires inspection/removal of that specific job directory;
it is not silently accepted. Changed plot settings need a new plot directory.
If density is unattainable within tolerance, computation stops rather than
manufacturing spectra at the nominal density.

## Presentation portability

The committed, previously verified seven-slide deck is
`response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9/analysis_ba_step0p2_eta0p1/presentation/spin_charge_spectral_10MMHK_eta0p1_7slides.pptx`.
It contains 42 editable native charts and the three spectral-map images.
Its original Codex-based builders are retained unchanged as historical sources.

The new public PptxGenJS exporter reproduces the **scientific plots** on seven
slides with editable headings/notes, embedding the figures as PNGs; it does
not require the private Codex chart renderer. Thus its editable object types,
fonts and layout are not byte-identical to the historical deck. Slide order:
velocity comparison; responses/spectra for n=.9; responses/spectra for n=.8;
responses/spectra for n=.7. Other density lists produce 1+2*density_count slides.
Scientific input hashes and export settings accompany each deck.

## Native setup (optional, Linux first)

Docker is the tested portable route. For a native Debian/Ubuntu build install
Python 3.12 with venv, GCC/GFortran, CMake >=3.24, Boost, OpenBLAS/LAPACKE and
Node/npm. For example:

```bash
sudo apt-get update
sudo apt-get install g++ gfortran cmake make libboost-dev libopenblas-dev liblapacke-dev python3-venv nodejs npm
PYTHON=python3.12 bash scripts/setup.sh --native
.repro-venv/bin/python reproduce.py doctor
.repro-venv/bin/python reproduce.py calculate --config configs/smoke.json --output runs/native-smoke
```

CMake finds system LAPACKE/OpenBLAS (LP64) and OpenMP. Dependencies are fetched
with SHA-256 checks. `CC`, `CXX`, `FC`, `CMAKE_PREFIX_PATH` and `BUILD_JOBS` can
be set by the caller. Native macOS needs a consistent GCC/GFortran/OpenMP and
LAPACKE installation; Apple Clang alone is insufficient. Native macOS/Windows
are **not** clean-tested here; use Docker rather than the historical MKL paths.
For plot-only native Python, install `requirements-repro.txt` in a Python 3.12
venv and run `npm ci --ignore-scripts` in `reproducibility/slides`.

## Tests

```bash
# inside the Docker full image, preserving its fresh solver installation
# (run from the clone root)
docker run --rm -v "$PWD:/workspace" --entrypoint python mmhk-1d:repro \
  -m pytest reproducibility/tests -q
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plot --output runs/test-archive-plots
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plan --config configs/study.json
# Verify the complete frozen archive, including the 7-slide deck:
docker run --rm -v "$PWD:/workspace" --entrypoint python mmhk-1d:plots \
  response/verify_repository_archive.py
```

The tests use independent Fock-space operators for neutral and charged spectra,
check phases, ground traces, source normalization and moments. Docker binds
only the clean checkout, never the original home/dependency directories.
Numerical reproducibility means agreement within documented tolerances, not
identical random ARPACK vectors, platform-dependent executables or PNG/PPT ZIP
bytes. The frozen inputs are authenticated without rewriting their old hashes.

Use `MMHK_NO_CACHE=1 bash scripts/clean_test.sh` for a cold-build check; the
U=8 50-twist regression configuration and independent L10 checks are described
in [the measured clean test report](CLEAN_TEST_REPORT.md).
