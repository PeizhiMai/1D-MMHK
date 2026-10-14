# Clean reproduction test (14 October 2026)

## Environment and scope

Tests ran on an Apple M4 Pro laptop (24 GiB system RAM), in **Linux/ARM64**
Docker containers with 14 available CPUs and about 7.65 GiB container RAM.
The checkout was an independent `git clone --no-hardlinks` under `/tmp`, not
an overlay of the original working directory. Both Docker targets were built
with `--no-cache`; no original Python venv, qbasis installation, book PDF,
installed skill or Codex runtime was mounted. The numerical candidate was
commit `4f6800be9ccb11577beff7e2332fe3a52677effc`.

The portable public slide exporter subsequently gained a native text legend
and was re-exported/checked; it does not change the scientific plot data.
Later documentation and validation records do not change the solver sources.

## Passed clean-clone checks

1. Pinned qbasis/ARPACK/OpenBLAS/LAPACKE source build, both native solvers,
   and dependency/version checks.
2. **17 tests**: physical momentum folding, invalid configuration rejection,
   six independent dense Fock-space neutral-response cases (including phase,
   doping and an operator-comparison check), and three independent dense
   particle/hole cases. No old native executable or numerical output was
   used as an exact-diagonalization shortcut.
3. End-to-end **L=4** smoke: four cases (U=0,4; n=.5,.75), common-mu scan,
   actual complete spin/momentum grounds, neutral and charged calculations,
   BA, response/spectral figures, tables and the five-slide test deck.
   The smoke's density tolerance is deliberately 0.13, not the study's 0.004.
4. Replot that new run in the **separate plot-only image**, which has no
   qbasis/native solver installation.
5. Replot **all 39 archived L=10 cases** in that solver-free image. The 39
   response/velocity rows agree exactly for mu, actual density and both peak
   velocity definitions. BA reference differences are at most
   **2.1e-13**. Reconstructed neutral low-window arrays agree exactly.
   All individual spectral maps and the density/BA panels are regenerated.
6. Full relative-file manifest verification: 92,334 files / 1,168,408,369
   bytes at the candidate commit; 39 spectra with 101 exact k points and
   eta=.1; 183 fine-grid BA rows; the unchanged old four-slide deck plus the
   seven-slide deck (**42 native charts and 42 embedded workbooks**).
7. Public exported deck: seven slides in the requested density order,
   scientific PNGs embedded, editable headings and speaker notes. All seven
   slides rendered and visually inspected. The original editable-chart
   seven-slide deck remains byte-identical (SHA-256
   `20187e3b17684c53e8235f61611274c97b87b7b0a70b803925764f7e7f02b045`).

Evidence: `validation/clean_clone.json`,
`validation/archive_verification_candidate.json`, and
`validation/versions_and_solver_hashes.txt`. Manifest file counts in these
records refer to the tested candidate, not a self-referential count after
adding these reports. The verifier calculates the final count on each clone.

## Additional fresh L=10 physics tests

- Independent **44,100-state sparse ED/Krylov** test at U=8, Ne=8, zero twist,
  q_min: ground-energy agreement at roundoff; maximum charge/spin spectral
  errors **4.3e-10 / 8.5e-8** against the depth-300 native calculation,
  eta=.1. This is a canonical-sector anchor, not a replacement for the full
  spin-multiplet trace. Report: `validation/l10_independent_sparse.json`.
- Fresh **U=8, 50 twists, all three target densities, k=101** complete
  common-mu pipeline, using only the new portable solvers; all response and
  spectral stages, plots and seven-slide export completed. The densities,
  mu and energies agree with the frozen study within the stated tolerances.
  Maximum single-particle old/new discrepancy is **5.43e-6**, or **3.0e-6 of
  the global spectral maximum**. Source normalization and first moments
  remain at roundoff. Report: `validation/l10_u8_twist50_regression.json`.
- The tiny finite-depth spectral differences are not silently hidden or
  attributed to changed physics. A targeted **depth-900** triplet-trace
  anchor at the largest n=.8 discrepancy gives new600/new900 difference
  **2.98e-6**, consistent with finite-depth Lanczos sensitivity. The archive
  already records depth-600/900 and duplicate-inversion differences. Report:
  `validation/u8_depth900_anchor.json`. This selected anchor is not a global
  all-case convergence proof. Regression uses 1e-4 of the global spectral
  maximum, not relative error on a vanishing tail; independent ED and native
  sum-rule gates remain separate and tighter.

## Completed full fresh study

The entire default study also finished in the clean portable container: **all
39 cases**, U=0..12, n=.9/.8/.7, 50 density twists, 101 momenta, eta=.1;
new scans and sources, all figures/tables, fine BA reference and seven slides.
No archived numerical data were used during this calculation. Measured wall
time was **112.8 minutes** with workers=4, threads=1 on the above laptop.

Comparison to the frozen archive after solving:

| Quantity | Maximum absolute difference |
| --- | --- |
| Ground energies over all number/twist tables | 4.31e-14 |
| Dominant charge / spin velocities | 2.48e-13 / 1.91e-13 |
| Charge / spin response intensity | 2.36e-10 / 3.68e-7 |
| Single-particle intensity | 7.34e-5 |
| New source normalization / first-moment errors | 9.77e-15 / 9.77e-14 |

The maximum single-particle difference is **4.67e-5 of the global spectral
maximum**, at U=12, n=.7, omega=13.58 (outside the displayed ymax=10).
A new depth-900 calculation at that exact momentum gives new600/new900 and
old600/new900 differences **3.85e-5 / 4.04e-5**, consistent with finite-depth
sensitivity. This remains an anchor, not a global convergence proof.
Neutral regression uses **1e-6 of the global response maximum** and spectral
regression **1e-4 of the global spectral maximum**. Using relative error on
small high-energy tails would mistake tiny finite-depth differences for an
operator change. Independent ED, density/energy, and sum-rule checks are
separate and tighter. No data were renormalized to force agreement. After moving the entire fresh
run from `/tmp` into the project, a solver-free replot at **eta=.2** regenerated
all 39 response and single-particle cases from the relative raw plans, verifying
relocation and rebroadening (`validation/relocation_rebroadening.json`).

Records: `validation/full_study_summary.json`,
`validation/full_study_regression.json`, `validation/full_study_native_timing.json`,
`validation/u12_depth900_worst_anchor.json`. The public exporter-generated
`validation/fresh_study_7slides.pptx` is included in addition to the unchanged
original native-chart deck. Full raw test outputs are retained locally under
`runs/reproducibility-validation-full-study` (ignored by Git), rather than
rewriting or duplicating the frozen archive.

## Repeat the checks

Run `MMHK_NO_CACHE=1 bash scripts/clean_test.sh` from a clone. Optional larger
checks (the latter needs the archive for comparison **after** solving):

```bash
docker run --rm -v "$PWD:/workspace" --entrypoint python mmhk-1d:repro \
  reproducibility/validate_l10.py --output runs/l10-independent
bash scripts/run_docker.sh calculate --config configs/l10-u8-validation.json --output runs/u8-validation
docker run --rm -v "$PWD:/workspace" --entrypoint python mmhk-1d:plots \
  reproducibility/compare_archive.py --input runs/u8-validation --couplings 8 --output runs/u8-regression.json
```

Linux/AMD64 CI smoke is configured but **not claimed as locally tested**.
Native macOS/Windows builds are not certified here. Docker is the tested laptop
route. Numerical agreement within tolerances is the reproducibility contract,
not identical ARPACK vectors, executables, PNGs or PPT ZIP bytes. No slow
student solver was run. The existing LaTeX note and numerical archive were
not edited.
