# 1D 10-MMHK spin and charge responses

**Current physical workflow:** `mmhk_density_workflow.py` now implements the
supported-transfer physical density response with common-mu selection and
reduced-BZ averaging. It is benchmarked against finite Hubbard TBC rings at
L=6, U=4,6,7,8,12. See `PHYSICAL_DENSITY_WORKFLOW.md` for commands, verified
outputs, and the distinction between a finite-ring off-grid window probe and
the not-yet-implemented generic-q inter-block MMHK response.

**2026-10-13 source-audit clarification:** these results are a validated physical
density response of one canonical Hubbard-ring/MMHK block, **not** a matched
reproduction of the student's Figure 3. The supplied student code uses a
reduced-BZ k sum, common-mu sector selection, and a different Fourier-mode
occupation operator. See `STUDENT_CODE_AUDIT.md` for exact source references.
The U=0 fixed-transfer bound does not apply to that different operator.

Implemented and tested locally on 2026-10-12. The existing Green-function sources
and notebooks were not modified. The Overleaf project was read/downloaded only.

## Model and conventions

The implementation uses the same dual-space Hubbard-ring Hamiltonian as the
existing twisted-boundary qbasis Green-function program:

```
H = -t sum_(j,sigma) [exp(+i phase) c^dagger_(j,sigma) c_(j+1,sigma) + h.c.]
    + U sum_j n_(j,up) n_(j,down).
```

- L=10 mixing orbitals by default, t=1, U=8.
- `filling = (Nup + Ndn)/(2 L)`. **0.5 is half-filling.**
- Figure 3 uses nominal filling 0.45: in a single canonical zero-twist block this
  means N=9, Nup=5, Ndn=4. It does not mean N=4 or N=5.
- `--phase` is a **per-bond** phase in radians; boundary twist is L times that
  phase. `--phase 0` is exactly zero; the old 1e-5 offset is not copied.
- Supported transfer momenta are q=2*pi*m/L, m=0,...,L-1. This is not an
  arbitrary continuous physical-k single-particle calculation.
- With c_K=sum_j exp(+i*K*j)c_j/sqrt(L), the bare energies are
  epsilon(K)=-2t*cos(K-phase). The Fourier mapping is independently checked.
- This solver evaluates **one canonical MMHK block**. It does not silently
  average different reduced-BZ k blocks, twists, particle sectors, or temperatures.
- Fractional filling that gives noninteger N is rejected. Arbitrary total
  fillings across many MMHK blocks require an explicit particle-sector/ensemble
  prescription, not rounding a single block's particle count.

## Physical observables (default)

```
rho_charge(q) = sum_j exp(+i q j) [n_(j,up)+n_(j,down)] / sqrt(L)
rho_spin(q)   = sum_j exp(+i q j) [n_(j,up)-n_(j,down)] / sqrt(L)
S_nu(q,w)    = sum_m |<m|rho_nu(q)|g>|^2 delta[w-(Em-Eg)].
```

Spin here is **2Sz**, not Sz; divide spin intensity by four for Sz normalization.
Unlike Gp/Gh, the target sector keeps Nup/Ndn unchanged and has Q_target=Q_ground+m.

The global ground energy is located by checking all translation sectors within
the chosen Nup/Ndn sector. Degenerate ground momenta are averaged with equal
weights. The default nine-electron cases have Q=2 and Q=8; spin exchange maps
Nup=5,Ndn=4 to its identical longitudinal-response partner Nup=4,Ndn=5.
The original 2026-10-12 version obtained one vector per tied Q. The updated
solver resolves degenerate ground subspaces within each Q using dense extraction
or an adaptive small IRAM prefix and orthonormalization. The new full workflow
also includes active spin sectors in an unpolarized ground trace. Historical
sources and binaries are preserved under `response/versions/`.

At q=0, charge and longitudinal spin are conserved totals. Their inelastic
responses are exactly zero. Analytic elastic weights are N^2/L and
(Nup-Ndn)^2/L, respectively. At other q, any poles within 1e-8 of zero are also
separated into elastic weight. `spectra.npz` contains only the inelastic part.

The positive Lehmann spectrum is broadened with eta/pi / [(w-Em+Eg)^2+eta^2].
No artificial chemical-potential shift is used: these are neutral excitations.
Negative-frequency retarded counterterms are not included in this one-sided S.

## Diagnostic alternative: the orbital-pair intensity sum

The note's Eq.(2) writes the orbital sum **outside** the squared matrix element.
That is different from the physical coherent density operator above. The
optional `--observable orbital-pairs` explicitly computes:

```
A_(K,nu)(q) = c^dagger_(K,up)c_(K+q,up)
              +/- c^dagger_(K,down)c_(K+q,down)
S_pair(q,w) = (1/L) sum_(K,m) |<m|A_(K,nu)(q)|g>|^2 delta[w-(Em-Eg)].
```

The required momentum-transfer constraint is kept. This mode is named separately
and is **not** presented as the physical density response. It can have finite
q=0 inelastic weight because each n_K is not conserved at U>0.
Both modes pass independent small-system ED tests.

This was an exploratory interpretation of the written equation. Inspection of
the actual student source on 2026-10-13 shows that it instead sums all orbital
amplitudes coherently with a phase-difference vertex, producing L*d_q^dagger*d_q.
Neither this optional pair-sum diagnostic nor the default physical rho(q) source
is the actual student source. The diagnostic is retained as a separate observable,
not as a reproduction of the student's implementation.

## Build and run

Use the installed qbasis, GCC 13, Intel MKL ILP64, and TBB. On this ARM Mac the
existing library and the linked executable are x86_64 and run through Rosetta.
The build script fixes the legacy bare qbasis dylib install-name in the new
executable, without modifying the installed library. Override `CXX`,
`QBASIS_ROOT`, `QBASIS_INSTALL`, `MKL_ROOT`, or `TBB_ROOT` for another installation.

```bash
bash /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/build_response.sh

OMP_NUM_THREADS=4 MKL_NUM_THREADS=1 \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/bin/chain_Hubbard_tbc_response.x \
  --L 10 --U 8 --filling 0.45 --phase 0 --depth 300 \
  --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/my_U8_run

MPLCONFIGDIR=/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/.mplconfig \
  ~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/postprocess_response.py \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/my_U8_run
```

Use `--filling 0.5` for half-filling, or explicit `--nup 5 --ndn 5`.
`--q 1` computes only q_min. Default is all ten transfer momenta.
`--ground-q Q` intentionally restricts the ground search; default scans every Q.
Completed output directories are protected from overwriting.

Run the requested five-coupling benchmark:

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/run_zero_twist_benchmark.py
```

The runner uses four OpenMP threads and one MKL thread, checks source/executable
hashes before reusing completed data, and produces all plots locally. Supply
`--output NEW_DIRECTORY` after rebuilding/changing the solver. It does not submit
Slurm jobs. Individual full-q L=10 runs took about 6 seconds here.

## Files produced

- `metadata.json`: filling, twist convention, ground momenta, energies/residuals,
  observable, and measured q indices.
- `ground_sectors.csv`: ground energy in every searched Q sector.
- `charge/spin_gsQ*_q*.dat`: raw qbasis norms and a/b Lanczos coefficients.
  Pair mode has an extra `_orbK` tag. q=0 physical elastic weight is analytic.
- `COMPLETE`: written only after the run finishes successfully.
- `spectra.npz`, `poles.csv`, `summary.json`: inelastic spectra, raw Ritz poles,
  elastic weights, mode candidates and numerical checks.
- `response_map.png`, `qmin_response.png`: response plots.
- `manifest.json`, `solver_source_snapshot.cc`: benchmark configuration/provenance.
- `figure3_comparison.csv/json/png`: reference peak and normalized-shape comparison.

Unmeasured q values in a partial run are NaN, not physical zeros. Postprocessing
checks that every required ground/orbital source exists. Ritz poles can be
duplicated by finite-precision unreorthogonalized Lanczos; do not count duplicate
poles as independent eigenstates. The native depth argument includes its unused
b0 slot: depth=300 usually returns 299 a coefficients.

## Verification and benchmark status

Six independent four-orbital dense Fock-space checks passed: half-filled U=8,
doped U=8, dilute U=0, nonzero-twist U=4, and two orbital-pair cases. They check
ground energies, all q spectra, total weights, and elastic weights. Maximum
spectral discrepancies were below 2e-13. Invalid/noninteger input is rejected.

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/tests/test_response.py
```

For the five requested L=10 runs, source norm / zeroth-moment and first-moment
errors are below 3e-14. Inversion symmetry differences are below 6e-7. For U=8,
doubling the q_min Lanczos depth from 300 to 600 changes charge/spin spectra by
less than 4e-11 relative to their maxima and leaves the first peaks unchanged.
Ground-vector residuals are about 2e-12. See the saved validation/convergence JSON.

An additional independent **L=10** check builds the unreduced 52,920-state
Hamiltonian (629,160 sparse entries), obtains its ground doublet with SciPy
ARPACK, and reconstructs q_min responses without qbasis operators or translation
reduction. It reproduces E0=-4.803879095562 and both first peaks. Maximum spectral
errors are 3.5e-15 (charge) and 1.3e-11 (spin). Its 600-step full-space chain is
needed to resolve sources that mix the two degenerate momentum sectors; a
200-step exploratory chain was insufficient in the weak high-energy spin tail.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/tests/test_L10_sparse.py
```

### Historical canonical comparison: not a matched Figure 3 benchmark

The downloaded original Figure 3 uses nominal filling 0.45 and q_min=2*pi/10.
The packaged NPZ arrays specify eta=0.2; the U=8 array's actual filling is
0.453125. U=4 and U=8 curves are read from those arrays; U=6,7,12 normalized
shapes are digitized from the original PNGs. The peak targets below come from
the **printed figure legends**, not graphical estimates. Absolute amplitudes
are not compared because the source normalization/aggregation is unspecified.

| U/t | qbasis first charge peak | Figure 3 charge | qbasis first spin peak | Figure 3 spin |
| --- | --- | --- | --- | --- |
| 4 | 0.648 | 1.343 | 0.683 | 0.792 |
| 6 | 0.490 | 1.095 | 0.526 | 0.654 |
| 7 | 0.444 | 0.952 | 0.468 | 1.454 |
| 8 | 0.411 | 0.847 | 0.421 | 1.068 |
| 12 | 0.341 | 0.639 | 0.299 | 0.789 |

No pair of first peaks agrees within two reference frequency-grid cells. The
charge global maximum is a different peak (e.g. 0.927 at U=8); that distinction
is preserved in the comparison CSV. It would be misleading to drop the lower
peak simply to match the reference. The optional Eq.(2) pair sum at U=8 also
does not reproduce the reference first peaks.

Additional controls:

- At half-filling U=8, the q_min spin pole is 0.451568; the first significant
  charge pole is 5.347434. These are finite-q neutral excitation energies, not
  thermodynamic gaps or fitted velocities.
- The U=0 reference supplies a useful definition check. For the supplied
  epsilon(k)=-2t*cos(k), all physical particle-hole excitations at q_min obey
  w <= 4t*sin(pi/10)=1.236068t, independently of filling or twist averaging.
  Its labelled q_min curves have extra local maxima at approximately 1.977,
  3.135 and 3.742. Broadening positive poles below the bound cannot create new
  local maxima beyond it. The 2026-10-13 source audit resolves this: the student
  uses a different Fourier-mode-occupation observable with no fixed momentum
  transfer. This bound does not apply to that observable and these extra peaks
  alone do not imply an error in its Lehmann calculation.

The student solver was subsequently obtained and inspected on 2026-10-13,
without execution. It sums over reduced-BZ k, chooses particle sectors using a
common chemical potential, and computes a coherent Fourier-mode-occupation
response rather than physical density at momentum transfer q. The supplied
driver defaults to four orbitals, so the exact historical ten-orbital Figure 3
launch settings still need confirmation. The source audit explains why the
historical comparison is not like-for-like; it does not establish a numerical
failure of the student's Lehmann solver.

The official five-coupling results are in
`/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/zero_twist_benchmark`.
Earlier exploratory/pilot directories are retained separately and are not the
authoritative benchmark. The half-filling, free-band, convergence, and operator
diagnostics are controls, not additional density scans or cluster submissions.
