# Physical MMHK charge/spin density response: implemented and tested

Updated 2026-10-13. **No student module was imported or executed. All new
many-body validation runs used L=6 (plus existing L=4 regression tests), not
L=10.** The implementation is parameterized for L=10 and other even ring sizes.

## Physical definition

At each MMHK reduced-BZ block, calculate the finite Hubbard-ring response using

```
rho_c(q) = sum_j exp(+i*q*j)*(n_j_up+n_j_dn)/sqrt(L)
rho_s(q) = sum_j exp(+i*q*j)*(n_j_up-n_j_dn)/sqrt(L)
q_m = 2*pi*m/L.
```

The orbital vertex is diagonal, delta_(j,j')*exp(+i*q*j)/sqrt(L).
It is not the student's Fourier-mode-occupation vertex exp[i*q*(j-j')].
The spin convention is **2Sz**; divide spin intensity/susceptibility by four
to use Sz. Filling is n=N/(2L), so **n=0.5 is half-filling**.

The primary output is the positive Lehmann dynamical structure factor with
Lorentzian broadening. Excitation energies are neutral E_m-E_0; no mu shift is
added to the frequency. Inelastic and zero-energy elastic weights are separate.

For the full retarded commutator response, the postprocessor computes

```
chi_R(q,z) = sum_m W_m(q)/(z-Omega_m(q))
             - sum_m W_m(-q)/(z+Omega_m(-q)).
```

The negative-frequency counterterm matters at finite eta. Therefore
`physical_density_response.npz` and `retarded_density_susceptibility.npz` are
separate artifacts; -Im(chi_R)/pi is not relabeled as the one-sided broadened S.
At nonzero twist, do not assume S(q)=S(-q) for an individual ring.

## TBC and gauge mapping

The student's orbital-basis hopping matrix is a Hubbard ring with interior
hopping -t and boundary matrix element

```
g_(0,L-1) = -t*exp(+i*theta), theta=L*k.
```

Our qbasis solver uses the uniform hopping gauge, with directed matrix element
`h_(j,j+1)=-t*exp(+i*f)`. For matching this boundary convention, **f=-k**.
The local density operators are invariant under the site-dependent gauge
transformation, so the same rho(q) is used in both gauges. This correspondence
was tested against an independently constructed boundary-bond Hamiltonian.

The new metadata field `boundary_phase_row0_last=-L*phase_per_bond` is the theta
used here. The older `boundary_twist=L*phase_per_bond` field is retained for
compatibility and uses the opposite orientation; do not confuse the two.

## Common chemical potential and ground ensemble

1. For every sampled k, scan E0(k,N) for **all N=0,...,2L**, including all
   translation sectors. SU(2) makes the minimum-|Sz| sector sufficient to find
   the ground energy for each N.
2. Choose a common mu, minimizing E0(k,N)-mu*N independently at each k.
3. Average the physical spectra with normalized k weights, sum(weights)=1.
   This is a per-site normalization. Multiply by the k count to reproduce an
   unnormalized equal-weight k sum.
4. The active spin sectors and complete ground subspaces are included in a
   zero-temperature **unpolarized ground-state trace**, not one arbitrary vector.
   Spin-flipped sectors have identical longitudinal intensities and are counted
   using their multiplicities.

Filling is a staircase on a finite mesh. Default mu selection chooses the stable
plateau closest to the requested filling and reports the achieved value. It does
not pretend that bisection solves a continuous density function or manufacture
fractional canonical particle numbers. At an explicitly supplied mu exactly on
a charge crossing, default behavior rejects ambiguity. Use
`--number-degeneracy uniform` for the equal-weight trace of all grand-ground
states at the crossing; this is not an adjustable mixture to force a density.

Small representative sectors use dense ground-subspace extraction; large ones
use an adaptive small IRAM prefix, with a non-ground sentinel, residual checks,
and orthonormalization of degenerate eigenvectors. General complex ARPACK can
return nonorthogonal right eigenvectors within a degenerate eigenspace; this is
handled explicitly. The forced-IRAM L=6 free-ring fourfold-ground test passes.

No 200-excited-state cutoff, four-orbital-index spectral tensor, or charged-sector
eigenvector cache is used for these neutral responses.

## Reduced-BZ mesh choices

`--nk` is the **actual** number of sampled k blocks.

- `--grid periodic` (production default): distinct periodic-system momenta,
  includes k=0, and counts equivalent RBZ endpoints only once.
- `--grid midpoint`: midpoint quadrature over the reduced BZ, without duplicate
  endpoints. Its finite mesh is offset relative to the periodic global mesh.
- `--grid student`: the supplied student's endpoint-inclusive linspace. Our
  16-point validation sweep deliberately uses this mesh to match its sampling
  while isolating the observable change. Both identified endpoints get equal
  weight here, so this is a discrete sample average, not a claim of converged
  continuum integration.

Increasing nk refines the reduced-BZ integral, but does not refine the pure
single-ring density-transfer grid q_m=2*pi*m/L.

## Off-grid q: distinguish two legitimate, different questions

### Arbitrary-q finite-ring window probe

On a finite ring, the vector sum_j exp(i*q*j)*n_j is definable for any real q.
Unless q is commensurate, it is not an eigenoperator of ring translation. It
decomposes into the supported momenta:

```
a_m(q) = sum_j exp[i*(q-q_m)*j]/L
S_window(q,w) = sum_m |a_m(q)|^2 S(q_m,w).
```

This identity holds for the translation-invariant ground ensemble used here,
including the gauge-equivalent TBC description. It was checked directly against
independent L=6 ED for q=pi/5. The optional `--finite-ring-probe-q` outputs this
well-defined windowed probe, using Dirichlet weights, not interpolation.

It is **not** the physical generic-q response of the full MMHK model, and it
does not create new low-q dispersion points for fitting velocities. Its apparent
continuous variation can be weight redistribution among existing ring channels.

### Physical MMHK transfer between different k blocks

The full density operator sum_p c^dagger_(p+q)*c_p is meaningful outside the
mixing grid. At such q, the two fermions belong to different reduced-BZ blocks.
The probe itself supplies the momentum; the Hamiltonian does not need a matching
scattering term for the response to exist (free fermions are a simple example).

For an independent-block grand-canonical product ensemble, this calculation
requires charged addition/removal spectra from both blocks, combined with the
correct folded/orbital vertices. It is not the one-block neutral resolvent.
The band-HK density response already illustrates nonzero-q response through a
Green-function convolution despite momentum factorization:
https://arxiv.org/abs/2409.07522 (Sec. III).

The missing inter-block interaction vertex can make such a response less useful
for diagnosing the intra-block collective spin/charge separation. This is a
model/interpretation limitation, not mathematical meaninglessness.
**Generic inter-block q is not implemented in this workflow.** It is never
silently replaced by the finite-ring window probe or set to zero.

## Commands

Build the updated qbasis solver:

```bash
bash /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/build_response.sh
```

Production-capable 10-orbital command (not executed during this validation):

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/mmhk_density_workflow.py \
  --L 10 --couplings 8 --nk 16 --grid periodic --target-filling 0.45 \
  --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/physical_mmhk_L10
```

Completed six-orbital validation configuration:

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/mmhk_density_workflow.py \
  --L 6 --couplings 4 6 7 8 12 --nk 16 --grid student --target-filling 0.45 \
  --depth 180 --threads 2 \
  --output /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/physical_mmhk_L6_verified_20261013
```

Completed raw runs are reused only when source/executable hashes and commands
match. After rebuilding, choose a fresh output directory. Original binaries and
sources are archived under `response/versions`; historical data are preserved.

For any completed U directory, create the full retarded susceptibility:

```bash
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/retarded_density_response.py \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/physical_mmhk_L6_verified_20261013/U8
```

The independent reference benchmark refuses L>6 and imports no student modules:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/tests/benchmark_physical_mmhk_L6.py \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/results/physical_mmhk_L6_verified_20261013/U8
```

## Verified results

All five couplings pass independent full-sector Hubbard ED at all 16 TBC values,
for every supported q and both channels. Across the sweep:

- 80 independently checked rings; 1,040 canonical E0(k,N) entries.
- Maximum per-ring spectral discrepancy about 1.1e-12.
- Maximum averaged spectral discrepancy below 1e-13.
- q=0 inelastic charge and longitudinal-spin responses exactly zero in the
  implementation, consistent with the independent reference.
- The arbitrary-q finite-ring window identity agrees with direct ED.
- All active ground subspaces and spin multiplicities agree with direct ED.

| U/t | selected mu/t | achieved n | first charge peak at q_min | first spin peak at q_min |
| --- | --- | --- | --- | --- |
| 4 | 0.922829 | 0.458333 | 2.03 | 1.11 |
| 6 | 1.190660 | 0.447917 | 1.70 | 0.83 |
| 7 | 1.291383 | 0.447917 | 1.63 | 0.74 |
| 8 | 1.369833 | 0.447917 | 1.56 | 0.66 |
| 12 | 1.565770 | 0.447917 | 1.38 | 0.47 |

These are broadened finite-q peak candidates at q_min=pi/3, eta=0.2, not
thermodynamic velocities or a proof of spin-charge separation. Note the achieved
density difference at U=4. The calculation is not claimed to reproduce the old
student Figure 3, which uses a different operator and unspecified historical
ten-orbital settings.

A separate U=0, mu=0 six-orbital test at boundary phases +/-pi includes all 16
grand-ground states, all charge-sector ties and spin multiplicities. It agrees
with independent ED and yields n=0.5. Connected zero-frequency charge/spin
elastic weights remain distinct from their zero inelastic responses.

The old fixed-nine-electron, zero-twist L=10 data remain a canonical single-ring
reference, not the newly implemented k-summed grand-canonical response. Failed
and superseded pilot directories are retained separately; use the verified
directory above for the current results.
