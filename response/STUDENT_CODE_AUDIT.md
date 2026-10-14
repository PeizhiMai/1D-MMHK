# Student-code audit: 2026-10-13

**No student solver was run or imported, at either 10 or 6 orbitals.** The
inspection used a read-only source download and AST/text inspection. An
independent 6x6 one-body algebra check was also performed; it is not a many-body
MMHK calculation. No new qbasis many-body runs were made during this audit.

Source: https://drive.google.com/drive/folders/1gzZ0ay-z-hGRC53dtZhbz-EPU4lWylX0

The downloaded seven Python source files and their hashes are recorded in
`/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/student_source_20261013/download_manifest.json`.
Only sources/helpers and two small pre-existing four-orbital outputs (one NPZ
and one PNG) were downloaded. No ten-orbital ED/cache data were downloaded.
Compiled caches were not imported or executed.

## 1. Reduced-BZ sum: present

In the supplied source:

- `params.py:8-9` creates Q=2*pi/N_ORBITALS and
  `K_GRID=linspace(-Q/2,Q/2,Nk//2,endpoint=True)`.
- Default input Nk=32 therefore produces **16** actual reduced-BZ points.
- `spin_charge_main_v10.py:62-64` obtains the ED cache on the whole k grid.
- `spin_charge/spin_charge_v5.py:141-143` computes responses at every k and q.
- `spin_charge_main_v10.py:88,108` explicitly sums the response over k.

This is not the single k=0, fixed-nine-electron ring calculated previously.
The main driver as supplied currently sets **N_ORBITALS=4**, U=8, x_target=0.45.
Thus the exact historical N=10 Figure 3 generating configuration is not included
as a separate launch file. Switching this workflow to N_ORBITALS=10 still retains
the k integration and the other conventions documented below.

The existing output in the shared data folder confirms the supplied four-orbital
configuration: N_ORB=4, U=8, x=0.453125, eta=0.2. Its saved figure title reports
mu approximately 1.209 (rounded). This is NOT a chemical potential for the
ten-orbital Figure 3 benchmark. Its files are saved under
`response/student_source_20261013/existing_data/`.

The raw k sum is not divided by Nk or multiplied by integration weights. This
affects amplitude comparisons. Both endpoints are included with equal weight;
the default even-length k grid contains no k=0 point. The endpoint blocks have
equivalent boundary phases and are counted twice in a simple rectangle sum.

## 2. Chemical potential: present, at zero temperature

`spin_charge/ED_cache_canonical_v2.py:48-74`:

1. Computes E0(k,Nup,Ndn) across all (Nup,Ndn) sectors.
2. Selects the sector minimizing E0(k,Nup,Ndn)-mu*(Nup+Ndn) independently at each k.
3. Tunes one common mu by bisection to target the average filling.

The achieved filling is

```
x = sum_k N_ground(k) / (2 * Nk * N_ORBITALS).
```

`ED_cache_canonical_v2.py:88-104,131-134` saves this achieved average. For 10
orbitals and 16 actual k points, the reference values 0.453125 and 0.446875
correspond to total counts 145/320 and 143/320. They are not single-block
noninteger electron counts and do not imply a thermal ensemble.

The supplied response uses ground states and averages degenerate ground vectors
within the selected (Nup,Ndn) sector (`spin_charge_v5.py:57,70,107`). There are no
Boltzmann weights or temperature parameters. It is a zero-temperature
grand-canonical sector-selection calculation, not finite-temperature ED.

For each neutral response, mu cancels from the excitation energy within its
fixed-N block. Its role is choosing N(k), not shifting omega after the response.

The filling function is a staircase; the bisection call can terminate at a
chemical-potential jump without achieving x_target exactly. The actual returned
x must be used. Comparing with a fixed N=9 at every k is not equivalent.

## 3. "Zero twist" is not a single untwisted ring

`fock_helpers/orbital_unitary.py:4-38` constructs

```
U_(a,j)(k) = exp[i*(k+2*pi*a/L)*j] / sqrt(L)
g_orb(k)  = U(k)^dagger diag[-2*cos(k+2*pi*a/L)] U(k).
```

This orbital gauge has interior nearest-neighbor matrix entries -t and boundary
entries

```
g_(0,L-1)(k) = -t*exp(+i*L*k),
g_(L-1,0)(k) = -t*exp(-i*L*k).
```

There is no additional external twist parameter in the supplied driver. However,
the reduced-BZ k scan corresponds to block boundary phases theta=L*k ranging
from -pi to +pi. Only k=0 is the genuinely untwisted 10-site ring. No k=0 point
is present in the default 16-point endpoint-inclusive grid.

The same block Hamiltonian can be put into the uniform-phase gauge of our qbasis
code with phase_per_bond=-k and c_uniform(j)=exp(+i*k*j)*c_orb(j). The physical
local-density transfer operator is invariant under this gauge change, whereas
the Fourier-mode occupation used below is not.

## 4. Actual response operator: not the physical density-transfer operator

`spin_charge/spin_charge_v5.py:64,74-103` uses every orbital pair, with

```
vertex(alpha,alpha') = exp[i*q*(alpha-alpha')],
operator term        = c^dagger_(alpha',sigma) c_(alpha,sigma).
```

Lines 105 and 171-173 retain all cross terms and sum them. The implementation is
**coherent**, contrary to the exploratory interpretation of the written note's
outside-the-modulus orbital sum used for our earlier optional diagnostic.

Define

```
d_(q,sigma) = sum_alpha exp(+i*q*alpha)*c_(alpha,sigma)/sqrt(L).
```

The actual operators implemented are exactly

```
R_charge(q) = L * [d^dagger_(q,up)d_(q,up) + d^dagger_(q,down)d_(q,down)],
R_spin(q)   = L * [d^dagger_(q,up)d_(q,up) - d^dagger_(q,down)d_(q,down)].
```

They probe **occupation fluctuations of one Fourier orbital mode**. Here q is
the mode label, not a fixed momentum transfer. The vertex is a rank-one
projector in orbital space.

For the physical density response at supported transfer q=2*pi*m/L, the orbital
vertex instead is

```
Gamma_(alpha',alpha)(q) = delta_(alpha',alpha)*exp(+i*q*alpha)/sqrt(L),
rho_charge/spin(q)     = sum_alpha exp(+i*q*alpha)
                         * [n_(alpha,up) +/- n_(alpha,down)]/sqrt(L).
```

In the folded-momentum basis this is the coherent sum of bilinears whose
momenta differ by q. Merely restricting the current phase-difference vertex to
alpha=alpha' would give phase=1 and would NOT implement the required operator.

At q=0, R(q) is the occupation of the zero Fourier mode, not conserved total N
or 2Sz. Finite-frequency q=0 weight is therefore possible for R. It does not by
itself demonstrate a particle-number conservation error. The student cannot
interpret q-dependent R-mode peak slopes as physical density-transfer spin/charge
velocities without an additional derivation.

In the uniform hopping gauge, this R source becomes the occupation of the
generally noncommensurate Fourier wave number q-k. It can involve many momentum
transfers. Consequently the physical fixed-transfer U=0 bound
4t*sin(|q|/2) used in our previous check does NOT apply to this student operator.
The extra U=0 reference peaks are now explicable as a different observable;
they are not evidence that this code violates that bound for its own observable.

## 5. Why the calculation is slow

For N_ORBITALS=10 with default Nk input 32:

- The initial sector scan calls **121*16=1936** ground-state solves.
- Each basis builder scans all 2^20 bit strings, repeatedly
  (`MMHK/block_hamiltonian_v2.py:37-50`).
- The active-sector solver asks ARPACK for up to **200 eigenvectors**, with no
  translation reduction (`ED_cache_canonical_v2.py:21-39`).
- It diagonalizes addition/removal sectors too, although applying a bilinear
  only requires intermediate sector bases, not their eigenvectors.
- Each channel allocates a complex tensor of shape
  (16,6,10,10,10,10,400), **6.144 GB** per channel, before eigenvector caches and
  intermediate tensors (`spin_charge_v5.py:134-137`). The driver retains both
  channel arrays, so these alone can occupy 12.288 GB.

The scalar total response can be computed directly from R(q)|ground> without
forming the four-orbital-index spectral tensor. For physical density responses,
use rho(q)|ground> instead. Lanczos continued fractions avoid computing 200
excited eigenvectors. Cache basis/Hamiltonian data, avoid duplicated endpoints,
and match the k/sector/ground-state prescription explicitly before comparing.

## Conclusion and next implementation target

The previously produced qbasis spectra remain validated physical fixed-N,
single-block density responses. They are **not a reproduction benchmark** of
this student's full k-summed, mu-selected, Fourier-mode-occupation spectra.

A matched comparison needs the exact Figure 3 launch settings and then either:

1. Reproduce R(q), its k sum, gauge, mu selection, state truncation and ground
   averaging faithfully; label it a Fourier-mode-occupation response.
2. Compute the physical rho(q) response with the same common mu and k sum, and
   use that quantity to investigate spin-charge separation.

Any execution of the student's reference solver should be limited to **6-MMHK**
as requested. This source inspection was sufficient to establish the present
differences; no such execution was necessary.
