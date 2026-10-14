# Dependencies, pins and compatibility scope

- Docker base: official Python 3.12 slim-bookworm, immutable multiarchitecture
  digest in `Dockerfile` (Python 3.12.15 at this update).
- Python numerical/plotting/test versions: `requirements-repro.txt`.
- Node slide library: [PptxGenJS](https://github.com/gitbrent/PptxGenJS) 4.0.1,
  MIT, with full transitive npm integrity lockfile.
- [Quantum Basis](https://github.com/wztzjhn/quantum_basis), pinned commit
  `33324ce6081145ea9af981c64b2a2cc32c7faf0a`, GPL-3.0. CMake fetches the original
  source archive with SHA-256, including its upstream license. Not vendored.
- [ARPACK-NG](https://github.com/opencollab/arpack-ng), 3.9.1 (BSD-style license),
  source archive SHA-256 pinned. C++ enums needed by this qbasis revision are
  not present in all distro ARPACK packages, so those packages are not relied on.
- [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) (BSD-3-Clause) and LAPACKE:
  Debian system packages in the image; native builds use their LP64 installations.
  Boost, GCC/OpenMP and GFortran are normal distro build dependencies.

The Python/base image and qbasis/ARPACK/npm sources are pinned. Apt package
repositories are not snapshot-frozen. Installed package versions are recorded
at `/opt/mmhk/system-packages.tsv` in the full image. Rebuilds at later dates
may pick newer security/system packages and produce different binary hashes.
The physics tests, source/build identity gates and numerical tolerances cover
that distinction; we do not promise bitwise cross-platform reproducibility.

`cmake/patch_qbasis.py` makes the CPU banner portable beyond x86, changes the
ARPACK include path to its exported target, and preserves the historical
`lanczos_precision=2e-12` change (upstream was 1e-12). The two physics solver
sources are unchanged. No Hamiltonian or response operator is patched.

`reproducibility/compat/mkl.h` is a narrow compatibility layer for the exact
subset used here: CBLAS and LAPACKE from OpenBLAS, and checked zero-index CSR
non-transposed multiply (general or upper symmetric/Hermitian, nonunit
 diagonal). Hermitian mirrors conjugate off-diagonal values. Unsupported sparse
operations fail explicitly. FEAST methods are intentionally unavailable and
throw; these solvers use dense/IRAM/Lanczos, not FEAST. No Intel MKL license or
TBB installation is needed. Its name allows unchanged upstream includes; it
is not a general MKL implementation.

The pinned qbasis GPL terms apply to redistribution of its source and linked
binaries. Retain upstream licenses/notices if redistributing a built image.
This dependency notice does not assign a new license to the authors' own
repository code/data. No student code, unpublished Overleaf materials, or
book PDF is downloaded, required or redistributed by this pipeline.
