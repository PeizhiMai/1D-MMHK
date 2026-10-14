#!/usr/bin/env bash
# Run from a fresh clone. Only this checkout is mounted into the containers.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
OUT="runs/clean-test-$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$OUT"
CACHE=(); [[ "${MMHK_NO_CACHE:-0}" == 1 ]] && CACHE=(--no-cache)
docker build "${CACHE[@]}" -t mmhk-1d:repro . 2>&1 | tee "$OUT/build-full.log"
docker build "${CACHE[@]}" --target plot-runtime -t mmhk-1d:plots . 2>&1 | tee "$OUT/build-plots.log"
bash scripts/run_docker.sh doctor | tee "$OUT/doctor.txt"
docker run --rm -v "$ROOT:/workspace" --entrypoint python mmhk-1d:repro \
  -m pytest reproducibility/tests -q 2>&1 | tee "$OUT/pytest.txt"
bash scripts/run_docker.sh calculate --config configs/smoke.json --output "$OUT/smoke" \
  2>&1 | tee "$OUT/smoke.log"
# Check plot-only without even a solver installation.
docker run --rm --entrypoint sh mmhk-1d:plots -c \
  'test ! -e /opt/mmhk/solver && ! command -v chain_Hubbard_tbc_response'
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plot \
  --input "$OUT/smoke" --output "$OUT/smoke-replot" 2>&1 | tee "$OUT/replot.log"
MMHK_IMAGE=mmhk-1d:plots bash scripts/run_docker.sh plot \
  --output "$OUT/archive-plots" 2>&1 | tee "$OUT/archive-plots.log"
# git is unnecessary: the relative manifest authenticates a source tar/clone.
docker run --rm -v "$ROOT:/workspace" --entrypoint python mmhk-1d:plots \
  response/verify_repository_archive.py | tee "$OUT/archive-verification.json"
docker run --rm -v "$ROOT:/workspace" --entrypoint python mmhk-1d:plots \
  reproducibility/check_outputs.py "$OUT" | tee "$OUT/clean-test.json"
echo "PASS clean test: $ROOT/$OUT"
