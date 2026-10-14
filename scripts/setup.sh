#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:---docker}"
if [[ "$MODE" == --docker || "$MODE" == --plots ]]; then
  command -v docker >/dev/null || { echo 'Install Docker Desktop (Mac/Windows) or Docker Engine (Linux) first.' >&2; exit 1; }
  if [[ "$MODE" == --plots ]]; then
    docker build --target plot-runtime -t mmhk-1d:plots "$ROOT"
  else
    docker build -t mmhk-1d:repro "$ROOT"
  fi
elif [[ "$MODE" == --native ]]; then
  PYTHON="${PYTHON:-python3}"
  "$PYTHON" -m venv "$ROOT/.repro-venv"
  "$ROOT/.repro-venv/bin/python" -m pip install -r "$ROOT/requirements-repro.txt"
  cmake -S "$ROOT" -B "$ROOT/.build/native" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$ROOT/.build/solver" -DPython3_EXECUTABLE="$ROOT/.repro-venv/bin/python"
  cmake --build "$ROOT/.build/native" -j "${BUILD_JOBS:-2}"
  cmake --install "$ROOT/.build/native"
  (cd "$ROOT/reproducibility/slides" && npm ci --ignore-scripts)
else
  echo 'Usage: bash scripts/setup.sh --docker | --plots | --native' >&2; exit 1
fi
