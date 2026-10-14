#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# No home directories, old executables, or installed Codex skills are mounted.
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp/mmhk-matplotlib \
  -v "$ROOT:/workspace" -w /workspace "${MMHK_IMAGE:-mmhk-1d:repro}" "$@"
