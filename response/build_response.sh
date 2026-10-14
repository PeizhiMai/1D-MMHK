#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QBASIS_ROOT="${QBASIS_ROOT:-/Users/cosdis/Downloads/ED/ED_Zhentao/quantum_basis}"
QBASIS_INSTALL="${QBASIS_INSTALL:-/Users/cosdis/Downloads/ED/ED_Zhentao/installs}"
MKL_ROOT="${MKL_ROOT:-/opt/intel/oneapi/mkl/2023.2.2}"
TBB_ROOT="${TBB_ROOT:-/opt/intel/oneapi/tbb/2021.10.0}"
CXX="${CXX:-/usr/local/bin/g++-13}"
mkdir -p "$ROOT/response/bin"
"$CXX" -std=c++17 -O3 -Wall -Wextra -DMKL_ILP64 -fopenmp \
  -I"$QBASIS_ROOT/src" -I"$QBASIS_INSTALL/include" -I"$MKL_ROOT/include" \
  "$ROOT/create_code/chain_Hubbard_tbc_response.cc" \
  -L"$QBASIS_INSTALL/lib" -L"$MKL_ROOT/lib" -L"$TBB_ROOT/lib" \
  -Wl,-rpath,"$QBASIS_INSTALL/lib" -Wl,-rpath,"$MKL_ROOT/lib" -Wl,-rpath,"$TBB_ROOT/lib" \
  -lqbasis -larpackILP64 -lmkl_intel_ilp64 -lmkl_tbb_thread -lmkl_core -ltbb -lpthread \
  -o "$ROOT/response/bin/chain_Hubbard_tbc_response.x"
if [[ "$(uname -s)" == Darwin ]]; then
  # This installed qbasis has the legacy bare install-name "libqbasis.so".
  install_name_tool -change libqbasis.so "$QBASIS_INSTALL/lib/libqbasis.so" \
    "$ROOT/response/bin/chain_Hubbard_tbc_response.x"
fi
echo "Built $ROOT/response/bin/chain_Hubbard_tbc_response.x"
