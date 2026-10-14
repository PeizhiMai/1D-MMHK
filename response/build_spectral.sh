#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
Q=/Users/cosdis/Downloads/ED/ED_Zhentao/quantum_basis
I=/Users/cosdis/Downloads/ED/ED_Zhentao/installs
M=/opt/intel/oneapi/mkl/2023.2.2
T=/opt/intel/oneapi/tbb/2021.10.0
/usr/local/bin/g++-13 -std=c++17 -O3 -Wall -DMKL_ILP64 -fopenmp \
 -I"$Q/src" -I"$I/include" -I"$M/include" "$ROOT/create_code/chain_Hubbard_tbc_spectral.cc" \
 -L"$I/lib" -L"$M/lib" -L"$T/lib" -Wl,-rpath,"$I/lib" -Wl,-rpath,"$M/lib" -Wl,-rpath,"$T/lib" \
 -lqbasis -larpackILP64 -lmkl_intel_ilp64 -lmkl_tbb_thread -lmkl_core -ltbb -lpthread \
 -o "$ROOT/response/bin/chain_Hubbard_tbc_spectral.x"
install_name_tool -change libqbasis.so "$I/lib/libqbasis.so" "$ROOT/response/bin/chain_Hubbard_tbc_spectral.x"
