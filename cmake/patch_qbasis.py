"""Minimal portability patch to the pinned external qbasis source tree."""
import sys
from pathlib import Path
p=Path(sys.argv[1])/'src'
f=p/'miscellaneous.cc';s=f.read_text()
if '// MMHK portable CPU banner' not in s:
    s=s.replace('#include <cpuid.h>', '#if defined(__i386__) || defined(__x86_64__)\n#include <cpuid.h>\n#endif',1)
    s=s.replace('        char CPUBrandString[0x40];','        // MMHK portable CPU banner\n#if defined(__i386__) || defined(__x86_64__)\n        char CPUBrandString[0x40];',1)
    s=s.replace('        std::cout << "CPU Type:', '#else\n        const char* CPUBrandString="non-x86 (OpenBLAS)";\n#endif\n        std::cout << "CPU Type:',1)
    # Preserve the original study's ONLY numerical qbasis source change.
    assert 'const double lanczos_precision = 1e-12;' in s
    s=s.replace('const double lanczos_precision = 1e-12;', 'const double lanczos_precision = 2e-12;',1)
    f.write_text(s)
f=p/'lanczos.cc';s=f.read_text().replace('"arpack-ng/arpack.hpp"','<arpack.hpp>');f.write_text(s)
