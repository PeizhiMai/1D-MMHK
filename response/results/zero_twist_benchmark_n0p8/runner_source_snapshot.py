#!/usr/bin/env python3
"""Run canonical zero-twist physical-density spectra at an exact density.

This uses the fast qbasis Hubbard-ring solver, never the student's solver.
Physical density is Ne/L, unlike the native solver's filling=Ne/(2L).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time

import numpy as np

from postprocess_response import reconstruct, save_result, plot_result

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tag(value):
    return f'{value:.8f}'.rstrip('0').rstrip('.').replace('.', 'p')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--density', type=float, default=.8)
    p.add_argument('--L', type=int, default=10)
    p.add_argument('--output', type=Path)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--depth', type=int, default=300)
    p.add_argument('--eta', type=float, default=.2)
    p.add_argument('--couplings', nargs='+', type=float, default=[4, 6, 7, 8, 12])
    args = p.parse_args()
    Ne = round(args.L*args.density)
    if (not 0 < args.density < 1 or abs(Ne-args.L*args.density) > 1e-9
            or args.L < 4 or args.L % 2 or args.threads < 1 or args.depth < 3
            or args.eta <= 0 or any(U < 0 for U in args.couplings)):
        p.error('Require even L, exact integer Ne=L*density, positive resources, U>=0')
    nup, ndn = (Ne+1)//2, Ne//2
    source = ROOT/'create_code/chain_Hubbard_tbc_response.cc'
    executable = HERE/'bin/chain_Hubbard_tbc_response.x'
    if not executable.is_file():
        p.error('Build the physical-density solver with response/build_response.sh')
    if executable.stat().st_mtime < source.stat().st_mtime:
        p.error('Solver source is newer than binary; rebuild before running')
    output = (args.output or HERE/f'results/zero_twist_benchmark_n{tag(args.density)}').resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = {
        'L': args.L, 'density_Ne_over_L': args.density, 'filling_Ne_over_2L': args.density/2,
        'Nup': nup, 'Ndn': ndn, 'phase_per_bond': 0, 'observable': 'density',
        'U_list': args.couplings, 'depth': args.depth, 'eta': args.eta,
        'OMP_NUM_THREADS': args.threads, 'MKL_NUM_THREADS': 1,
        'source_sha256': sha(source), 'executable_sha256': sha(executable),
        'runner_sha256': sha(__file__), 'postprocess_sha256': sha(HERE/'postprocess_response.py'),
        'host': platform.platform(), 'student_solver_executed': False, 'runs': [],
    }
    manifest_path = output/'manifest.json'
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if any(previous.get(k) != v for k, v in manifest.items() if k != 'runs'):
            raise RuntimeError('Existing output has different provenance/configuration; choose a new directory')
        manifest = previous
    else:
        shutil.copy2(source, output/'solver_source_snapshot.cc')
        shutil.copy2(executable, output/'solver_executable_snapshot.x')
        shutil.copy2(__file__, output/'runner_source_snapshot.py')
        manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    # Use an immutable run-specific executable; do not rebuild/overwrite shared tools.
    executable = output/'solver_executable_snapshot.x'
    if sha(executable) != manifest['executable_sha256']:
        raise RuntimeError('Executable snapshot hash mismatch')
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), MKL_NUM_THREADS='1',
               MPLCONFIGDIR=str(HERE/'.mplconfig'))
    os.environ['MPLCONFIGDIR'] = env['MPLCONFIGDIR']
    omega = np.linspace(-1, 16, 1701)
    for U in args.couplings:
        run = output/f'L{args.L}_U{U:g}_n{tag(args.density/2)}_phi0'
        run.mkdir(exist_ok=True)
        command = [str(executable), '--L', str(args.L), '--nup', str(nup), '--ndn', str(ndn),
                   '--U', str(U), '--phase', '0', '--observable', 'density',
                   '--depth', str(args.depth), '--output', str(run)]
        if not (run/'COMPLETE').exists():
            start = time.monotonic()
            print(f'RUN U/t={U:g}, L={args.L}, Ne={Ne}, physical n={args.density:g}', flush=True)
            with (run/'run.log').open('w') as log:
                subprocess.run(command, env=env, check=True, stdout=log, stderr=subprocess.STDOUT)
            manifest['runs'].append({'U': U, 'command': command,
                                     'elapsed_seconds': time.monotonic()-start})
            manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
        result = reconstruct(run, omega, args.eta)
        m = result['metadata']
        if (m['L'], m['Nup'], m['Ndn'], m['U'], m['phase_per_bond']) != (args.L, nup, ndn, U, 0):
            raise RuntimeError('Run metadata disagrees with requested system')
        summary = save_result(run, result, args.eta)
        plot_result(run, result, args.eta)
        print('COMPLETE', U, 'grounds:', m['ground_states'], 'checks:', summary['checks'], flush=True)
    print('Density benchmark output:', output, flush=True)


if __name__ == '__main__':
    main()
