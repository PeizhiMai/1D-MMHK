#!/usr/bin/env python3
"""Relative-path archive integrity and final-output checks; launches no solver.

The file manifest deliberately excludes itself. Git authenticates that manifest
as part of the commit; its SHA-256 entries authenticate the remaining files.
Original campaign manifests/absolute provenance paths are not rewritten.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'repository_manifest.tsv'
CAMPAIGN = ROOT / 'response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def tracked_files():
    result = subprocess.run(['git', 'ls-files', '-z'], cwd=ROOT,
                            check=True, capture_output=True)
    return sorted(p.decode('utf-8') for p in result.stdout.split(b'\0')
                  if p and p.decode('utf-8') != MANIFEST.name)


def write_manifest():
    names = tracked_files()
    if not names:
        raise RuntimeError('Stage the intended archive files before making its manifest')
    with MANIFEST.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['sha256', 'bytes', 'path'], delimiter='\t')
        w.writeheader()
        for name in names:
            p = ROOT / name
            if not p.is_file() or p.is_symlink():
                raise RuntimeError(f'Archive requires regular files: {name}')
            w.writerow(dict(sha256=digest(p), bytes=p.stat().st_size, path=name))
    print(f'Wrote manifest for {len(names)} tracked files: {MANIFEST}')


def verify_manifest():
    with MANIFEST.open(newline='') as f:
        records = list(csv.DictReader(f, delimiter='\t'))
    names = [r['path'] for r in records]
    if len(set(names)) != len(names):
        raise RuntimeError('Duplicate manifest paths')
    for r in records:
        relative = Path(r['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise RuntimeError(f'Unsafe manifest path: {relative}')
        p = ROOT / relative
        if not p.is_file() or p.is_symlink():
            raise RuntimeError(f'Missing/unsupported archive file: {relative}')
        if p.stat().st_size != int(r['bytes']) or digest(p) != r['sha256']:
            raise RuntimeError(f'Archive hash/size mismatch: {relative}')
    # Git is optional for a downloaded source archive, but if present check
    # the index so additions/deletions cannot silently fall outside the manifest.
    if (ROOT / '.git').exists():
        if set(names) != set(tracked_files()):
            raise RuntimeError('Tracked-file set differs from repository_manifest.tsv')
    return dict(files=len(records), bytes=sum(int(r['bytes']) for r in records))


def verify_current_outputs():
    import numpy as np
    spectra = CAMPAIGN / 'single_particle_eta0p1_k101_depth600'
    config = json.loads((spectra / 'config.json').read_text())
    validation = json.loads((spectra / 'validation.json').read_text())
    if config['k_count'] != 101 or config['interpolate']:
        raise RuntimeError('Unexpected final momentum-grid contract')
    if not validation['passed'] or validation['plot_y_limits_omega_over_t'] != [-6, 10]:
        raise RuntimeError('Unexpected final spectral audit/display range')
    cases = 0
    for U in range(13):
        for tag in ['0p7', '0p8', '0p9']:
            p = spectra / f'U{U}' / f'n{tag}'
            if not (p / 'COMPLETE').is_file():
                raise RuntimeError(f'Incomplete final spectrum: {p}')
            with np.load(p / 'Akw.npz', allow_pickle=False) as d:
                np.testing.assert_array_equal(d['k'], np.linspace(0, np.pi, 101))
                np.testing.assert_allclose(d['omega'], np.linspace(-6, 16, 2201), atol=1e-12)
                np.testing.assert_allclose(d['eta'], .1, atol=1e-12)
                A = d['A_per_spin']
                if A.shape != (101, 2201) or not np.isfinite(A).all() or (A < 0).any():
                    raise RuntimeError(f'Invalid final spectrum: {p}')
            for name in ['spectral_Akw.png', 'spectral_Akw.svg']:
                if not (p / name).is_file():
                    raise RuntimeError(f'Missing figure: {p / name}')
            cases += 1
    analysis = CAMPAIGN / 'analysis_ba_step0p2_eta0p1'
    with (analysis / 'BA_velocities_fine_grid.csv').open(newline='') as f:
        ba = list(csv.DictReader(f))
    for n in [.7, .8, .9]:
        values = sorted(float(r['U_over_t']) for r in ba if abs(float(r['density']) - n) < 1e-12)
        np.testing.assert_allclose(values, np.linspace(0, 12, 61), atol=1e-12)
    for tag in ['0p7', '0p8', '0p9']:
        for p in [analysis / f'response_panel_n{tag}.png',
                  analysis / f'response_panel_n{tag}_full_range.png',
                  spectra / f'spectral_panel_n{tag}.png']:
            if not p.is_file():
                raise RuntimeError(f'Missing panel: {p}')
    deck = analysis / 'presentation/spin_charge_10MMHK_50twists_eta0p1_v2.pptx'
    import zipfile
    import re
    with zipfile.ZipFile(deck) as z:
        slides = [p for p in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml', p)]
        charts = [p for p in z.namelist() if re.fullmatch(r'ppt/slides/charts/chart\d+\.xml', p)]
        if len(slides) != 4 or len(charts) != 42:
            raise RuntimeError('Unexpected slide/native-chart count in final deck')
    return dict(spectral_cases=cases, momenta_per_case=101, stored_omega_range=[-6, 16],
                plotted_omega_range=[-6, 10], ba_grid_step=.2, ba_rows=len(ba),
                final_ppt_slides=len(slides), final_native_charts=len(charts))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write-manifest', action='store_true',
                   help='After staging archive files, generate its SHA-256 manifest')
    a = p.parse_args()
    if a.write_manifest:
        write_manifest()
        return
    result = dict(passed=True, archive=verify_manifest(), outputs=verify_current_outputs(),
                  many_body_solver_launched=False)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
