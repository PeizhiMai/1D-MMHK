#!/usr/bin/env python3
"""Portable entry point. Calculate into a new run or replot archived data."""
import argparse
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
ARCHIVE=ROOT/'response/results/twist50_common_mu_U0to12_n0p7_0p8_0p9'
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ.setdefault(key,'1')


def slide_export(output):
    node=shutil.which('node')
    if not node:raise RuntimeError('Node.js is required for slide export; use the Docker image or install the documented public slide dependencies')
    subprocess.run([node,str(ROOT/'reproducibility/slides/export.mjs'),str(output/'presentation_data.json'),str(output)],check=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    for name in ['plan','calculate']:
        p=commands.add_parser(name)
        p.add_argument('--config',type=Path,default=ROOT/'configs/study.json')
        p.add_argument('--output',type=Path,default=ROOT/'runs/study')
        p.add_argument('--solver-dir',type=Path)
        if name=='calculate':
            p.add_argument('--stage',choices=['all','response','spectral'],default='all')
            p.add_argument('--no-slides',action='store_true')
    p=commands.add_parser('plot');p.add_argument('--input',type=Path,default=ARCHIVE)
    p.add_argument('--output',type=Path,default=ROOT/'runs/archive-plots');p.add_argument('--config',type=Path)
    p.add_argument('--eta',type=float);p.add_argument('--no-slides',action='store_true')
    p=commands.add_parser('doctor');p.add_argument('--solver-dir',type=Path)
    args=parser.parse_args()
    if args.command=='doctor':
        for mod in ['numpy','scipy','matplotlib']:print(mod,importlib.import_module(mod).__version__)
        from reproducibility.workflow import solvers,sha
        for k,p in solvers(args.solver_dir).items():
            subprocess.run([str(p),'--help'],check=True,stdout=subprocess.DEVNULL);print(k,p,sha(p))
        return
    from reproducibility.workflow import load_config
    if args.command in ['plan','calculate']:
        c=load_config(args.config)
        if args.command=='plan':print(json.dumps(c,indent=2));print('Plan only: no solver launched.');return
        from reproducibility.workflow import solvers,initialize,responses,spectra
        exe=solvers(args.solver_dir);out=initialize(c,args.output,exe)
        if args.stage in ['all','response']:responses(c,out,exe)
        if args.stage in ['all','spectral']:spectra(c,out,exe)
        if args.stage=='all':
            from reproducibility.plots import plot_all
            report=plot_all(out,out/'plots',c)
            if not args.no_slides:slide_export(out/'plots')
            print(json.dumps(report,indent=2))
        return
    root=args.input.resolve()
    cfg=args.config or (root/'config.json' if 'twists' in json.loads((root/'config.json').read_text()) else ROOT/'configs/study.json')
    c=load_config(cfg)
    if args.eta is not None:c['eta']=args.eta
    from reproducibility.plots import plot_all
    report=plot_all(root,args.output.resolve(),c)
    if not args.no_slides:slide_export(args.output.resolve())
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
