#!/usr/bin/env python3
"""Run the five requested couplings locally, postprocess and compare Figure 3."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import numpy as np

from postprocess_response import reconstruct,save_result,plot_result
from compare_figure3 import compare

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=HERE/"results/zero_twist_benchmark")
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--depth",type=int,default=300)
    p.add_argument("--eta",type=float,default=.2)
    p.add_argument("--couplings",nargs="+",type=float,default=[4,6,7,8,12])
    args=p.parse_args()
    if args.threads<1 or args.depth<3 or args.eta<=0:p.error("Invalid resources or broadening")
    exe=HERE/"bin/chain_Hubbard_tbc_response.x"
    if not exe.is_file():p.error("Build first with bash response/build_response.sh")
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    manifest={"L":10,"filling":.45,"Nup":5,"Ndn":4,"phase_per_bond":0,
              "observable":"density","U_list":args.couplings,"depth":args.depth,"eta":args.eta,
              "OMP_NUM_THREADS":args.threads,"MKL_NUM_THREADS":1,
              "source_sha256":sha(ROOT/"create_code/chain_Hubbard_tbc_response.cc"),
              "executable_sha256":sha(exe),"host":platform.platform(),"runs":[]}
    old=output/"manifest.json"
    if old.exists():
        previous=json.loads(old.read_text())
        if any(previous.get(k)!=manifest[k] for k in manifest if k!="runs"):
            raise RuntimeError("Existing benchmark has different provenance/configuration; choose another output directory")
        manifest=previous
    else:
        (output/"solver_source_snapshot.cc").write_bytes((ROOT/"create_code/chain_Hubbard_tbc_response.cc").read_bytes())
    env=dict(os.environ,OMP_NUM_THREADS=str(args.threads),MKL_NUM_THREADS="1",
             MPLCONFIGDIR=str(HERE/".mplconfig"))
    os.environ["MPLCONFIGDIR"]=env["MPLCONFIGDIR"]
    directories=[];omega=np.linspace(-1,16,1701)
    for U in args.couplings:
        run=output/f"L10_U{U:g}_n0p45_phi0";run.mkdir(exist_ok=True)
        command=[str(exe),"--L","10","--filling","0.45","--U",str(U),"--phase","0",
                 "--depth",str(args.depth),"--output",str(run)]
        if not (run/"COMPLETE").exists():
            start=time.monotonic()
            with (run/"run.log").open("w") as log:
                subprocess.run(command,env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
            manifest["runs"].append({"U":U,"command":command,"elapsed_seconds":time.monotonic()-start})
            old.write_text(json.dumps(manifest,indent=2)+"\n")
        result=reconstruct(run,omega,args.eta);summary=save_result(run,result,args.eta);plot_result(run,result,args.eta)
        print("COMPLETE",U,summary["checks"],flush=True);directories.append(run)
    records=compare(directories,output,args.eta)
    for r in records:print(f"U={r['U']:g} {r['channel']}: qbasis={r['qbasis_first_peak']:.3f}, Figure3={r['reference_labelled_peak']:.3f}, delta={r['peak_difference']:+.3f}")
    print("Benchmark output:",output)


if __name__=="__main__":main()
