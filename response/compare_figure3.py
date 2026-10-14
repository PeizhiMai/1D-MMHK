#!/usr/bin/env python3
"""Compare zero-twist qbasis responses with downloaded Overleaf Figure 3."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from PIL import Image

import postprocess_response as post

HERE=Path(__file__).resolve().parent
REFERENCE=HERE/"benchmark_reference"


def reference_curve(U,channel):
    # Prefer actual arrays when available, verifying they match the Figure 3 label.
    candidates=list((REFERENCE/"source/N=10_interaction_sweep").glob(f"structure_factors_Norb10_U{U:.2f}_*_separated.npz"))
    if candidates:
        with np.load(candidates[0],allow_pickle=False) as d:
            w=d["omega_vals"].copy();s=d[f"S_{channel}_q_w"][1].copy()
            mask=(w>=0)&(w<=5)
            return w[mask],s[mask]/s[mask].max(),str(candidates[0]),float(d["x"])
    # Other U values: read the original PNGs in the LaTeX source. These have
    # identical 1800x1200 canvas and axes x=186..1742, spanning omega=0..5.
    # The zero y-axis tick is at row 980/981; only normalized shapes are compared.
    path=REFERENCE/f"source/velocity_figures/peak_verify_U={U:.2f}_{channel}.png"
    a=np.asarray(Image.open(path).convert("RGB"))
    if a.shape!=(1200,1800,3):raise ValueError(f"Unexpected reference image size: {path}")
    R,G,B=a[:,:,0],a[:,:,1],a[:,:,2]
    color=((R>180)&(G<80)&(B<160)) if channel=="charge" else ((B>90)&(R<80)&(G<80))
    color[:113]=False;color[1021:]=False
    color[120:295,1090:1730]=False  # legend, not data
    xs=np.arange(188,1741);ys=np.full(len(xs),np.nan)
    for i,x in enumerate(xs):
        y=np.flatnonzero(color[:,x]);
        if len(y):ys[i]=np.median(y)
    good=np.isfinite(ys)
    if good.sum()<0.9*len(xs):raise ValueError(f"Insufficient digitized curve pixels: {path}")
    ys=np.interp(xs,xs[good],ys[good])
    s=np.maximum(0,980.5-ys);s/=s.max()
    w=(xs-186)/((1742-186)/5)
    output=REFERENCE/"digitized";output.mkdir(exist_ok=True)
    np.savetxt(output/f"U{U:g}_{channel}.csv",np.column_stack([w,s]),delimiter=",",header="omega,normalized_intensity",comments="")
    return w,s,str(path),None


def compare(directories,output,eta=0.2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    targets=json.loads((REFERENCE/"figure3_targets.json").read_text())
    omega=np.linspace(0,5,5001)
    records=[];output=Path(output);output.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(len(directories),2,figsize=(9,2.1*len(directories)),layout="constrained",squeeze=False)
    for i,directory in enumerate(directories):
        result=post.reconstruct(directory,omega,eta);U=result["metadata"]["U"]
        for j,channel in enumerate(("charge","spin")):
            wref,sref,source,actual_filling=reference_curve(U,channel)
            s=result[channel][1];peak=post.first_peak(omega,s)
            labelled=targets["peaks"][f"{U:g}"][channel]
            normalized=s/s.max();rmse=float(np.sqrt(np.mean((np.interp(omega,wref,sref)-normalized)**2)))
            delta=None if peak is None else peak-labelled
            record={"U":U,"channel":channel,"filling_qbasis":result["metadata"]["filling"],
                    "filling_reference":actual_filling,"reference_labelled_peak":labelled,
                    "qbasis_first_peak":peak,"qbasis_dominant_peak":float(omega[np.argmax(s)]),
                    "peak_difference":delta,"normalized_shape_rmse":rmse,
                    "within_two_reference_grid_cells":bool(delta is not None and abs(delta)<=2*targets["reference_grid_spacing"]),
                    "reference_curve_source":source}
            records.append(record)
            ax=axes[i,j]
            ax.plot(wref,sref,color="crimson",label="Figure 3 reference",lw=1.3)
            ax.plot(omega,normalized,color="navy",label="qbasis: coherent density",lw=1.1)
            ax.axvline(labelled,color="crimson",ls=":",lw=.9)
            if peak is not None:ax.axvline(peak,color="navy",ls="--",lw=.9)
            ax.set(xlim=(0,5),ylim=(0,1.07),title=f"U/t={U:g}, {channel}",ylabel="Normalized intensity")
            if i==0:ax.legend(fontsize=8,loc="upper right")
            if i==len(directories)-1:ax.set_xlabel(r"$\omega/t$")
    fig.suptitle(r"Zero twist, 10 orbitals, $n=0.45$, $q=2\pi/10$, $\eta/t=0.2$"+"\nIndependent peak/shape comparison; amplitude conventions not matched",fontsize=11)
    fig.savefig(output/"figure3_comparison.png",dpi=170);plt.close(fig)
    (output/"figure3_comparison.json").write_text(json.dumps(records,indent=2)+"\n")
    with (output/"figure3_comparison.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    return records


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("directories",nargs="+",type=Path)
    p.add_argument("--output",type=Path,default=HERE/"results/zero_twist_benchmark")
    p.add_argument("--eta",type=float,default=.2)
    args=p.parse_args()
    for row in compare(args.directories,args.output,args.eta):print(row)
