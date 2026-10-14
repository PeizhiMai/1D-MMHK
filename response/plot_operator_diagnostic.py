#!/usr/bin/env python3
"""Plot the saved U=8 coherent/pair response comparison and q=0 diagnostic."""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from postprocess_response import reconstruct


def main():
    here=Path(__file__).resolve().parent;base=here/"results/zero_twist_benchmark"
    omega=np.linspace(0,5,5001)
    physical=reconstruct(base/"L10_U8_n0p45_phi0",omega,.2)
    pairs=reconstruct(here/"results/diagnostic_U8_pairs",omega,.2)
    pairs0=reconstruct(here/"results/diagnostic_U8_pairs_q0",omega,.2)
    ref=np.load(here/"benchmark_reference/source/N=10_interaction_sweep/structure_factors_Norb10_U8.00_x0.453_separated.npz",allow_pickle=False)
    fig,axes=plt.subplots(2,2,figsize=(8,5.5),layout="constrained")
    for j,channel in enumerate(["charge","spin"]):
        ax=axes[0,j]
        for data,label,color,style in [(physical,"Coherent density","navy","-"),(pairs,"Orbital-pair sum","darkorange","--")]:
            s=data[channel][1];ax.plot(omega,s/s.max(),label=label,color=color,ls=style,lw=1.2)
        w=ref["omega_vals"];s=ref[f"S_{channel}_q_w"][1];valid=(w>=0)&(w<=5)
        ax.plot(w[valid],s[valid]/s[valid].max(),color="crimson",label="Figure 3",lw=1.2)
        ax.set(title=f"{channel.capitalize()}, q=2pi/10",ylabel="Normalized intensity",xlabel=r"$\omega/t$",xlim=(0,5),ylim=(0,1.07))
        if j==0:ax.legend(fontsize=8)
        ax=axes[1,j]
        ax.plot(omega,pairs0[channel][0],color="darkorange",label="Orbital-pair sum")
        ax.plot(omega,physical[channel][0],color="navy",label="Conserved total density")
        ax.set(title=f"{channel.capitalize()}, q=0",ylabel=r"$S(q=0,\omega)$",xlabel=r"$\omega/t$",xlim=(0,5))
        if j==0:ax.legend(fontsize=8)
    ref.close()
    fig.suptitle("Operator-definition diagnostic: L=10, U=8, n=0.45, zero twist",fontsize=11)
    fig.savefig(base/"operator_diagnostic_U8.png",dpi=170);plt.close(fig)


if __name__=="__main__":main()
