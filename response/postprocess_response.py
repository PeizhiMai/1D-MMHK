#!/usr/bin/env python3
"""Reconstruct physical density spectra, poles and checks from qbasis fractions."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.signal import find_peaks


def read_fraction(path):
    data = {}
    for line in Path(path).read_text().splitlines():
        key, *values = line.split()
        data[key] = np.asarray(values, dtype=float) if key in ("a", "b") else float(values[0])
    return data


def poles_from_fraction(data):
    if data["steps"] == 0:
        return np.array([]), np.array([])
    a, b = data["a"], data["b"]
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Nonfinite Lanczos coefficients")
    energies, vectors = eigh_tridiagonal(a, b[1:], lapack_driver="stev")
    return energies - data["E0"], data["norm"]**2 * vectors[0]**2


def lorentzian(omega, energies, weights, eta):
    return np.sum(weights[:, None] * eta / np.pi /
                  ((omega[None, :] - energies[:, None])**2 + eta**2), axis=0)


def first_peak(omega, spectrum, cutoff=5.0, prominence=0.0):
    valid = (omega > 0) & (omega <= cutoff)
    w, s = omega[valid], spectrum[valid]
    if len(w) < 3 or np.max(s) <= 1e-12:
        return None
    indices, _ = find_peaks(s, height=np.max(s)*1e-3, prominence=prominence*np.max(s))
    if not len(indices):
        return None  # Do not quietly replace "no local peak" by a global maximum.
    i = indices[0]
    return float(w[i])


def reconstruct(directory, omega, eta, elastic_tol=1e-8):
    directory = Path(directory)
    if not (directory / "COMPLETE").is_file():
        raise ValueError(f"Incomplete calculation: {directory}")
    metadata = json.loads((directory / "metadata.json").read_text())
    if metadata.get("ground_only",False):raise ValueError("Ground-only run has no spectra")
    L = metadata["L"]
    result = {"omega": omega, "q": 2*np.pi*np.arange(L)/L,
              "charge": np.zeros((L, len(omega))), "spin": np.zeros((L, len(omega))),
              "elastic_charge": np.zeros(L), "elastic_spin": np.zeros(L),
              "weight_charge": np.zeros(L), "weight_spin": np.zeros(L)}
    rows, diagnostics = [], []
    seen={"charge":set(),"spin":set()}
    components={"charge":set(),"spin":set()}
    for channel in ("charge", "spin"):
        for path in sorted(directory.glob(f"{channel}_gsQ*_q*.dat")):
            d = read_fraction(path)
            q, factor = int(d["q_index"]), d["ensemble_weight"]
            seen[channel].add(q)
            component=(int(d["ground_Q"]),int(d.get("ground_id",0)),q,int(d.get("source_orbital",-1)))
            if component in components[channel]:raise ValueError(f"Duplicate source {component}")
            components[channel].add(component)
            if abs(factor-1/len(metadata["ground_states"]))>1e-12:
                raise ValueError("Incorrect ensemble weight")
            energies, weights = poles_from_fraction(d)
            if np.any((energies < -elastic_tol) & (weights > 1e-8)):
                raise ValueError(f"Significant negative excitation energy in {path}")
            elastic = np.abs(energies) <= elastic_tol
            positive = energies > elastic_tol
            result[channel][q] += factor * lorentzian(omega, energies[positive], weights[positive], eta)
            result[f"elastic_{channel}"][q] += factor*(d["elastic_weight"] + weights[elastic].sum())
            result[f"weight_{channel}"][q] += factor*(d["elastic_weight"] + weights.sum())
            for energy, weight in zip(energies, weights):
                rows.append({"channel": channel, "q_index": q, "ground_Q": int(d["ground_Q"]),
                             "ground_id":int(d.get("ground_id",0)),
                             "source_orbital":int(d.get("source_orbital",-1)),
                             "omega": float(energy), "weight": float(weight*factor),
                             "elastic": bool(abs(energy) <= elastic_tol)})
            if "source_weight" in d:
                diagnostics.append({"channel": channel, "q_index": q, "ground_Q": int(d["ground_Q"]),
                                    "weight_error": float(abs(weights.sum()-d["source_weight"])),
                                    "first_moment_error": float(abs(energies@weights-d["first_moment"]))})
    result["metadata"] = metadata
    if seen["charge"]!=seen["spin"] or not seen["charge"]:
        raise ValueError("Missing charge/spin response files")
    expected=set(metadata.get("q_indices",seen["charge"]))
    if seen["charge"]!=expected:
        raise ValueError("Incomplete set of q indices")
    orbital_indices=range(L) if metadata.get("observable","density")=="orbital-pairs" else [-1]
    required={(g["Q"],g.get("id",0),q,orbital) for g in metadata["ground_states"] for q in expected for orbital in orbital_indices}
    if any(components[channel]!=required for channel in ("charge","spin")):
        raise ValueError("Missing or unexpected ground-state/orbital response components")
    result["computed_q_indices"]=np.array(sorted(expected),dtype=int)
    for channel in ("charge","spin"):
        for q in set(range(L))-expected:
            result[channel][q]=np.nan
            result[f"weight_{channel}"][q]=np.nan
            result[f"elastic_{channel}"][q]=np.nan
    result["poles"] = rows
    result["moment_checks"] = diagnostics
    return result


def save_result(directory, result, eta):
    directory = Path(directory)
    np.savez_compressed(directory / "spectra.npz", **{k: v for k,v in result.items() if isinstance(v,np.ndarray)}, eta=eta)
    with (directory / "poles.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["channel","q_index","ground_Q","ground_id","source_orbital","omega","weight","elastic"])
        writer.writeheader();writer.writerows(result["poles"])
    summary={"metadata":result["metadata"],"eta":eta,"moment_checks":result["moment_checks"],"modes":[]}
    for q in result["computed_q_indices"]:
        q=int(q)
        row={"q_index":q,"q":float(result["q"][q])}
        for channel in ("charge","spin"):
            total=result[f"weight_{channel}"][q]
            significant=[p["omega"] for p in result["poles"] if p["q_index"]==q
                         and p["channel"]==channel and not p["elastic"] and p["omega"]>0
                         and p["weight"]>max(1e-10,1e-4*total)]
            row[f"{channel}_first_peak"]=first_peak(result["omega"],result[channel][q])
            row[f"{channel}_lowest_pole"]=float(min(significant)) if significant else None
            row[f"{channel}_total_weight"]=float(total)
            row[f"{channel}_elastic_weight"]=float(result[f"elastic_{channel}"][q])
        summary["modes"].append(row)
    checks={
        "q0_charge_inelastic_max":float(np.max(result["charge"][0])) if 0 in result["computed_q_indices"] else None,
        "q0_spin_inelastic_max":float(np.max(result["spin"][0])) if 0 in result["computed_q_indices"] else None,
        "max_weight_error":max((x["weight_error"] for x in result["moment_checks"]),default=0),
        "max_first_moment_error":max((x["first_moment_error"] for x in result["moment_checks"]),default=0),
    }
    if len(result["metadata"]["ground_states"])>1 or result["metadata"]["phase_per_bond"]==0:
        pairs=[int(q) for q in result["computed_q_indices"] if (-q)%len(result["q"]) in result["computed_q_indices"]]
        for channel in ("charge","spin"):
            checks[f"{channel}_q_minus_q_max_difference"]=float(max((np.max(np.abs(result[channel][q]-result[channel][(-q)%len(result["q"])])) for q in pairs),default=0)) if pairs else None
    summary["checks"]=checks
    (directory / "summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    return summary


def plot_result(directory,result,eta):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    metadata=result["metadata"];L=metadata["L"]
    fig,axes=plt.subplots(1,2,figsize=(8,3.2),layout="constrained")
    for ax,channel in zip(axes,("charge","spin")):
        w=result["omega"]/metadata["t"];dq=2*np.pi/L;dw=w[1]-w[0]
        im=ax.imshow(result[channel][:L//2+1].T,origin="lower",aspect="auto",cmap="turbo",interpolation="none",
                     extent=[-dq/2,np.pi+dq/2,w[0]-dw/2,w[-1]+dw/2])
        ax.set(xlabel=r"$q$",ylabel=r"$\omega/t$",title=channel.capitalize(),xticks=[0,np.pi/2,np.pi],xticklabels=["0",r"$\pi/2$",r"$\pi$"])
        ax.set_xlim(0,np.pi);ax.set_ylim(w[0],w[-1])
        fig.colorbar(im,ax=ax,label=r"$S(q,\omega)$")
    fig.suptitle(f"L={L}, U/t={metadata['U']/metadata['t']:g}, n={metadata['filling']:g}, phase={metadata['phase_per_bond']:g}, eta={eta:g}")
    if set(range(L//2+1)).issubset(result["computed_q_indices"]):
        fig.savefig(Path(directory)/"response_map.png",dpi=160)
    plt.close(fig)
    if 1 not in result["computed_q_indices"]:return
    fig,ax=plt.subplots(figsize=(5,3.2),layout="constrained")
    for channel,color in [("charge","crimson"),("spin","navy")]:
        ax.plot(result["omega"]/metadata["t"],result[channel][1],color=color,label=channel)
    ax.set(xlabel=r"$\omega/t$",ylabel=r"$S(q_{\min},\omega)$",xlim=(0,5),title=rf"$q_{{\min}}=2\pi/{L}$, $U/t={metadata['U']/metadata['t']:g}$")
    ax.legend();fig.savefig(Path(directory)/"qmin_response.png",dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directories",nargs="+",type=Path)
    p.add_argument("--eta",type=float,default=0.2)
    p.add_argument("--omega-min",type=float,default=-1)
    p.add_argument("--omega-max",type=float,default=16)
    p.add_argument("--points",type=int,default=1701)
    p.add_argument("--no-plots",action="store_true")
    args=p.parse_args()
    if args.eta<=0 or args.points<3 or args.omega_max<=args.omega_min:p.error("Invalid frequency grid or broadening")
    omega=np.linspace(args.omega_min,args.omega_max,args.points)
    for d in args.directories:
        result=reconstruct(d,omega,args.eta);summary=save_result(d,result,args.eta)
        if not args.no_plots:plot_result(d,result,args.eta)
        print(d,summary["checks"],next((m for m in summary["modes"] if m["q_index"]==1),None))


if __name__=="__main__":main()
