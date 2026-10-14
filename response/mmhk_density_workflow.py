#!/usr/bin/env python3
"""Physical MMHK density response = common-mu, reduced-BZ sum of TBC rings.

This module does not import or execute any student module. Raw calculations use
the qbasis C++ solver. All neutral spectra use rho(q)=sum_j e^(iqj)n_j/sqrt(L).
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from postprocess_response import reconstruct, save_result, first_peak

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DEGENERACY_TOL=1e-8


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def k_grid(L,nk,kind):
    half=np.pi/L
    if kind=="student":
        k=np.linspace(-half,half,nk,endpoint=True)
    elif kind=="periodic":
        # Distinct momenta of a zero-external-twist L*nk-site periodic system.
        # Includes k=0, counts the identified RBZ endpoints only once.
        k=2*half*(np.arange(nk)-nk//2)/nk
    elif kind=="midpoint":
        k=-half+(np.arange(nk)+.5)*2*half/nk
    else:raise ValueError("Unknown grid kind")
    return k,np.ones(nk)/nk


def number_minimizers(energies,mu,tol=DEGENERACY_TOL):
    energies=np.asarray(energies)
    omega=energies-mu*np.arange(energies.shape[-1])
    return [np.flatnonzero(row-row.min()<tol).tolist() for row in omega]


def choose_mu(energies,weights,target,L):
    """Choose a stable nearest-filling PLATEAU, not a fake continuous root.

    Obtain all lower-convex-envelope transition chemical potentials; cluster
    numerical copies, then examine finite intervals between true transitions.
    If the target lies in a density jump, report the closest attainable plateau.
    No implicit mixture of particle-number sectors is used to manufacture it.
    """
    energies=np.asarray(energies);thresholds=[]
    for row in energies:
        hull=[]
        for N,E in enumerate(row):
            while len(hull)>=2:
                a,b=hull[-2:]
                slope1=(row[b]-row[a])/(b-a)
                slope2=(E-row[b])/(N-b)
                if slope2<=slope1+1e-10:hull.pop()
                else:break
            hull.append(N)
        thresholds.extend((row[b]-row[a])/(b-a) for a,b in zip(hull[:-1],hull[1:]))
    groups=[]
    for value in sorted(thresholds):
        if not groups or value-groups[-1][-1]>1e-8:groups.append([float(value)])
        else:groups[-1].append(float(value))
    thresholds=np.array([np.mean(group) for group in groups])
    edges=np.r_[thresholds[0]-10,thresholds,thresholds[-1]+10]
    plateaus=[]
    for left,right in zip(edges[:-1],edges[1:]):
        if right-left<1e-7:continue
        mu=float((left+right)/2)
        minima=number_minimizers(energies,mu)
        if any(len(row)!=1 for row in minima):continue
        numbers=np.array([row[0] for row in minima])
        filling=float(weights@numbers/(2*L))
        plateaus.append({"mu":mu,"mu_left":float(left),"mu_right":float(right),
                         "filling":filling,"numbers":numbers.tolist()})
    if not plateaus:raise RuntimeError("No stable grand-canonical plateau found")
    chosen=min(plateaus,key=lambda row:(round(abs(row["filling"]-target),12),abs(row["mu"])))
    return dict(chosen, target=float(target),filling_error=chosen["filling"]-target,
                exact_target_attained=abs(chosen["filling"]-target)<1e-10,
                policy="nearest attainable stable zero-temperature plateau",plateaus=plateaus)


def finite_ring_probe_weights(q,L):
    """Windowed finite-ring probe, NOT the generic-q physical MMHK response.

    For a translation-invariant ground ensemble its inelastic spectrum is
    sum_m |a_m|^2 S(q_m), a_m=L^-1 sum_j exp[i(q-q_m)j].
    """
    q_grid=2*np.pi*np.arange(L)/L
    a=np.mean(np.exp(1j*(float(q)-q_grid[:,None])*np.arange(L)),axis=1)
    weights=abs(a)**2
    np.testing.assert_allclose(weights.sum(),1.,atol=1e-12)
    return weights


class QbasisRuns:
    def __init__(self,root,L,U,depth,threads,solver_sha,source_sha):
        self.root=Path(root);self.L=L;self.U=U;self.depth=depth
        self.exe=HERE/"bin/chain_Hubbard_tbc_response.x"
        self.env=dict(os.environ,OMP_NUM_THREADS=str(threads),MKL_NUM_THREADS="1")
        self.solver_sha=solver_sha;self.source_sha=source_sha;self.new_runs=0;self.elapsed=0.

    def run(self,k_index,k,Nup,Ndn,ground_only):
        kind="energy_scan" if ground_only else "canonical_responses"
        output=self.root/kind/f"k{k_index:03d}"/f"Nup{Nup}_Ndn{Ndn}"
        command=[str(self.exe),"--L",str(self.L),"--U",str(self.U),"--nup",str(Nup),
                 "--ndn",str(Ndn),"--phase",repr(float(-k)),"--depth",str(self.depth),
                 "--observable","density","--output",str(output.resolve())]
        if ground_only:command.append("--ground-only")
        identity={"command":command,"source_sha256":self.source_sha,"executable_sha256":self.solver_sha}
        lineage=output/"lineage.json"
        if (output/"COMPLETE").exists():
            if not lineage.exists() or json.loads(lineage.read_text())!=identity:
                raise RuntimeError(f"Incompatible cached calculation: {output}")
        else:
            if (output/"metadata.json").exists():
                raise RuntimeError(f"Interrupted run needs a fresh output directory: {output}")
            output.mkdir(parents=True,exist_ok=True);lineage.write_text(json.dumps(identity,indent=2)+"\n")
            start=time.monotonic()
            with (output/"run.log").open("w") as log:
                subprocess.run(command,env=self.env,check=True,stdout=log,stderr=subprocess.STDOUT)
            self.elapsed+=time.monotonic()-start;self.new_runs+=1
        metadata=json.loads((output/"metadata.json").read_text())
        if metadata["L"]!=self.L or metadata["U"]!=self.U or metadata["Nup"]!=Nup or metadata["Ndn"]!=Ndn:
            raise RuntimeError("Canonical run metadata mismatch")
        if abs(metadata["phase_per_bond"]+k)>1e-12:raise RuntimeError("Wrong block/twist gauge")
        return output,metadata


def calculate_coupling(root,L,U,k_values,weights,omega,eta,target,mu_override,depth,threads,number_policy,solver_sha,source_sha):
    runs=QbasisRuns(root,L,U,depth,threads,solver_sha,source_sha)
    energies=np.empty((len(k_values),2*L+1))
    for ik,k in enumerate(k_values):
        for N in range(2*L+1):
            _,m=runs.run(ik,k,(N+1)//2,N//2,True)
            energies[ik,N]=m["E0"]
        print(f"U={U:g} ground table k={ik+1}/{len(k_values)}",flush=True)
    if mu_override is None:
        selection=choose_mu(energies,weights,target,L);mu=selection["mu"]
    else:
        mu=float(mu_override);selection={"mu":mu,"target":target,"policy":"explicit chemical potential"}
    chosen_numbers=number_minimizers(energies,mu)
    if number_policy=="reject" and any(len(ns)>1 for ns in chosen_numbers):
        raise RuntimeError("mu lies at a particle-number crossing; choose an interior mu or --number-degeneracy uniform")
    output=Path(root);output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/"ground_energy_table.npz",k=k_values,weights=weights,E0_by_N=energies)
    per_k={ch:np.zeros((len(k_values),L,len(omega))) for ch in ["charge","spin"]}
    total_weight={ch:np.zeros((len(k_values),L)) for ch in ["charge","spin"]}
    elastic={ch:np.zeros((len(k_values),L)) for ch in ["charge","spin"]}
    first_moment={ch:np.zeros((len(k_values),L)) for ch in ["charge","spin"]}
    k_report=[];counts=[];poles=[]
    max_weight_error=0.;max_first_error=0.;max_residual=0.
    for ik,(k,wk,numbers) in enumerate(zip(k_values,weights,chosen_numbers)):
        active=[]
        for N in numbers:
            E_min=energies[ik,N]
            for Nu in range((N+1)//2,min(L,N)+1):
                Nd=N-Nu
                if Nd<0 or Nd>L:continue
                _,m=runs.run(ik,k,Nu,Nd,True)
                if m["E0"]>E_min+DEGENERACY_TOL:break  # SU(2): restricted-Sz energies are monotone
                if abs(m["E0"]-E_min)>DEGENERACY_TOL:raise RuntimeError("Minimum-Sz scan missed a lower energy")
                path,m=runs.run(ik,k,Nu,Nd,False)
                if not m.get("ground_subspace_complete"):raise RuntimeError("Uncertified ground subspace")
                data=reconstruct(path,omega,eta);save_result(path,data,eta)
                active.append({"N":N,"Nup":Nu,"Ndn":Nd,"spin_flip_multiplicity":1 if Nu==Nd else 2,
                               "ground_dimension":len(m["ground_states"]),"directory":str(path.resolve()),"data":data})
        dimension=sum(row["spin_flip_multiplicity"]*row["ground_dimension"] for row in active)
        if not dimension:raise RuntimeError("Empty grand-canonical ground subspace")
        mean_N=0.;mean_N_squared=0.;mean_spin_squared=0.;components=[]
        for row in active:
            probability=row["spin_flip_multiplicity"]*row["ground_dimension"]/dimension
            data=row["data"];mean_N+=probability*row["N"];mean_N_squared+=probability*row["N"]**2
            mean_spin_squared+=probability*(row["Nup"]-row["Ndn"])**2
            components.append({key:value for key,value in row.items() if key!="data"}|{"probability":probability})
            for ch in per_k:
                per_k[ch][ik]+=probability*data[ch]
                total_weight[ch][ik]+=probability*data[f"weight_{ch}"]
                elastic[ch][ik]+=probability*data[f"elastic_{ch}"]
            for moment in data["moment_checks"]:
                max_weight_error=max(max_weight_error,moment["weight_error"])
                max_first_error=max(max_first_error,moment["first_moment_error"])
            max_residual=max(max_residual,max(g["residual"] for g in data["metadata"]["ground_states"]))
            for pole in data["poles"]:
                ch=pole["channel"];q=pole["q_index"]
                first_moment[ch][ik,q]+=probability*pole["omega"]*pole["weight"]
                poles.append(pole|{"k_index":ik,"Nup":row["Nup"],"Ndn":row["Ndn"],
                                  "grand_probability":probability,"averaged_weight":float(wk*probability*pole["weight"])})
        counts.append(mean_N)
        variance_N=mean_N_squared-mean_N**2
        if variance_N < -1e-9:raise RuntimeError("Negative particle-number variance")
        connected_charge=elastic["charge"][ik].copy();connected_charge[0]=max(0.,variance_N)/L
        connected_spin=elastic["spin"][ik].copy();connected_spin[0]=mean_spin_squared/L
        k_report.append({"k_index":ik,"k":float(k),"weight":float(wk),"theta_boundary":float(L*k),
                         "phase_per_bond":float(-k),"mean_N":mean_N,"grand_ground_dimension":dimension,
                         "components":components,"elastic_charge_connected":connected_charge.tolist(),
                         "elastic_spin_connected":connected_spin.tolist()})
        print(f"U={U:g} response k={ik+1}/{len(k_values)}, mean N={mean_N:g}, ground states={dimension}",flush=True)
    averaged={ch:np.einsum('k,kqw->qw',weights,per_k[ch]) for ch in per_k}
    actual=float(weights@np.asarray(counts)/(2*L))
    selection.update(actual_filling=actual,filling_error=actual-target,exact_target_attained=abs(actual-target)<1e-10)
    connected={ch:np.einsum('k,kq->q',weights,np.asarray([r[f"elastic_{ch}_connected"] for r in k_report])) for ch in per_k}
    np.savez_compressed(output/"physical_density_response.npz",omega=omega,q=2*np.pi*np.arange(L)/L,
                        k=k_values,k_weights=weights,mean_N_per_k=np.asarray(counts),
                        S_charge=averaged["charge"],S_spin=averaged["spin"],
                        S_charge_per_k=per_k["charge"],S_spin_per_k=per_k["spin"],
                        elastic_charge_connected=connected["charge"],elastic_spin_connected=connected["spin"],
                        total_weight_charge=np.einsum('k,kq->q',weights,total_weight["charge"]),
                        total_weight_spin=np.einsum('k,kq->q',weights,total_weight["spin"]),
                        first_moment_charge=np.einsum('k,kq->q',weights,first_moment["charge"]),
                        first_moment_spin=np.einsum('k,kq->q',weights,first_moment["spin"]),
                        U=U,L=L,mu=mu,target_filling=target,actual_filling=actual,eta=eta)
    summary={"L":L,"U":U,"mu_selection":selection,"eta":eta,"k_count":len(k_values),
             "observable":"physical charge density and longitudinal 2Sz density",
             "normalization":"1/sqrt(L) source; normalized reduced-BZ weights sum to one",
             "ensemble":"T=0, complete unpolarized grand-canonical ground trace at each k",
             "q_scope":"only q=2*pi*m/L intra-block transfers; no generic inter-block response",
             "k_blocks":k_report,"checks":{"max_source_weight_error":max_weight_error,
             "max_first_moment_error":max_first_error,"max_ground_residual":max_residual,
             "q0_charge_inelastic_max":float(averaged["charge"][0].max()),
             "q0_spin_inelastic_max":float(averaged["spin"][0].max())},
             "new_solver_calls":runs.new_runs,"new_solver_elapsed_seconds":runs.elapsed,
             "qmin_first_charge_peak":first_peak(omega,averaged["charge"][1]),
             "qmin_first_spin_peak":first_peak(omega,averaged["spin"][1])}
    previous_file=output/"summary.json"
    previous=json.loads(previous_file.read_text()) if previous_file.exists() else {}
    summary['initial_solver_calls']=previous.get('initial_solver_calls',previous.get('new_solver_calls',runs.new_runs))
    summary['initial_solver_elapsed_seconds']=previous.get('initial_solver_elapsed_seconds',previous.get('new_solver_elapsed_seconds',runs.elapsed))
    reverse=(-np.arange(L))%L
    for ch in per_k:
        summary['checks'][f'{ch}_q_minus_q_average_max_difference']=float(np.max(abs(averaged[ch]-averaged[ch][reverse])))
    (output/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    with (output/"grand_canonical_poles.csv").open('w',newline='') as f:
        if poles:
            writer=csv.DictWriter(f,fieldnames=list(poles[0]));writer.writeheader();writer.writerows(poles)
    return summary


def plot_coupling(directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    directory=Path(directory)
    with np.load(directory/"physical_density_response.npz",allow_pickle=False) as d:
        L=int(d["L"]);U=float(d["U"]);n=float(d["actual_filling"]);mu=float(d["mu"])
        omega=d["omega"];q=d["q"];dw=omega[1]-omega[0];dq=2*np.pi/L
        fig,axes=plt.subplots(1,2,figsize=(8,3.2),layout="constrained")
        for ax,ch in zip(axes,["charge","spin"]):
            im=ax.imshow(d[f'S_{ch}'][:L//2+1].T,origin='lower',aspect='auto',interpolation='none',cmap='turbo',
                         extent=[-dq/2,np.pi+dq/2,omega[0]-dw/2,omega[-1]+dw/2])
            ax.set(title=ch.capitalize(),xlabel=r'$q$',ylabel=r'$\omega/t$',xlim=(0,np.pi),ylim=(omega[0],omega[-1]),
                   xticks=[0,np.pi/2,np.pi],xticklabels=['0',r'$\pi/2$',r'$\pi$'])
            fig.colorbar(im,ax=ax,label=r'$S(q,\omega)$')
        fig.suptitle(f'Physical density: L={L}, U/t={U:g}, n={n:.5f}, mu={mu:.5f}',fontsize=10)
        fig.savefig(directory/'physical_density_map.png',dpi=160);plt.close(fig)
        fig,ax=plt.subplots(figsize=(5,3.2),layout='constrained')
        for ch,color in [('charge','crimson'),('spin','navy')]:ax.plot(omega,d[f'S_{ch}'][1],label=ch,color=color)
        ax.set(xlabel=r'$\omega/t$',ylabel=r'$S(q_{\min},\omega)$',xlim=(0,5),title=f'L={L}, U={U:g}, common-mu k average')
        ax.legend();fig.savefig(directory/'physical_qmin_response.png',dpi=160);plt.close(fig)
        fig,ax=plt.subplots(figsize=(5,2.8),layout='constrained')
        ax.plot(L*d['k'],d['mean_N_per_k'],'o-',color='navy',ms=4)
        ax.set(xlabel=r'Boundary twist $\theta=Lk$',ylabel=r'$N_0(\theta)$',xticks=[-np.pi,0,np.pi],xticklabels=[r'$-\pi$','0',r'$\pi$'])
        fig.savefig(directory/'selected_particle_numbers.png',dpi=160);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--L',type=int,default=10)
    parser.add_argument('--couplings',type=float,nargs='+',default=[8.])
    parser.add_argument('--nk',type=int,default=16,help='Actual reduced-BZ point count, not the student input Nk')
    parser.add_argument('--grid',choices=['periodic','student','midpoint'],default='periodic')
    parser.add_argument('--target-filling',type=float,default=.45)
    parser.add_argument('--mu',type=float)
    parser.add_argument('--number-degeneracy',choices=['reject','uniform'],default='reject')
    parser.add_argument('--depth',type=int,default=180)
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--eta',type=float,default=.2)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--no-plots',action='store_true')
    parser.add_argument('--finite-ring-probe-q',type=float,nargs='*',default=[])
    args=parser.parse_args()
    if args.L<4 or args.L%2 or args.nk<2 or args.depth<3 or args.threads<1 or not np.isfinite(args.eta) or args.eta<=0 or not 0<=args.target_filling<=1:
        parser.error('Invalid ring, grid, filling or solver configuration')
    if any(U<0 or not np.isfinite(U) for U in args.couplings):parser.error('Require finite nonnegative U')
    if args.mu is not None and not np.isfinite(args.mu):parser.error('Require finite mu')
    if any(not np.isfinite(q) for q in args.finite_ring_probe_q):parser.error('Require finite probe momenta')
    output=(args.output or HERE/f'results/physical_mmhk_L{args.L}').resolve();output.mkdir(parents=True,exist_ok=True)
    exe=HERE/'bin/chain_Hubbard_tbc_response.x'
    if not exe.exists():parser.error('Build qbasis response solver first')
    source_sha=digest(ROOT/'create_code/chain_Hubbard_tbc_response.cc');solver_sha=digest(exe)
    config={"L":args.L,"couplings":args.couplings,"nk":args.nk,"grid":args.grid,"target_filling":args.target_filling,
            "mu":args.mu,"number_degeneracy":args.number_degeneracy,"depth":args.depth,"eta":args.eta,
            "source_sha256":source_sha,"executable_sha256":solver_sha}
    manifest=output/'config.json'
    if manifest.exists() and json.loads(manifest.read_text())!=config:raise RuntimeError('Existing output configuration differs; choose a new output directory')
    manifest.write_text(json.dumps(config,indent=2)+'\n')
    (output/'solver_source_snapshot.cc').write_bytes((ROOT/'create_code/chain_Hubbard_tbc_response.cc').read_bytes())
    os.environ.setdefault('MPLCONFIGDIR',str(HERE/'.mplconfig'))
    k,weights=k_grid(args.L,args.nk,args.grid);omega=np.linspace(-1,16,1701)
    for U in args.couplings:
        directory=output/f'U{U:g}'
        summary=calculate_coupling(directory,args.L,U,k,weights,omega,args.eta,args.target_filling,args.mu,args.depth,args.threads,
                                   args.number_degeneracy,solver_sha,source_sha)
        if not args.no_plots:plot_coupling(directory)
        if args.finite_ring_probe_q:
            with np.load(directory/'physical_density_response.npz',allow_pickle=False) as d:
                mixing=np.array([finite_ring_probe_weights(q,args.L) for q in args.finite_ring_probe_q])
                np.savez_compressed(directory/'finite_ring_window_probes.npz',probe_q=args.finite_ring_probe_q,mixing_weights=mixing,
                                    omega=d['omega'],charge=mixing@d['S_charge'],spin=mixing@d['S_spin'],
                                    interpretation='Finite-ring windowed density probes; NOT physical MMHK generic-q inter-block response')
        print('COMPLETE',U,'mu=',summary['mu_selection']['mu'],'achieved filling=',summary['mu_selection']['actual_filling'],summary['checks'],flush=True)
    analysis={p.name:digest(p) for p in [Path(__file__),HERE/'postprocess_response.py']}
    (output/'analysis_hashes.json').write_text(json.dumps(analysis,indent=2)+'\n')


if __name__=='__main__':main()
