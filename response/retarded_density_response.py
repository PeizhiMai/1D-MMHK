#!/usr/bin/env python3
"""Full retarded susceptibility from saved physical density Lehmann poles.

The positive-frequency structure factor and the retarded commutator response
are different at finite eta. This tool preserves them as separate artifacts.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def evaluate_retarded(directory,omega,eta):
    directory=Path(directory)
    summary=json.loads((directory/'summary.json').read_text());L=summary['L']
    with (directory/'grand_canonical_poles.csv').open() as f:rows=list(csv.DictReader(f))
    result={ch:np.zeros((L,len(omega)),complex) for ch in ['charge','spin']}
    z=np.asarray(omega)[None,:]+1j*eta
    for ch in result:
        for q in range(L):
            for selected,sign in [(q,1),((-q)%L,-1)]:
                poles=[(float(row['omega']),float(row['averaged_weight'])) for row in rows
                       if row['channel']==ch and int(row['q_index'])==selected and float(row['omega'])>1e-8]
                for start in range(0,len(poles),128):
                    block=np.asarray(poles[start:start+128])
                    result[ch][q]+=np.sum(sign*block[:,1,None]/(z-sign*block[:,0,None]),axis=0)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directories',nargs='+',type=Path)
    p.add_argument('--omega-min',type=float,default=-16);p.add_argument('--omega-max',type=float,default=16)
    p.add_argument('--points',type=int,default=3201);p.add_argument('--eta',type=float,default=.2)
    args=p.parse_args()
    if not np.isfinite([args.eta,args.omega_min,args.omega_max]).all() or args.eta<=0 or args.points<3 or args.omega_max<=args.omega_min:
        p.error('Invalid frequency grid or eta')
    omega=np.linspace(args.omega_min,args.omega_max,args.points)
    for directory in args.directories:
        chi=evaluate_retarded(directory,omega,args.eta)
        np.savez_compressed(directory/'retarded_density_susceptibility.npz',omega=omega,eta=args.eta,
                            chi_charge=chi['charge'],chi_spin=chi['spin'],
                            loss_charge=-chi['charge'].imag/np.pi,loss_spin=-chi['spin'].imag/np.pi,
                            convention='chi_R= sum positive W_q/(z-Omega) - sum W_minus_q/(z+Omega)')
        checks={'q0_charge_max':float(np.max(abs(chi['charge'][0]))),'q0_spin_max':float(np.max(abs(chi['spin'][0])))}
        if np.allclose(omega,-omega[::-1]):
            L=chi['charge'].shape[0]
            for ch in chi:
                checks[f'{ch}_hermitian_frequency_relation_error']=float(np.max(abs(chi[ch]-chi[ch][(-np.arange(L))%L,::-1].conj())))
        (directory/'retarded_response_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        print(directory,checks)


if __name__=='__main__':main()
