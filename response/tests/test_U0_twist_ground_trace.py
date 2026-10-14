#!/usr/bin/env python3
"""Compare exact free particle-hole response to qbasis full ground traces."""
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from twist_averaged_velocity_scan import Runs, ground_components, fractions
from postprocess_twist_velocity_scan import spectrum

ROOT = HERE/'results/twist20_common_mu_U0to12_n0p8'


def main():
    config = json.loads((ROOT/'config.json').read_text())
    source = ROOT/'U0/n0p8'
    d = np.load(source/'exact_free_reference.npz')
    k, epsilon, f, E0, mu = [d[key] for key in ['k','epsilon','f','E0_by_N','mu']]
    m = json.loads((source/'metadata.json').read_text())
    runs = Runs(ROOT/'U0/validation/qbasis', 0., 10, 300, 1, ROOT/'solver_executable_snapshot.x',
                config['solver_source_sha256'], config['solver_executable_sha256'])
    omega = np.linspace(0,2,4001); reports = []
    for i in [0,5,10]:
        support = [int(N) for N in m['twists'][i]['N_support']]
        components = []
        for N in support:
            components += ground_components(runs,i,k[i],N,E0[i,N])
        dimension = sum(c['spin_flip_multiplicity']*c['ground_dimension'] for c in components)
        assert dimension == m['twists'][i]['ground_dimension']
        mean_N = sum(c['N']*c['spin_flip_multiplicity']*c['ground_dimension']/dimension for c in components)
        np.testing.assert_allclose(mean_N,m['twists'][i]['N'],atol=1e-12)
        ep, wp = [], []
        for dest in range(10):
            src = (dest+1)%10
            prob = f[i,src]*(1-f[i,dest])
            if prob>0:
                ep.append(epsilon[i,dest]-epsilon[i,src]);wp.append(2*prob/10)
        expected = spectrum(omega,np.asarray(ep),np.asarray(wp),.2)
        case = {'twist_index':i,'theta':float(-10*k[i]),'ground_dimension':dimension,
                'mean_N':mean_N,'N_support':support,'channel_errors':{}}
        actual = {}
        for ch in ['charge','spin']:
            actual[ch] = np.zeros_like(omega)
            for c in components:
                E,W=fractions(c,ch,{},[])
                factor=c['spin_flip_multiplicity']*c['ground_dimension']/dimension
                actual[ch]+=factor*spectrum(omega,E,W,.2)
            error=float(np.max(abs(actual[ch]-expected)))
            assert error<1e-7,(i,ch,error)
            case['channel_errors'][ch]=error
        np.testing.assert_allclose(actual['charge'],actual['spin'],atol=1e-7,rtol=0)
        reports.append(case)
        print('PASS U0 actual-grand-ground trace:',case,flush=True)
    report={'passed':True,'mu':float(mu),'reference':'exact momentum-space free particle-hole sum with half-filled degenerate Fermi modes',
            'cases':reports,'zero_twist_grand_ground_dimension':16,'zero_twist_number_variance':1.,
            'spin_equals_charge':True,'actual_density':.8,'new_qbasis_calls_for_validation':runs.new_calls}
    (source/'free_qbasis_validation.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
