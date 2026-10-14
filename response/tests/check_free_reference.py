#!/usr/bin/env python3
"""Free-band kinematic bound: a definition check on the supplied reference."""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.signal import find_peaks

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import postprocess_response as post


def main():
    L=10;t=1.;q=2*np.pi/L;bound=4*t*np.sin(q/2)
    # Verify the supplied per-bond phase convention by Fourier transformation.
    phase=.137;h=np.zeros((L,L),dtype=complex)
    for j in range(L):
        h[j,(j+1)%L]+=-t*np.exp(1j*phase)
        h[(j+1)%L,j]+=-t*np.exp(-1j*phase)
    K=2*np.pi*np.arange(L)/L;F=np.exp(-1j*np.arange(L)[:,None]*K)/np.sqrt(L)
    expected=np.diag(-2*t*np.cos(K-phase))
    bloch_error=float(np.max(abs(F.conj().T@h@F-expected)))
    assert bloch_error<1e-12
    run=HERE/"results/free_L10_U0_n0p45_phi0"
    omega=np.linspace(0,5,5001);ours=post.reconstruct(run,omega,.2)
    reference_path=next((HERE/"benchmark_reference/source/N=10_interaction_sweep").glob("structure_factors_Norb10_U0.00*_separated.npz"))
    result={"L":L,"t":t,"q":q,"physical_free_band_upper_bound":float(bound),
            "bloch_transform_max_error":bloch_error,"channels":{},
            "interpretation":"Inconsistent with the physical coherent density of the supplied epsilon(k)=-2t*cos(k) band; generating reference solver/settings unavailable."}
    with np.load(reference_path,allow_pickle=False) as ref:
        for channel in ["charge","spin"]:
            poles=[p["omega"] for p in ours["poles"] if p["channel"]==channel
                   and p["q_index"]==1 and p["weight"]>1e-8]
            assert max(poles)<=bound+1e-8
            np.testing.assert_allclose(ours["charge"][1],ours["spin"][1],atol=1e-10,rtol=1e-10)
            w=ref["omega_vals"];s=ref[f"S_{channel}_q_w"][1]
            mask=(w>0)&(w<5);indices,_=find_peaks(s[mask],height=.001*s[mask].max())
            peaks=w[mask][indices]
            result["channels"][channel]={"qbasis_distinct_poles":sorted(set(round(x,10) for x in poles)),
                                          "reference_local_maxima":peaks.tolist(),
                                          "reference_maxima_above_bound":peaks[peaks>bound].tolist()}
    (run/"free_reference_check.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
