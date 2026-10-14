#!/usr/bin/env python3
"""Zero-field thermodynamic Hubbard velocities from Essler et al.'s equations.

Book: (5.102)--(5.105), with velocities from (6.44)--(6.47).
Density n=N_e/L (half filling=1), u=U/(4t), t=a0=hbar=1 here.
This is an integral-equation calculation, NOT a finite-cluster/MMHK solver.
"""
import argparse
import csv
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import brentq
from scipy.special import digamma, expit, i0e, i1e, roots_legendre

HERE=Path(__file__).resolve().parent
BOOK=Path('/Users/cosdis/Downloads/books/Essler, Frahm, Gohmann, Klumper, Korepin - The One-Dimensional Hubbard Model.pdf')


def kernel_R(x,u):
    """Eq. (5.104), exact digamma representation of its Fourier integral."""
    x=np.asarray(x)
    z=1j*x/(4*u)
    return np.real(digamma(1+z)-digamma(.5+z))/(4*np.pi*u)


@lru_cache(maxsize=100)
def kernel_interpolant(u):
    # |sin(k')-sin(k)|<=2. Evaluate the special functions once per coupling,
    # then interpolate with spacing <= u/200. Independently checked below.
    samples=max(2049,int(np.ceil(400/u))+1)
    grid=np.linspace(0,2,samples)
    return CubicSpline(grid,kernel_R(grid,u),extrapolate=False)


def kernel_fast(x,u):
    return kernel_interpolant(u)(np.clip(np.abs(x),0,2))


@lru_cache(maxsize=16)
def quadrature(nodes):
    return roots_legendre(nodes)


def charge_system(Q,u,nodes):
    x,w0=quadrature(nodes);k=Q*x;w=Q*w0
    R=kernel_fast(np.sin(k)[None,:]-np.sin(k)[:,None],u)
    A=np.eye(nodes)-np.cos(k)[:,None]*R*w[None,:]
    lu=lu_factor(A,check_finite=False)
    rho=lu_solve(lu,np.full(nodes,1/(2*np.pi)),check_finite=False)
    return k,w,R,lu,rho


def solve(U_over_t,density=.9,nodes=400,return_state=False):
    U=float(U_over_t);n=float(density)
    if not 0<n<=1 or U<0:raise ValueError('Require 0<n<=1 and U/t>=0')
    if U==0:
        vf=2*np.sin(np.pi*n/2)
        return {'U_over_t':0.,'u':0.,'density':n,'Q':np.pi*n/2,
                'v_charge':float(vf),'v_spin':float(vf),'nodes':0,
                'density_error':0.,'rho_residual':0.,'energy_derivative_residual':0.,
                'kappa_residual':0.,'mu_book':float(-2*np.cos(np.pi*n/2)),
                'mu_bare_Hubbard':float(-2*np.cos(np.pi*n/2))}
    u=U/4
    # More nodes resolve the narrow kernel when U/t is small.
    N=int(np.ceil(nodes*max(1,np.sqrt(.5/U))))
    N+=N%2
    last={}
    def density_residual(Q):
        k,w,R,lu,rho=charge_system(Q,u,N)
        last.update(Q=Q,k=k,w=w,R=R,lu=lu,rho=rho)
        return float(w@rho-n)
    if n==1:Q=np.pi
    else:
        # For a positive interaction the endpoint Q=pi has exactly n=1.
        # Evaluate the physical endpoint to detect insufficient quadrature.
        if density_residual(np.pi)<0:raise ValueError('Cannot bracket density; increase quadrature nodes')
        Q=brentq(lambda q:-n if q==0 else density_residual(q),0,np.pi,xtol=2e-12,rtol=1e-13)
    if last.get('Q')!=Q:density_residual(Q)
    k,w,R,lu,rho=(last[key] for key in ['k','w','R','lu','rho'])
    # Differentiating Eq. (5.103), integrating by parts, using kappa(+/-Q)=0:
    # y(k)=kappa'(k)=2 sin(k)+cos(k)*int R(sin k'-sin k)y(k')dk'.
    y=lu_solve(lu,2*np.sin(k),check_finite=False)
    rb=kernel_fast(np.sin(k)-np.sin(Q),u)
    rhoQ=1/(2*np.pi)+np.cos(Q)*(rb*w)@rho
    yQ=2*np.sin(Q)+np.cos(Q)*(rb*w)@y
    vc=float(yQ/(2*np.pi*rhoQ))
    # A->infinity in (6.45): use asymptotic sech kernel, cancelling its
    # exponentially small factor analytically instead of dividing two zeros.
    # epsilon1'(Lambda)=int s(Lambda-sin k)*kappa'(k)dk.
    a=np.pi/(2*u);smax=1. if Q>=np.pi/2 else np.sin(Q)
    e=np.exp(a*(np.sin(k)-smax))
    vs=float((w*e)@y/(2*np.pi*((w*e)@rho)))
    # Reconstruct kappa and mu explicitly as an independent consistency check
    # on the dressed-energy equation (5.103). A_kappa = W^-1 A_rho^T W.
    driving=np.column_stack([-2*np.cos(k),np.ones(N)])
    fg=lu_solve(lu,w[:,None]*driving,trans=1,check_finite=False)/w[:,None]
    f,g=fg[:,0],fg[:,1]
    C=float((-2*np.cos(Q)+(rb*w*np.cos(k))@f)/(1+(rb*w*np.cos(k))@g))
    kappa=f-C*g
    kres=kappa-(-2*np.cos(k)-C+R@(w*np.cos(k)*kappa))
    rres=rho-(1/(2*np.pi)+np.cos(k)*(R@(w*rho)))
    yres=y-(2*np.sin(k)+np.cos(k)*(R@(w*y)))
    if np.min(rho)<-1e-9 or np.max(kappa)>1e-8:
        raise ValueError('Unphysical charge density or occupied-sea dressed energy')
    if vc<-1e-8 or vs<0:raise ValueError('Negative velocity')
    result={'U_over_t':U,'u':u,'density':n,'Q':float(Q),
            'v_charge':max(vc,0.),'v_spin':vs,'nodes':N,
            'density_error':float(w@rho-n),
            'rho_residual':float(np.max(abs(rres))),
            'energy_derivative_residual':float(np.max(abs(yres))),
            'kappa_residual':float(np.max(abs(kres))),
            'mu_book':C-2*u,'mu_bare_Hubbard':C}
    if return_state:
        return result,{'k':k,'w':w,'rho':rho,'kappa':kappa,'C':C,'u':u,'Q':Q}
    return result


def validation(density=.9,nodes=400):
    report={'kernel_integral':[],'kernel_interpolation':[],'mesh_convergence':[],
            'dressed_energy_velocity_checks':[],'half_filled_spin':[]}
    sample_x=np.random.default_rng(123).uniform(-2,2,3001)
    for u in [.025,.0625,.25,1.,2.,4.]:
        error=float(np.max(abs(kernel_fast(sample_x,u)-kernel_R(sample_x,u))))
        assert error<1e-9,(u,error)
        report['kernel_interpolation'].append({'u':u,'max_absolute_error':error})
    for u in [.05,.25,1.,3.]:
        for x in [0.,.3,1.7]:
            integral,error=quad(lambda w:np.cos(w*x)*expit(-2*u*w),0,np.inf,epsabs=2e-11,limit=800)
            numeric=integral/np.pi;analytic=float(kernel_R(x,u));delta=abs(analytic-numeric)
            assert delta<1e-8,(u,x,delta,error)
            report['kernel_integral'].append({'u':u,'x':x,'absolute_error':delta})
    for U in [.1,.25,1.,4.,8.,16.]:
        low=solve(U,density,nodes);high=solve(U,density,int(nodes*1.5))
        errors={c:abs(high[c]-low[c]) for c in ['v_charge','v_spin']}
        assert max(errors.values())<2e-7,(U,errors)
        report['mesh_convergence'].append({'U_over_t':U,'nodes':low['nodes'],
                                         'refined_nodes':high['nodes'],'absolute_errors':errors})
    for U in [.1,.25,4.,8.]:
        row,state=solve(U,density,nodes,return_state=True)
        k,w,rho,kappa,C,u,Q=(state[key] for key in ['k','w','rho','kappa','C','u','Q'])
        def energy_at(p):
            return -2*np.cos(p)-C+(w*np.cos(k)*kernel_R(np.sin(k)-np.sin(p),u))@kappa
        h=1e-5
        slope=(energy_at(Q+h)-energy_at(Q-h))/(2*h)
        rhoQ=1/(2*np.pi)+np.cos(Q)*(w*kernel_R(np.sin(k)-np.sin(Q),u))@rho
        direct_vc=float(slope/(2*np.pi*rhoQ))
        a=np.pi/(2*u);smax=1. if Q>=np.pi/2 else np.sin(Q)
        Lambda=smax+20/a;x=a*(Lambda-np.sin(k))
        sech=2*np.exp(-x)/(1+np.exp(-2*x))
        s=sech/(4*u)
        sigma=float((w*s)@rho)
        eps_prime=float((w*np.cos(k)*s*(-a*np.tanh(x)))@kappa)
        direct_vs=eps_prime/(2*np.pi*sigma)
        dc=abs(direct_vc-row['v_charge']);ds=abs(direct_vs-row['v_spin'])
        assert dc<2e-7 and ds<2e-7,(U,dc,ds)
        report['dressed_energy_velocity_checks'].append({'U_over_t':U,
                'charge_finite_difference_error':dc,'spin_direct_asymptotic_error':ds})
    for U in [4.,8.,12.]:
        row=solve(U,1.,nodes);a=2*np.pi/U;exact=2*i1e(a)/i0e(a)
        delta=abs(row['v_spin']-exact);assert delta<2e-8,(U,delta)
        report['half_filled_spin'].append({'U_over_t':U,'computed':row['v_spin'],
                                        'exact_Bessel':float(exact),'absolute_error':delta})
    free=solve(0.,density,nodes);report['noninteracting_velocity']=free['v_charge']
    report['strong_coupling']={
        'charge_limit':float(2*np.sin(np.pi*density)),
        'spin_leading_U_times_velocity':float(2*np.pi*(1-np.sin(2*np.pi*density)/(2*np.pi*density)))}
    strong=solve(1000.,density,nodes)
    assert abs(strong['v_charge']-report['strong_coupling']['charge_limit'])<.01
    assert abs(1000*strong['v_spin']-report['strong_coupling']['spin_leading_U_times_velocity'])<.1
    report['strong_coupling']['check_U_over_t']=1000.
    report['strong_coupling']['check_v_charge']=strong['v_charge']
    report['strong_coupling']['check_U_times_v_spin']=1000*strong['v_spin']
    report['passed']=True
    return report


def density_tag(density):
    return f'{density:.6f}'.rstrip('0').rstrip('.').replace('.','p')


def plot(rows,density,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    U=np.array([r['U_over_t'] for r in rows]);vc=np.array([r['v_charge'] for r in rows]);vs=np.array([r['v_spin'] for r in rows])
    plt.rcParams.update({'font.size':11})
    fig,ax=plt.subplots(figsize=(6.1,4.1),layout='constrained')
    ax.plot(U,vc,color='#c43c64',lw=2,label=r'Charge $v_c$')
    ax.plot(U,vs,color='#234f91',lw=2,label=r'Spin $v_s$',ls='--')
    selected=np.isin(U,[4,6,7,8,12])
    ax.plot(U[selected],vc[selected],'o',color='#c43c64',ms=4)
    ax.plot(U[selected],vs[selected],'s',color='#234f91',ms=4)
    ax.plot([0],[vc[0]],'o',color='#444444',ms=5)
    ax.set(xlabel=r'$U/t$',ylabel=r'Velocity $v/(t a_0/\hbar)$',xlim=(0,float(U.max())),ylim=(0,max(vc.max(),vs.max())*1.08))
    ax.set_title(rf'1D Hubbard model: $n=N_e/L={density:.2f}$, $T=0$, $B=0$',fontsize=12)
    ax.text(.97,.94,'Thermodynamic Bethe ansatz\nEssler et al., Secs. 5.5.4 and 6.5',transform=ax.transAxes,ha='right',va='top',fontsize=9,color='#555555')
    ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper right',bbox_to_anchor=(1,.79),frameon=False)
    fig.savefig(output/f'velocities_n{density_tag(density)}.png',dpi=180)
    fig.savefig(output/f'velocities_n{density_tag(density)}.svg')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--density',type=float,default=.9)
    p.add_argument('--U-max',type=float,default=16.)
    p.add_argument('--U-step',type=float,default=.25)
    p.add_argument('--nodes',type=int,default=400)
    p.add_argument('--output',type=Path)
    p.add_argument('--skip-validation',action='store_true')
    args=p.parse_args()
    if not 0<args.density<1 or args.U_max<=0 or args.U_step<=0 or args.nodes<50:p.error('Invalid density, U range or quadrature')
    output=(args.output or HERE/f'results_n{density_tag(args.density)}').resolve()
    output.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    couplings=np.unique(np.r_[np.arange(0,args.U_max+1e-9,args.U_step),[.1] if args.U_max>=.1 else [],args.U_max])
    rows=[]
    for U in couplings:
        row=solve(U,args.density,args.nodes);rows.append(row)
        print(f"U/t={U:5.2f} Q={row['Q']:.9f} vc={row['v_charge']:.9f} vs={row['v_spin']:.9f}",flush=True)
    with (output/'velocities.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    checks={} if args.skip_validation else validation(args.density,args.nodes)
    checks['maximum_equation_residuals']={k:max(abs(r[k]) for r in rows) for k in ['density_error','rho_residual','energy_derivative_residual','kappa_residual']}
    (output/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    plot(rows,args.density,output)
    (output/'metadata.json').write_text(json.dumps({
        'density_Ne_over_L':args.density,'earlier_filling_x':args.density/2,
        'magnetic_field':0,'temperature':0,'t':1,'a0':1,'hbar':1,
        'regime':'repulsive Hubbard thermodynamic limit; no finite-cluster solver executed',
        'book_path':str(BOOK),'book_sha256':hashlib.sha256(BOOK.read_bytes()).hexdigest(),
        'book_equations':['5.102','5.103','5.104','5.105','6.44','6.45','6.46','6.47'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'U_over_t_range':[0,args.U_max],'elapsed_seconds':time.monotonic()-start,
        'validated':bool(checks.get('passed',False))},indent=2)+'\n')
    print('Output:',output,'elapsed_seconds:',time.monotonic()-start,flush=True)


if __name__=='__main__':main()
