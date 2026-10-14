# Essler Bethe-ansatz velocities at density n=0.9

This is a **thermodynamic-limit** calculation for the repulsive 1D Hubbard
model at T=0, B=0 and n=N_e/L=0.9. In the earlier MMHK filling convention this
is x=n/2=0.45. The book's coupling is u=U/(4t), not U/t. The plot uses U/t.
No student-code, finite-cluster, or MMHK solver was executed. The open LaTeX
response-definition note was not changed.

## Equations and numerical method

Source: *The One-Dimensional Hubbard Model*, Essler, Frahm, Gohmann, Klumper
and Korepin. The relevant printed pages are 166-167 (PDF pages 184-185) and
187-188 (PDF pages 205-206).

- Eq. (5.102) fixes Q through `n = integral[-Q,Q] rho(k) dk`.
- Eqs. (5.103)-(5.105) give dressed energies, R and root densities.
- Eqs. (6.44)-(6.47) give velocities as **dressed energy / dressed momentum
  derivatives**; bare rapidity derivatives alone are not velocities.

Use

```
rho(k) = 1/(2*pi) + cos(k)*integral R(sin(k')-sin(k))*rho(k') dk'
y(k)   = 2*sin(k) + cos(k)*integral R(sin(k')-sin(k))*y(k') dk'
v_c    = y(Q)/(2*pi*rho(Q))
```

Here y=kappa'. Its equation follows by differentiating (5.103) and integrating
by parts with kappa(+/-Q)=0. The spin Fermi boundary is A=infinity at zero field.
To avoid numerically dividing two exponentially small functions, take the
asymptotic limit analytically:

```
v_s = integral exp[pi*sin(k)/(2u)]*y(k) dk
      / {2*pi*integral exp[pi*sin(k)/(2u)]*rho(k) dk}.
```

This follows from epsilon_1'(Lambda)=integral s(Lambda-sin(k))*kappa'(k) dk and
sigma_1(Lambda)=integral s(Lambda-sin(k))*rho(k) dk, with
s(x)=1/[4u*cosh(pi*x/(2u))]. A common exponential factor is removed for stability.

R is evaluated using its exact digamma representation:

```
R(x) = Re[psi(1+i*x/(4u)) - psi(1/2+i*x/(4u))]/(4*pi*u).
```

The solver uses Gauss-Legendre quadrature, a scalar root search for Q, and
linear solves. A dense one-dimensional interpolation of R accelerates the
matrix assembly; its error is checked against the exact special function.
Kappa and the chemical potential are reconstructed independently from (5.103).
The CSV's `mu_book` is the parameter in (5.103); for an unshifted onsite-U
Hamiltonian its effective chemical potential is `mu_bare_Hubbard=mu_book+2u`.

## Results

Velocities are in units ta0/hbar; t=a0=hbar=1 in the calculation.

| U/t | v_c | v_s |
| --- | --- | --- |
| 0 | 1.975377 | 1.975377 |
| 4 | 1.749791 | 1.272087 |
| 6 | 1.358547 | 0.985824 |
| 7 | 1.232506 | 0.877538 |
| 8 | 1.138850 | 0.787995 |
| 12 | 0.933448 | 0.552389 |

The output covers 0<=U/t<=16 in steps of 0.25, with an extra U/t=0.1 point.
The U=0 result is supplied analytically: v_c=v_s=2*sin(pi*n/2).

## Reproduce

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/.mplconfig \
~/.venvs/myenv/bin/python \
  /Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/essler_bethe_ansatz/hubbard_velocities.py \
  --density 0.9 --U-max 16 --U-step 0.25 --nodes 400
```

The default output directory is
`/Users/cosdis/Desktop/projects/HK_orbital_cluster/1d/response/essler_bethe_ansatz/results_n0p9`:

- `velocities_n0p9.png` and `.svg`: plot.
- `velocities.csv`: complete numerical data, Q, effective mu and residuals.
- `validation.json`: kernel, mesh, dressed-derivative and limiting checks.
- `metadata.json`: density convention, units, source/book hashes and equation list.

Validation checks include independent Fourier integration of R, higher-order
quadrature, direct finite differences of kappa at Q, the direct large-Lambda
spin-energy derivative, and the half-filled Bessel result
v_s(n=1)=2*I1(2*pi/(U/t))/I0(2*pi/(U/t)). Also checked the strong-coupling limits
v_c -> 2*sin(pi*n) and U*v_s -> 2*pi*[1-sin(2*pi*n)/(2*pi*n)].

Numerical equation residuals are below 4e-15 and density residuals below 2e-13
on the plotted grid. Quadrature refinement changes tested velocities by less
than 2e-13. The tested kernel interpolation error is below 3e-11. These are
numerical checks, not uncertainty estimates for finite-MMHK peak extractions.
