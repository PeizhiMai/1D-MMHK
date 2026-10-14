# Operator, momentum, and ensemble contract

## Hubbard block / twist map

For t=1, use the unshifted interaction U*n_up*n_down and

H(phi) = -sum_(j,sigma) [exp(i*phi) c†_(j,sigma)c_(j+1,sigma) + h.c.]
         + U*sum_j n_(j,up)n_(j,down).

The L-site dual-space ring is one reduced-BZ block of this 1D-MMHK construction.
The orbital index j here is the dual-space ring coordinate. Check the mapping
if a different mixing geometry or physical embedding is introduced.

- c_K = sum_j exp(+i*K*j)c_j/sqrt(L), K=2*pi*m/L.
- epsilon(K)=-2*cos(K-phi).
- Reduced-BZ coordinate k maps to phi=-k; native boundary twist theta=L*phi=-L*k.
- Row(last,0) boundary hopping has exp(+i*theta); row(0,last) has exp(-i*theta).
- The periodic grid k=2*pi/L*(arange(Nk)-Nk//2)/Nk counts only one identified
  RBZ endpoint. Even Nk contains zero. Nk=1 means k=phi=theta=0.
- The qbasis translation label for a free Slater state is Q=-sum(occupied K
  indices) mod L in this Fourier convention, not the opposite sign.

## Physical coherent density

rho_c/s(q) = sum_j exp(+i*q*j)(n_(j,up) +/- n_(j,down))/sqrt(L)
           = sum_K(c†_(K,up)c_(K+q,up) +/- c†_(K,down)c_(K+q,down))/sqrt(L).

S_nu(q,omega) = sum_(g,m) p_g |<m|rho_nu(q)|g>|² delta(omega-(Em-Eg)).

The orbital sum is **inside** the matrix element before squaring. The native
spin source is 2Sz; convert S_2Sz to S_Sz by dividing intensity by four.
The neutral target keeps Nup/Ndn unchanged and has Q_target=Q_ground+m for
q_m=2*pi*m/L. q=0 is a conserved total; its inelastic spectrum must vanish.

Supported intra-block transfers are q_m=2*pi*m/L. An off-grid Fourier probe
on a finite ring is a windowed mixture of these responses, **not** the generic
MMHK inter-block density response. Such a generic q may be well-defined, but
needs the appropriate inter-block operators and joint dynamics, not relabelled
intra-block spectra. The present calculation does not implement that extension.

The student's written orbital-pair intensity sum and the student's actual
Fourier-mode occupation source differ from this physical source and from each
other. For an audit, read the project's `response/STUDENT_CODE_AUDIT.md`; do not
claim Figure 3 reproduction by changing only normalization or peak choices.

## Actual ground ensemble and common chemical potential

At each twist choose all minimizers of E0(theta,N)-mu*N; one mu is shared by
all blocks for a given U and target density. The average physical density is
sum_theta w_theta*<N>_theta/L. Do not constrain each block separately unless a
canonical calculation is specifically requested.

Minimum-|Sz| scans contain every SU(2) multiplet of that N. Resolve all ground
Q and within-Q degeneracies, then include higher-|Sz| sectors until their
ground energies exceed the global minimum. Each ground vector has equal
probability in the identity-projector trace. Negative M can be represented by
spin-flip multiplicity two: the longitudinal charge/spin intensity is identical.
Never equal-weight magnetization sectors of different dimensions.

At a number crossing the ensemble is not unique. A full identity ground trace,
a specified statistical mixture, and a selected pure degenerate ground state
can give different densities/responses. Distinguish them explicitly; a target
at a crossing is not a stable mu plateau. The U=0 helper uses independent
occupation f=1/2 for exactly degenerate Fermi modes, i.e. the full grand ground
projector, and accepts it only when its mean density matches the target.

As a regression example at L=10, n=.8, U=4..12: zero-twist Ne=8 grounds can be
triplets. For 20 twists a common mu happened to select N=8 at every twist; that
was a result, not an imposed constraint. At U=7 the n=.9 stable plateau is
absent (jump .89→.91); it is a finite-grid step, not evidence by itself for
thermodynamic phase separation.

## Broadening and response definition

Use delta_eta(x)=eta/[pi*(x²+eta²)] for the plotted positive Lehmann spectrum.
Keep elastic weight separate. Finite-eta -Im(chi_R)/pi is not this one-sided
Lorentzian sum: it also contains the negative-frequency commutator term.
For the convention where the positive source is rho(q),

chi_R(q,z) = sum W_q/(z-Omega) - sum W_(-q)/(z+Omega).

At T=0 and eta→0 their positive-frequency inelastic parts agree. Both density
operators preserve N, so grand energy differences are Em-Eg with no mu shift.
Do not reuse the single-particle particle/hole frequency prescription.
