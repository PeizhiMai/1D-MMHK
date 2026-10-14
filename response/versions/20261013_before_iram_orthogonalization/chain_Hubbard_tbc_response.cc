// Physical longitudinal spin/charge response of one 1D n-MMHK block.
// The dual-space block is an L-site Hubbard ring with a uniform per-bond phase.
// rho_c/s(q) = sum_j exp(+i*q*j) (n_up +/- n_dn) / sqrt(L).
// This is a coherent density operator, not a sum of individual orbital-pair
// intensities. n = (Nup+Ndn)/(2L), so n=0.5 is half-filling.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "qbasis.h"

using C = std::complex<double>;
using Model = qbasis::model<C>;
using Operator = qbasis::mopr<C>;
namespace fs = std::filesystem;
constexpr double pi = 3.1415926535897932384626433832795;

struct Config {
    int L=10, nup=-1, ndn=-1, q=-1, ground_q=-1;
    int dense_ground_limit=128;
    double U=8.0, t=1.0, filling=0.45, phase=0.0;
    MKL_INT depth=300;
    bool ground_only=false;
    std::string observable="density";
    fs::path output="response_output";
};

Config parse(int argc, char** argv) {
    Config c;
    for (int i=1; i<argc; ++i) {
        const std::string key=argv[i];
        if (key=="--help") {
            std::cout << "Usage: chain_Hubbard_tbc_response.x [--L 10] [--U 8] [--t 1]\n"
              "  [--filling 0.45 | --nup 5 --ndn 4] [--phase 0]\n"
              "  [--q INDEX (default all 0..L-1)] [--ground-q INDEX (default scan)]\n"
              "  [--depth 300] [--output DIR]\n"
              "  [--observable density|orbital-pairs]\n"
              "  [--ground-only] (energies only, for chemical-potential selection)\n"
              "  [--dense-ground-limit 128] (ground-prefix solver crossover)\n"
              "phase is radians PER BOND; boundary twist theta=L*phase. No 1e-5 offset.\n"
              "filling = N/(2L), half-filling=0.5. Noninteger N is rejected.\n";
            std::exit(0);
        }
        if (key=="--ground-only") {c.ground_only=true;continue;}
        if (++i==argc) throw std::invalid_argument("Missing value for "+key);
        const std::string value=argv[i];
        if (key=="--output") { c.output=value; continue; }
        if (key=="--observable") { c.observable=value; continue; }
        // Check full conversion: do not silently truncate malformed arguments.
        size_t used=0;
        const double v=std::stod(value, &used);
        if (used!=value.size() || !std::isfinite(v)) throw std::invalid_argument("Invalid "+key);
        auto integer=[&]() {
            if (v!=std::round(v) || std::abs(v)>1000000)
                throw std::invalid_argument("Expected integer for "+key);
            return int(v);
        };
        if (key=="--L") c.L=integer();
        else if (key=="--nup") c.nup=integer();
        else if (key=="--ndn") c.ndn=integer();
        else if (key=="--q") c.q=integer();
        else if (key=="--ground-q") c.ground_q=integer();
        else if (key=="--dense-ground-limit") c.dense_ground_limit=integer();
        else if (key=="--depth") c.depth=integer();
        else if (key=="--U") c.U=v;
        else if (key=="--t") c.t=v;
        else if (key=="--filling") c.filling=v;
        else if (key=="--phase") c.phase=v;
        else throw std::invalid_argument("Unknown option "+key);
    }
    if (c.L<4 || c.L%2 || c.depth<3 || c.t<=0 || c.U<0 || c.filling<0 || c.filling>1)
        throw std::invalid_argument("Require even L>=4, depth>=3, t>0, U>=0, 0<=filling<=1");
    if (c.observable!="density" && c.observable!="orbital-pairs")
        throw std::invalid_argument("Observable must be density or orbital-pairs");
    if(c.dense_ground_limit<4 || c.dense_ground_limit>1024)
        throw std::invalid_argument("Dense ground limit must be 4..1024");
    if (c.nup < -1 || c.ndn < -1) throw std::invalid_argument("Negative Nup/Ndn");
    if ((c.nup<0)!=(c.ndn<0)) throw std::invalid_argument("Specify both --nup and --ndn");
    if (c.nup<0) {
        const double N=2*c.L*c.filling;
        if (std::abs(N-std::round(N))>1e-9)
            throw std::invalid_argument("Filling gives noninteger N in a single canonical block; supply exact Nup/Ndn");
        c.nup=(int(std::round(N))+1)/2;
        c.ndn=int(std::round(N))/2;
    }
    if (c.nup>c.L || c.ndn>c.L) throw std::invalid_argument("Nup/Ndn exceed L");
    if (c.q<-1 || c.q>=c.L || c.ground_q<-1 || c.ground_q>=c.L)
        throw std::invalid_argument("q and ground-q must be -1 or 0..L-1");
    c.filling=double(c.nup+c.ndn)/(2*c.L);
    return c;
}

struct Ground { int Q; double E, residual; std::vector<C> vector; int id=0; };

// Avoid qbasis simple-Lanczos's initial zero-beta issue in tiny sectors.
double solve_ground(Model& h, bool vector, std::vector<C>* saved=nullptr) {
    const MKL_INT dim=h.dim_repr[0];
    if (dim<=128) {
        auto a=h.HamMat_csr_repr[0].to_dense();
        std::vector<double> energies(dim);
        const MKL_INT info=LAPACKE_zheev(LAPACK_COL_MAJOR, vector ? 'V' : 'N', 'U', dim, a.data(), dim, energies.data());
        if (info) throw std::runtime_error("Small-sector diagonalization failed");
        if (vector && saved) saved->assign(a.begin(),a.begin()+dim);
        return energies[0];
    }
    h.locate_E0_lanczos(qbasis::which_sym::repr, 1, vector ? 1 : 0, 1000);
    if (vector && saved) *saved=h.eigenvecs_repr;
    return h.energy_min();
}

// Obtain the COMPLETE ground subspace within a translation sector, not one
// arbitrary vector. Dense sectors cover all L<=6 tests. Larger sectors request
// only a small IRAM prefix, enlarging it until a non-ground sentinel is found.
std::vector<Ground> ground_subspace(Model& h,int Q,double E0,double tolerance,int dense_limit) {
    const MKL_INT dim=h.dim_repr[0];
    std::vector<double> ev;
    std::vector<C> vectors;
    MKL_INT count=dim;
    if (dim<=std::max(6,dense_limit)) {
        vectors=h.HamMat_csr_repr[0].to_dense();ev.resize(dim);
        const MKL_INT info=LAPACKE_zheev(LAPACK_COL_MAJOR,'V','U',dim,vectors.data(),dim,ev.data());
        if(info)throw std::runtime_error("Ground subspace diagonalization failed");
    } else {
        for(MKL_INT requested=4;;requested*=2) {
            if(requested>64 || requested+2>dim)throw std::runtime_error("Ground degeneracy prefix exhausted");
            const MKL_INT ncv=std::min(dim,std::max(static_cast<MKL_INT>(20),2*requested+2));
            std::vector<C> start(dim);std::mt19937 engine(177+Q);
            std::normal_distribution<double> normal;
            for(auto& z:start)z=C(normal(engine),normal(engine));
            ev.assign(requested,0);vectors.assign(dim*requested,C(0));count=0;
            qbasis::iram(dim,h.HamMat_csr_repr[0],start.data(),requested,ncv,
                        static_cast<MKL_INT>(3000),"sr",count,ev.data(),vectors.data());
            if(count!=requested)throw std::runtime_error("Ground IRAM prefix did not converge");
            if(std::abs(*std::min_element(ev.begin(),ev.end())-E0)>tolerance)
                throw std::runtime_error("Ground energy disagrees between scan and IRAM");
            if(std::any_of(ev.begin(),ev.end(),[&](double e){return e>E0+tolerance;}))break;
        }
    }
    std::vector<Ground> result;
    for(MKL_INT i=0;i<count;++i)if(std::abs(ev[i]-E0)<tolerance) {
        Ground g{Q,ev[i],0,{}};
        g.vector.assign(vectors.begin()+i*dim,vectors.begin()+(i+1)*dim);
        std::vector<C> hv(dim);h.HamMat_csr_repr[0].MultMv(g.vector.data(),hv.data());
        double norm=0,error=0;
        for(MKL_INT j=0;j<dim;++j) {norm+=std::norm(g.vector[j]);error+=std::norm(hv[j]-g.E*g.vector[j]);}
        g.residual=std::sqrt(error);
        if(g.residual>1e-8 || std::abs(norm-1)>1e-9)throw std::runtime_error("Invalid ground vector");
        for(const auto& previous:result) {
            C overlap=0;for(MKL_INT j=0;j<dim;++j)overlap+=std::conj(previous.vector[j])*g.vector[j];
            if(std::abs(overlap)>1e-7)throw std::runtime_error("Ground vectors are not orthogonal");
        }
        result.push_back(std::move(g));
    }
    if(result.empty())throw std::runtime_error("Empty ground subspace");
    return result;
}

int main(int argc, char** argv) {
  try {
    const Config c=parse(argc,argv);
    if (fs::exists(c.output/"COMPLETE"))
        throw std::invalid_argument("Output directory contains a completed run; choose a new directory");
    fs::create_directories(c.output);
    // qbasis writes log_*.txt in cwd; isolate them from other runs.
    fs::current_path(fs::absolute(c.output));
    qbasis::initialize(false);
    std::cout << std::setprecision(17);
    qbasis::lattice lattice("chain",{uint32_t(c.L)},{"pbc"});
    Model h(lattice,5,std::max(100.0,c.U*c.L+4*c.t*c.L+10));
    h.add_orbital(c.L,"electron");
    std::vector<std::vector<C>> cup(4,std::vector<C>(4)),cdn=cup;
    cup[0][1]=1.0; cup[2][3]=1.0; cdn[0][2]=1.0; cdn[1][3]=-1.0;
    std::vector<Operator> up(c.L),dn(c.L);
    std::vector<qbasis::opr<C>> uops,dops;
    Operator Nup,Ndn;
    for (int j=0;j<c.L;++j) {
        uint32_t site,next; std::vector<int> work(lattice.dim);
        lattice.coor2site({j},0,site,work); lattice.coor2site({j+1},0,next,work);
        auto u=qbasis::opr<C>(site,0,true,cup),d=qbasis::opr<C>(site,0,true,cdn);
        uops.push_back(u);dops.push_back(d);
        auto ud=u,dd=d;ud.dagger();dd.dagger();
        up[j]+=ud*u;dn[j]+=dd*d;
        Nup+=up[j];Ndn+=dn[j];
        auto un=qbasis::opr<C>(next,0,true,cup),dnxt=qbasis::opr<C>(next,0,true,cdn);
        auto und=un,dnd=dnxt;und.dagger();dnd.dagger();
        const C hop=-c.t*std::exp(C(0,c.phase));
        h.add_Ham(hop*(ud*un)); h.add_Ham(std::conj(hop)*(und*u));
        h.add_Ham(hop*(dd*dnxt)); h.add_Ham(std::conj(hop)*(dnd*d));
        h.add_Ham(C(c.U,0)*(up[j]*dn[j]));
    }
    h.fill_Weisse_table();
    std::vector<double> energies(c.L,std::numeric_limits<double>::infinity());
    std::ofstream sectors("ground_sectors.csv");sectors.exceptions(std::ios::badbit|std::ios::failbit);
    sectors<<std::setprecision(18)<<"Q,E0,repr_dimension\n";
    for (int Q=0;Q<c.L;++Q) {
        if (c.ground_q>=0 && Q!=c.ground_q) continue;
        h.enumerate_basis_repr({Q},{Nup,Ndn},{double(c.nup),double(c.ndn)},0);
        if (!h.dim_repr[0] || std::none_of(h.norm_repr[0].begin(),h.norm_repr[0].end(),
                [](double v){return std::abs(v)>qbasis::lanczos_precision;})) continue;
        h.generate_Ham_sparse_repr(0);h.switch_sec_mat(0);
        energies[Q]=solve_ground(h,false);
        sectors<<Q<<","<<energies[Q]<<","<<h.dim_repr[0]<<"\n";sectors.flush();
    }
    const double E0=*std::min_element(energies.begin(),energies.end());
    if (!std::isfinite(E0)) throw std::runtime_error("No nonempty ground sector");
    std::vector<Ground> grounds;
    constexpr double ground_tol=1e-9;
    for (int Q=0;Q<c.L;++Q) if (std::abs(energies[Q]-E0)<ground_tol) {
        if(c.ground_only) {grounds.push_back(Ground{Q,energies[Q],0,{}});continue;}
        h.enumerate_basis_repr({Q},{Nup,Ndn},{double(c.nup),double(c.ndn)},0);
        h.generate_Ham_sparse_repr(0);h.switch_sec_mat(0);
        for(auto& g:ground_subspace(h,Q,E0,ground_tol,c.dense_ground_limit))grounds.push_back(std::move(g));
    }
    for(size_t i=0;i<grounds.size();++i)grounds[i].id=static_cast<int>(i);
    std::ofstream meta("metadata.json");meta.exceptions(std::ios::badbit|std::ios::failbit);meta<<std::setprecision(18);
    meta<<"{\n\"L\":"<<c.L<<",\"Nup\":"<<c.nup<<",\"Ndn\":"<<c.ndn
        <<",\"filling\":"<<c.filling<<",\"U\":"<<c.U<<",\"t\":"<<c.t
        <<",\"phase_per_bond\":"<<c.phase<<",\"boundary_twist\":"<<c.L*c.phase
        <<",\"boundary_phase_row0_last\":"<<-c.L*c.phase
        <<",\"ground_only\":"<<(c.ground_only?"true":"false")
        <<",\"ground_subspace_complete\":"<<(c.ground_only?"false":"true")
        <<",\"dense_ground_limit\":"<<c.dense_ground_limit
        <<",\"depth\":"<<c.depth<<",\"ground_scan_all_Q\":"<<(c.ground_q<0?"true":"false")
        <<",\"E0\":"<<E0<<",\"density_normalization\":\"1/sqrt(L)\""
        <<",\"observable\":\""<<c.observable<<"\""
        <<",\"spin_operator\":\"n_up-n_down (2Sz)\",\"q_indices\":[";
    for(int q=0;q<c.L && !c.ground_only;++q)if(c.q<0 || q==c.q) {
        if(c.q<0 && q)meta<<",";
        meta<<q;
    }
    meta<<"],\"ground_states\":[";
    for (size_t i=0;i<grounds.size();++i) {
        const auto& g=grounds[i];if(i)meta<<",";
        meta<<"{\"Q\":"<<g.Q<<",\"E0\":"<<g.E<<",\"id\":"<<g.id<<",\"residual\":";
        if(c.ground_only)meta<<"null";else meta<<g.residual;
        meta<<"}";
    }
    meta<<"]}\n";meta.close();
    // Each component is written separately; postprocessing averages equally
    // over exactly degenerate ground momenta. Spin-flipped sectors give the
    // same longitudinal spectra (global up/down exchange symmetry).
    if(!c.ground_only)for (const auto& g:grounds) {
      h.enumerate_basis_repr({g.Q},{Nup,Ndn},{double(c.nup),double(c.ndn)},0);
      h.eigenvecs_repr=g.vector;h.eigenvals_repr={g.E};h.nconv=1;
      for (int q=0;q<c.L;++q) {
        if (c.q>=0 && q!=c.q) continue;
        const int target=(g.Q+q)%c.L;
        if (q || c.observable=="orbital-pairs") {
            h.enumerate_basis_repr({target},{Nup,Ndn},{double(c.nup),double(c.ndn)},1);
            h.generate_Ham_sparse_repr(1);h.switch_sec_mat(1);
        }
        for (int spin=0;spin<2;++spin) {
          for (int orbital=0;orbital<(c.observable=="density"?1:c.L);++orbital) {
            const std::string channel=spin?"spin":"charge";
            const std::string tag=c.observable=="density"?"":"_orb"+std::to_string(orbital);
            std::ofstream out(channel+"_gsQ"+std::to_string(g.Q)+"_gsi"+std::to_string(g.id)+"_q"+std::to_string(q)+tag+".dat");
            out.exceptions(std::ios::badbit|std::ios::failbit);
            out<<std::setprecision(18)<<"E0 "<<g.E<<"\nq_index "<<q<<"\nq "<<2*pi*q/c.L
               <<"\nground_Q "<<g.Q<<"\ntarget_Q "<<target<<"\nensemble_weight "<<1.0/grounds.size()<<"\n";
            out<<"ground_id "<<g.id<<"\n";
            out<<"source_orbital "<<(c.observable=="density"?-1:orbital)<<"\n";
            if (!q && c.observable=="density") {
                // Fixed N and Sz: only elastic weight, exactly zero connected source.
                const double total=spin?c.nup-c.ndn:c.nup+c.ndn;
                out<<"norm 0\nelastic_weight "<<total*total/c.L<<"\nsteps 0\na\nb\n";
                continue;
            }
            Operator rho;
            if (c.observable=="density") {
                for(int j=0;j<c.L;++j) {
                    const C coeff=std::exp(C(0,2*pi*q*j/c.L))/std::sqrt(double(c.L));
                    rho+=coeff*up[j];rho+=(spin?-coeff:coeff)*dn[j];
                }
            } else {
                // Diagnostic for the OUTSIDE-the-modulus orbital sum printed
                // in Eq.(2) of the student's note: sum_K |<m|A_K|g>|^2 / L.
                // A_K=c_K^dagger c_(K+q) with c_K=sum_j exp(+iKj)c_j/sqrt(L).
                // NOT the physical coherent |<m|sum_K A_K|g>|^2/L.
                for(int i=0;i<c.L;++i)for(int j=0;j<c.L;++j) {
                    const C coeff=std::exp(C(0,2*pi*(-orbital*i+(orbital+q)*j)/c.L)) /
                                  (c.L*std::sqrt(double(c.L)));
                    auto ud=uops[i],dd=dops[i];ud.dagger();dd.dagger();
                    rho+=coeff*(ud*uops[j]);rho+=(spin?-coeff:coeff)*(dd*dops[j]);
                }
            }
            // Verify norm and first moment using an independent source application.
            std::vector<C> v(h.dim_repr[1]),hv(h.dim_repr[1]);
            h.moprXvec_repr(rho,0,1,static_cast<MKL_INT>(0),v.data());
            double w=0,first=0;h.HamMat_csr_repr[1].MultMv(v.data(),hv.data());
            for(size_t j=0;j<v.size();++j) {w+=std::norm(v[j]);first+=std::real(std::conj(v[j])*(hv[j]-g.E*v[j]));}
            std::vector<double> ab(2*c.depth,0);MKL_INT steps=0;double norm=0;
            // Native qbasis dynamic implementation, not a substitute solver.
            double variance=0;
            if (w>1e-24) {
                const double mean=g.E+first/w;
                for(size_t j=0;j<v.size();++j) variance+=std::norm(hv[j]-mean*v[j]);
            }
            if (w>1e-24 && variance/w<1e-22) {
                // An exactly one-pole source needs no recurrence. The existing
                // qbasis simple Lanczos divides by b1 before checking for zero.
                norm=std::sqrt(w);steps=1;ab[c.depth]=g.E+first/w;
            } else {
                h.measure_repr_dynamic(rho,0,1,c.depth,steps,norm,ab.data());
            }
            if(std::abs(norm*norm-w)>1e-9) throw std::runtime_error("Source norm mismatch");
            if(norm>qbasis::lanczos_precision && steps==0) throw std::runtime_error("Empty dynamic Krylov chain");
            out<<"norm "<<norm<<"\nelastic_weight 0\nsource_weight "<<w<<"\nfirst_moment "<<first
               <<"\nsteps "<<steps<<"\na";
            for(MKL_INT j=0;j<steps;++j)out<<" "<<ab[c.depth+j];
            out<<"\nb";for(MKL_INT j=0;j<steps;++j)out<<" "<<ab[j];out<<"\n";
            std::cout<<"RESPONSE gsQ="<<g.Q<<" q="<<q<<" channel="<<channel<<" orbital="<<orbital<<" weight="<<w<<" steps="<<steps<<"\n";
          }
        }
      }
    }
    std::ofstream complete("COMPLETE");complete.exceptions(std::ios::badbit|std::ios::failbit);complete<<"success\n";
    std::cout<<"RESPONSE_COMPLETE L="<<c.L<<" U="<<c.U<<" filling="<<c.filling<<" E0="<<E0<<"\n";
    return 0;
  } catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
