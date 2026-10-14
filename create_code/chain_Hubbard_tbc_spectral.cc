// Single-particle sources on the same TBC Hamiltonian/complete ground projector.
// Reuse the validated Hamiltonian ground-subspace functions without changing
// the neutral-response source or its executable.
#define main neutral_response_unused_main
#include "chain_Hubbard_tbc_response.cc"
#undef main

int main(int argc,char** argv) {
 try {
    const Config c=parse(argc,argv);
    if(c.ground_q<0)throw std::invalid_argument("Supply certified --ground-q from existing ground table");
    if(fs::exists(c.output/"COMPLETE"))throw std::invalid_argument("Completed output exists");
    fs::create_directories(c.output);fs::current_path(fs::absolute(c.output));
    qbasis::initialize(false);std::cout<<std::setprecision(17);
    qbasis::lattice lattice("chain",{uint32_t(c.L)},{"pbc"});
    Model h(lattice,5,std::max(100.0,c.U*c.L+4*c.t*c.L+10));h.add_orbital(c.L,"electron");
    std::vector<std::vector<C>> cup(4,std::vector<C>(4)),cdn=cup;
    cup[0][1]=1.;cup[2][3]=1.;cdn[0][2]=1.;cdn[1][3]=-1.;
    std::vector<Operator> up(c.L),dn(c.L);std::vector<qbasis::opr<C>> uops,dops;Operator Nup,Ndn;
    for(int j=0;j<c.L;++j) {
      uint32_t site,next;std::vector<int> work(lattice.dim);
      lattice.coor2site({j},0,site,work);lattice.coor2site({j+1},0,next,work);
      auto u=qbasis::opr<C>(site,0,true,cup),d=qbasis::opr<C>(site,0,true,cdn);
      uops.push_back(u);dops.push_back(d);auto ud=u,dd=d;ud.dagger();dd.dagger();
      up[j]+=ud*u;dn[j]+=dd*d;Nup+=up[j];Ndn+=dn[j];
      auto un=qbasis::opr<C>(next,0,true,cup),dnext=qbasis::opr<C>(next,0,true,cdn);
      auto und=un,dnd=dnext;und.dagger();dnd.dagger();const C hop=-c.t*std::exp(C(0,c.phase));
      h.add_Ham(hop*(ud*un));h.add_Ham(std::conj(hop)*(und*u));
      h.add_Ham(hop*(dd*dnext));h.add_Ham(std::conj(hop)*(dnd*d));h.add_Ham(C(c.U)*(up[j]*dn[j]));
    }
    h.fill_Weisse_table();
    h.enumerate_basis_repr({c.ground_q},{Nup,Ndn},{double(c.nup),double(c.ndn)},0);
    h.generate_Ham_sparse_repr(0);h.switch_sec_mat(0);const double E0=solve_ground(h,false);
    auto grounds=ground_subspace(h,c.ground_q,E0,1e-9,c.dense_ground_limit);
    for(size_t i=0;i<grounds.size();++i)grounds[i].id=int(i);
    std::ofstream meta("metadata.json");meta<<std::setprecision(18)
      <<"{\"L\":"<<c.L<<",\"U\":"<<c.U<<",\"Nup\":"<<c.nup<<",\"Ndn\":"<<c.ndn
      <<",\"phase_per_bond\":"<<c.phase<<",\"ground_Q\":"<<c.ground_q<<",\"E0\":"<<E0
      <<",\"depth\":"<<c.depth<<",\"ground_dimension\":"<<grounds.size()<<",\"grounds\":[";
    for(size_t i=0;i<grounds.size();++i){if(i)meta<<",";meta<<"{\"id\":"<<i<<",\"residual\":"<<grounds[i].residual<<"}";}
    meta<<"]}";meta.close();
    std::ofstream checks("sum_rules.csv");checks<<"ground_id,spin,q,particle_weight,hole_weight,weight_error,moment_error\n";
    for(const auto& g:grounds) {
      h.enumerate_basis_repr({g.Q},{Nup,Ndn},{double(c.nup),double(c.ndn)},0);
      h.eigenvecs_repr=g.vector;h.eigenvals_repr={g.E};h.nconv=1;
      for(int spin=0;spin<2;++spin) {
        // Balanced-sector COMPLETE ground projector has up/down exchange symmetry.
        if(spin && c.nup==c.ndn)continue;
        double total_hole=0;
        for(int q=0;q<c.L;++q) {
          if(c.q>=0 && c.q!=q)continue;
          double weights[2]={0,0},moments[2]={0,0};
          for(int branch=0;branch<2;++branch) {
            const bool particle=branch==0;int Nu=c.nup,Nd=c.ndn;
            (spin?Nd:Nu)+=particle?1:-1;
            const int target=(g.Q+(particle?-q:q)+c.L)%c.L;
            std::ofstream out(std::string(particle?"particle":"hole")+"_spin"+std::to_string(spin)+"_g"+std::to_string(g.id)+"_q"+std::to_string(q)+".dat");
            out<<std::setprecision(18)<<"E0 "<<g.E<<"\nq_index "<<q<<"\nspin "<<spin<<"\nground_id "<<g.id<<"\ntarget_Q "<<target<<"\n";
            if(Nu<0 || Nd<0 || Nu>c.L || Nd>c.L){out<<"norm 0\nfirst_moment 0\nsteps 0\na\nb\n";continue;}
            h.enumerate_basis_repr({target},{Nup,Ndn},{double(Nu),double(Nd)},1);
            h.generate_Ham_sparse_repr(1);h.switch_sec_mat(1);Operator source;
            for(int j=0;j<c.L;++j) {
              auto op=spin?dops[j]:uops[j];if(particle)op.dagger();
              source+=std::exp(C(0,(particle?-1:1)*2*pi*q*j/c.L))/std::sqrt(double(c.L))*op;
            }
            std::vector<C> v(h.dim_repr[1]),hv(v.size());h.moprXvec_repr(source,0,1,static_cast<MKL_INT>(0),v.data());
            h.HamMat_csr_repr[1].MultMv(v.data(),hv.data());double w=0,first=0,var=0;
            for(size_t j=0;j<v.size();++j){w+=std::norm(v[j]);first+=std::real(std::conj(v[j])*(hv[j]-g.E*v[j]));}
            if(w>1e-24)for(size_t j=0;j<v.size();++j)var+=std::norm(hv[j]-(g.E+first/w)*v[j]);
            std::vector<double> ab(2*c.depth);MKL_INT steps=0;double norm=0;
            if(w<1e-24){}else if(var/w<1e-22){norm=std::sqrt(w);steps=1;ab[c.depth]=g.E+first/w;}
            else h.measure_repr_dynamic(source,0,1,c.depth,steps,norm,ab.data());
            if(std::abs(norm*norm-w)>1e-9)throw std::runtime_error("Charged source norm mismatch");
            weights[branch]=w;moments[branch]=first;
            out<<"norm "<<norm<<"\nfirst_moment "<<first<<"\nsteps "<<steps<<"\na";
            for(MKL_INT i=0;i<steps;++i)out<<" "<<ab[c.depth+i];out<<"\nb";
            for(MKL_INT i=0;i<steps;++i)out<<" "<<ab[i];out<<"\n";
          }
          const double dw=std::abs(weights[0]+weights[1]-1),dm=std::abs(moments[0]-moments[1]-(-2*c.t*std::cos(2*pi*q/c.L-c.phase)+c.U*(spin?c.nup:c.ndn)/double(c.L)));
          checks<<g.id<<","<<spin<<","<<q<<","<<weights[0]<<","<<weights[1]<<","<<dw<<","<<dm<<"\n";
          if(dw>1e-8 || dm>1e-7)throw std::runtime_error("Charged fermion sum rule failed: q="+std::to_string(q)+" weight="+std::to_string(dw)+" moment="+std::to_string(dm));
          total_hole+=weights[1];
        }
        if(c.q<0 && std::abs(total_hole-(spin?c.ndn:c.nup))>1e-8)throw std::runtime_error("Integrated momentum occupation sum failed");
      }
    }
    std::ofstream("COMPLETE")<<"success\n";
    return 0;
 }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
