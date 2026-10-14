#include <cassert>
#include <iostream>
#include <iomanip>
#include "qbasis.h"
#include <cmath>
#include <fstream>
#include <cstdlib>  // For std::atoi

#define PI 3.1415926535897932

// Fermi-Hubbard model on square lattice
int main(int argc, char *argv[]) {
    // Check if at least two arguments are provided (not counting the program name)
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << "N_up N_dn <kx> <ky> qx qy" << std::endl;
        return 1;
    }

    // Read the first two command-line arguments into kx and ky
    int Nupinput = std::atoi(argv[1]);
    int Ndninput = std::atoi(argv[2]);
    int kx = std::atoi(argv[3]);  // Convert the first argument to an integer
    int ky = std::atoi(argv[4]);  // Convert the second argument to an integer
    int Qm = std::atoi(argv[5]);  // Convert the third argument to an integer
    int Qn = std::atoi(argv[6]);  // Convert the fourth argument to an integer

    qbasis::initialize();
    std::cout << std::setprecision(10);
    // parameters
    double t = 1;
    double tp = 0.0;
    double U = 6.0;
    int Lx = 4;
    int Ly = 4;
    double Nup_total = 1.0*Nupinput;
    double Ndn_total = 1.0*Ndninput;

    double fx = kx * PI/40+0.00001;
    double fy = ky * PI/40+0.00001;
    double sx = sin(fx);
    double cx = cos(fx);
    double sy = sin(fy);
    double cy = cos(fy);
    double cxpy = cos(fx+fy);
    double sxpy = sin(fx+fy);
    double cxmy = cos(fx-fy);
    double sxmy = sin(fx-fy);

    //int step = 150;

    std::cout << "Lx =      " << Lx << std::endl;
    std::cout << "Ly =      " << Ly << std::endl;
    std::cout << "t =       " << t << std::endl;
    std::cout << "U =       " << U << std::endl;
    std::cout << "N_up =    " << Nup_total << std::endl;
    std::cout << "N_dn =    " << Ndn_total << std::endl << std::endl;

    // lattice object
    std::vector<std::string> bc{"pbc", "pbc"};
    qbasis::lattice lattice("square",{static_cast<uint32_t>(Lx), static_cast<uint32_t>(Ly)},bc);


    // local matrix representation
    auto c_up = std::vector<std::vector<std::complex<double>>>(4,std::vector<std::complex<double>>(4, 0.0));
    auto c_dn = std::vector<std::vector<std::complex<double>>>(4,std::vector<std::complex<double>>(4, 0.0));
    c_up[0][1] = std::complex<double>(1.0,0.0);
    c_up[2][3] = std::complex<double>(1.0,0.0);
    c_dn[0][2] = std::complex<double>(1.0,0.0);
    c_dn[1][3] = std::complex<double>(-1.0,0.0);


    // constructing the Hamiltonian in operator representation
    qbasis::model<std::complex<double>> Hubbard(lattice);
    Hubbard.add_orbital(lattice.Nsites, "electron");
    qbasis::mopr<std::complex<double>> Nup;   // an operator representating total electron number
    qbasis::mopr<std::complex<double>> Ndown;
    for (int x = 0; x < Lx; x++) {
        for (int y = 0; y < Ly; y++) {
            uint32_t site_i, site_j;
            std::vector<int> work(lattice.dim);
            lattice.coor2site({x,y}, 0, site_i, work); // obtain site label of (x,y)
            // construct operators on each site
            auto c_up_i    = qbasis::opr<std::complex<double>>(site_i,0,true,c_up);
            auto c_dn_i    = qbasis::opr<std::complex<double>>(site_i,0,true,c_dn);
            auto c_up_dg_i = c_up_i; c_up_dg_i.dagger();
            auto c_dn_dg_i = c_dn_i; c_dn_dg_i.dagger();
            auto n_up_i    = c_up_dg_i * c_up_i;
            auto n_dn_i    = c_dn_dg_i * c_dn_i;

            // hopping to neighbor (x+1, y)
            if (bc[0] == "pbc" || (bc[0] == "obc" && x < Lx - 1)) {
                lattice.coor2site({x+1,y}, 0, site_j, work);
                auto c_up_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_up);
                auto c_dn_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_dn);
                auto c_up_dg_j = c_up_j; c_up_dg_j.dagger();
                auto c_dn_dg_j = c_dn_j; c_dn_dg_j.dagger();
                Hubbard.add_Ham(std::complex<double>(-t*cx,-t*sx) * ( c_up_dg_i * c_up_j ));
                Hubbard.add_Ham(std::complex<double>(-t*cx,t*sx ) * ( c_up_dg_j * c_up_i ));
                Hubbard.add_Ham(std::complex<double>(-t*cx,-t*sx) * ( c_dn_dg_i * c_dn_j ));
                Hubbard.add_Ham(std::complex<double>(-t*cx,t*sx ) * ( c_dn_dg_j * c_dn_i ));
            }

            // hopping to neighbor (x, y+1)
            if (bc[1] == "pbc" || (bc[1] == "obc" && y < Ly - 1)) {
                lattice.coor2site({x,y+1}, 0, site_j, work);
                auto c_up_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_up);
                auto c_dn_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_dn);
                auto c_up_dg_j = c_up_j; c_up_dg_j.dagger();
                auto c_dn_dg_j = c_dn_j; c_dn_dg_j.dagger();
                Hubbard.add_Ham(std::complex<double>(-t*cy,-t*sy) * ( c_up_dg_i * c_up_j ));
                Hubbard.add_Ham(std::complex<double>(-t*cy,t*sy ) * ( c_up_dg_j * c_up_i ));
                Hubbard.add_Ham(std::complex<double>(-t*cy,-t*sy) * ( c_dn_dg_i * c_dn_j ));
                Hubbard.add_Ham(std::complex<double>(-t*cy,t*sy ) * ( c_dn_dg_j * c_dn_i ));
            }

            // hopping to second neighbor (x+1, y+1)
            if (bc[0] == "pbc" || (bc[0] == "obc" && x < Lx - 1 && y < Ly - 1)) {
                lattice.coor2site({x+1,y+1}, 0, site_j, work);
                auto c_up_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_up);
                auto c_dn_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_dn);
                auto c_up_dg_j = c_up_j; c_up_dg_j.dagger();
                auto c_dn_dg_j = c_dn_j; c_dn_dg_j.dagger();
                Hubbard.add_Ham(std::complex<double>(-tp*cxpy,-tp*sxpy) * ( c_up_dg_i * c_up_j ));
                Hubbard.add_Ham(std::complex<double>(-tp*cxpy,tp*sxpy ) * ( c_up_dg_j * c_up_i ));
                Hubbard.add_Ham(std::complex<double>(-tp*cxpy,-tp*sxpy) * ( c_dn_dg_i * c_dn_j ));
                Hubbard.add_Ham(std::complex<double>(-tp*cxpy,tp*sxpy ) * ( c_dn_dg_j * c_dn_i ));
            }

            if (bc[0] == "pbc" || (bc[0] == "obc" && x < Lx - 1 && y < Ly - 1 && y > 0)) {
                lattice.coor2site({x+1,y-1}, 0, site_j, work);
                auto c_up_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_up);
                auto c_dn_j    = qbasis::opr<std::complex<double>>(site_j,0,true,c_dn);
                auto c_up_dg_j = c_up_j; c_up_dg_j.dagger();
                auto c_dn_dg_j = c_dn_j; c_dn_dg_j.dagger();
                Hubbard.add_Ham(std::complex<double>(-tp*cxmy,-tp*sxmy) * ( c_up_dg_i * c_up_j ));
                Hubbard.add_Ham(std::complex<double>(-tp*cxmy,tp*sxmy ) * ( c_up_dg_j * c_up_i ));
                Hubbard.add_Ham(std::complex<double>(-tp*cxmy,-tp*sxmy) * ( c_dn_dg_i * c_dn_j ));
                Hubbard.add_Ham(std::complex<double>(-tp*cxmy,tp*sxmy ) * ( c_dn_dg_j * c_dn_i ));
            }

            // Hubbard repulsion, note that this operator is a sum (over sites) of diagonal matrices
            Hubbard.add_Ham(std::complex<double>(U,0.0) * (n_up_i * n_dn_i));

            // total electron operator
            Nup   += n_up_i;
            Ndown += n_dn_i;
        }
    }


    // to use translational symmetry, we first fill the Weisse tables
    Hubbard.fill_Weisse_table();


    // measure operators
    //qbasis::mopr<std::complex<double>> cupdg1cup5(qbasis::opr<std::complex<double>>(1,0,true,c_up).dagger() *
    //                  qbasis::opr<std::complex<double>>(5,0,true,c_up));
    //std::complex<double> m1;

    std::vector<double> E0_list;
    for (int m = 0; m < 1; m++) {
        for (int n = 0; n < 1; n++) {
            
            // constructing the Hilbert space basis
            Hubbard.enumerate_basis_repr({Qm,Qn}, {Nup,Ndown}, {Nup_total,Ndn_total});

            // generating matrix of the Hamiltonian in the subspace
            Hubbard.generate_Ham_sparse_repr();
            std::cout << std::endl;

            // obtaining the lowest eigenvals of the matrix
            Hubbard.locate_E0_lanczos(qbasis::which_sym::repr);
            std::cout << std::endl;

            E0_list.push_back(Hubbard.eigenvals_repr[0]);

            //if (m == 0 && n == 0) {
            //    m1 = Hubbard.measure_repr_static(cupdg1cup5, 0, 0);
            //}
        }
    }

    // prepare list of momentum points for measurement
    std::vector<std::vector<int>> q_list;
    // along q in the upper triangle
    //for (int m = 0; m <= Lx/2; m++) {
    //    for (int n = 0; n <= Ly/2; n++) {
    for (int m = 0; m < 1; m++) {
        for (int n = 0; n < 1; n++) {
            q_list.push_back({m,n});
        }
    }
    //std::cout << "q_list = " << q_list << std::endl;

    // calculate single-particle Green function G(q,w)
    std::string output_name_p = "L"+std::to_string(Lx)+"x"+std::to_string(Ly)+"Nup"+std::to_string(int(Nup_total))+"Ndn"+std::to_string(int(Ndn_total))+"f"+std::to_string(kx)+"x"+std::to_string(ky)+"Q"+std::to_string(Qm)+"x"+std::to_string(Qn)+"_Gp.dat";
    std::ofstream fout_p(output_name_p, std::ios::out);
    fout_p << std::scientific << std::right;
    fout_p << std::setprecision(18);
    fout_p << "Lx\t" << Lx << std::endl;
    fout_p << "Ly\t" << Ly << std::endl;
    fout_p << "U\t" << U << std::endl;
    fout_p << "E0\t" << Hubbard.energy_min() << std::endl;
    fout_p << "Gap\t" << Hubbard.energy_gap() << std::endl << std::endl;

    std::string output_name_h = "L"+std::to_string(Lx)+"x"+std::to_string(Ly)+"Nup"+std::to_string(int(Nup_total))+"Ndn"+std::to_string(int(Ndn_total))+"f"+std::to_string(kx)+"x"+std::to_string(ky)+"Q"+std::to_string(Qm)+"x"+std::to_string(Qn)+"_Gh.dat";
    std::ofstream fout_h(output_name_h, std::ios::out);
    fout_h << std::scientific << std::right;
    fout_h << std::setprecision(18);
    fout_h << "Lx\t" << Lx << std::endl;
    fout_h << "Ly\t" << Ly << std::endl;
    fout_h << "U\t" << U << std::endl;
    fout_h << "E0\t" << Hubbard.energy_min() << std::endl;
    fout_h << "Gap\t" << Hubbard.energy_gap() << std::endl << std::endl;

    const MKL_INT maxit=1000;
    MKL_INT step_used_p, step_used_h;
    double norm_p, norm_h;
    //std::vector<double> Hessenberg(2*maxit),Hessenberg_h(2*maxit), Hessenberg_p(2*maxit);

    // loop over all the q points
    for (auto q : q_list) {
        int m = q[0];
        int n = q[1];
        std::cout << "k\t" << m << "," << n << std::endl;
        std::vector<double> Hessenberg(2*maxit),Hessenberg_h(2*maxit), Hessenberg_p(2*maxit);

        // construct c_up(q) and c_up^\dagger(q)
        qbasis::mopr<std::complex<double>> c_up_q, c_up_dg_q;
        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                double qdotr = 2.0 * PI * (m * x / static_cast<double>(Lx) + n * y / static_cast<double>(Ly));
                auto coeff = 1.0 / sqrt(static_cast<double>(Lx * Ly)) * std::exp(std::complex<double>{0.0,qdotr});
                auto coeffm = 1.0 / sqrt(static_cast<double>(Lx * Ly)) * std::exp(std::complex<double>{0.0,-qdotr});
                uint32_t site;
                std::vector<int> work(lattice.dim);
                lattice.coor2site({x,y}, 0, site, work);
                auto c_up_i    = qbasis::opr<std::complex<double>>(site,0,true,c_up);
                auto c_up_dg_i = c_up_i; c_up_dg_i.dagger();
                c_up_q += coeff * c_up_i;
                c_up_dg_q += coeffm * c_up_dg_i;
            }
        }


        // calculate Gp
        Hubbard.enumerate_basis_repr({Qm+m,Qn+n}, {Nup,Ndown}, {Nup_total+1,Ndn_total},1);
        Hubbard.generate_Ham_sparse_repr(1);
        Hubbard.switch_sec_mat(1);
        
        Hubbard.measure_repr_dynamic(c_up_dg_q, 0, 1, maxit, step_used_p, norm_p, Hessenberg_p.data());
        
        std::cout << "norm_p = " << norm_p << std::endl;
        fout_p << "Q\t" << m << "," << n << std::endl;
        fout_p << "nrm2\t" << norm_p << std::endl;
        if (norm_p > qbasis::lanczos_precision) {
            fout_p << "b\t";
            for (MKL_INT i = 0; i < step_used_p; i++) {
                fout_p << std::setw(30) << std::setprecision(18) << Hessenberg_p[i];
            }
            fout_p << std::endl;
            fout_p << "a\t";
            for (MKL_INT i = maxit; i < maxit + step_used_p; i++) {
                fout_p << std::setw(30) << std::setprecision(18) << Hessenberg_p[i];
            }
            fout_p << std::endl;
        } else {
            fout_p << "b\t" << std::endl;
            fout_p << "a\t" << std::endl;
        }
        std::cout << std::endl;

        // calculate Gh
        Hubbard.enumerate_basis_repr({Qm-m,Qn-n}, {Nup,Ndown}, {Nup_total-1,Ndn_total},2);
        Hubbard.generate_Ham_sparse_repr(2);
        Hubbard.switch_sec_mat(2);
        
        Hubbard.measure_repr_dynamic(c_up_q, 0, 2, maxit, step_used_h, norm_h, Hessenberg_h.data());

        std::cout << "norm_h = " << norm_h << std::endl;
        fout_h << "Q\t" << m << "," << n << std::endl;
        fout_h << "nrm2\t" << norm_h << std::endl;
        if (norm_h > qbasis::lanczos_precision) {
            fout_h << "b\t";
            for (MKL_INT i = 0; i < step_used_h; i++) {
                fout_h << std::setw(30) << std::setprecision(18) << Hessenberg_h[i];
            }
            fout_h << std::endl;
            fout_h << "a\t";
            for (MKL_INT i = maxit; i < maxit + step_used_h; i++) {
                fout_h << std::setw(30) << std::setprecision(18) << Hessenberg_h[i];
            }
            fout_h << std::endl;
        } else {
            fout_h << "b\t" << std::endl;
            fout_h << "a\t" << std::endl;
        }
        std::cout << std::endl;
        
    }
    fout_p.close();
    fout_h.close();

}
