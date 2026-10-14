#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "qbasis.h"

#define PI 3.1415926535897932

// Fermi-Hubbard model on a 1D chain with twisted boundary condition
int main(int argc, char *argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " N_up N_dn <k> Q" << std::endl;
        return 1;
    }

    int Nupinput = std::atoi(argv[1]);
    int Ndninput = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);
    int Q = std::atoi(argv[4]);

    qbasis::initialize();
    std::cout << std::setprecision(10);

    double t = 1.0;
    double U = 8.0;
    int L = 4;
    double Nup_total = 1.0 * Nupinput;
    double Ndn_total = 1.0 * Ndninput;

    double f = k * PI / 40.0 + 0.00001;
    double sf = std::sin(f);
    double cf = std::cos(f);

    std::cout << "L =       " << L << std::endl;
    std::cout << "t =       " << t << std::endl;
    std::cout << "U =       " << U << std::endl;
    std::cout << "N_up =    " << Nup_total << std::endl;
    std::cout << "N_dn =    " << Ndn_total << std::endl << std::endl;

    std::vector<std::string> bc{"pbc"};
    qbasis::lattice lattice("chain", {static_cast<uint32_t>(L)}, bc);

    auto c_up = std::vector<std::vector<std::complex<double>>>(4, std::vector<std::complex<double>>(4, 0.0));
    auto c_dn = std::vector<std::vector<std::complex<double>>>(4, std::vector<std::complex<double>>(4, 0.0));
    c_up[0][1] = std::complex<double>(1.0, 0.0);
    c_up[2][3] = std::complex<double>(1.0, 0.0);
    c_dn[0][2] = std::complex<double>(1.0, 0.0);
    c_dn[1][3] = std::complex<double>(-1.0, 0.0);

    qbasis::model<std::complex<double>> Hubbard(lattice);
    Hubbard.add_orbital(lattice.Nsites, "electron");
    qbasis::mopr<std::complex<double>> Nup;
    qbasis::mopr<std::complex<double>> Ndown;

    for (int x = 0; x < L; x++) {
        uint32_t site_i, site_j;
        std::vector<int> work(lattice.dim);
        lattice.coor2site({x}, 0, site_i, work);

        auto c_up_i = qbasis::opr<std::complex<double>>(site_i, 0, true, c_up);
        auto c_dn_i = qbasis::opr<std::complex<double>>(site_i, 0, true, c_dn);
        auto c_up_dg_i = c_up_i;
        auto c_dn_dg_i = c_dn_i;
        c_up_dg_i.dagger();
        c_dn_dg_i.dagger();
        auto n_up_i = c_up_dg_i * c_up_i;
        auto n_dn_i = c_dn_dg_i * c_dn_i;

        if (bc[0] == "pbc" || (bc[0] == "obc" && x < L - 1)) {
            lattice.coor2site({x + 1}, 0, site_j, work);
            auto c_up_j = qbasis::opr<std::complex<double>>(site_j, 0, true, c_up);
            auto c_dn_j = qbasis::opr<std::complex<double>>(site_j, 0, true, c_dn);
            auto c_up_dg_j = c_up_j;
            auto c_dn_dg_j = c_dn_j;
            c_up_dg_j.dagger();
            c_dn_dg_j.dagger();

            Hubbard.add_Ham(std::complex<double>(-t * cf, -t * sf) * (c_up_dg_i * c_up_j));
            Hubbard.add_Ham(std::complex<double>(-t * cf, t * sf) * (c_up_dg_j * c_up_i));
            Hubbard.add_Ham(std::complex<double>(-t * cf, -t * sf) * (c_dn_dg_i * c_dn_j));
            Hubbard.add_Ham(std::complex<double>(-t * cf, t * sf) * (c_dn_dg_j * c_dn_i));
        }

        Hubbard.add_Ham(std::complex<double>(U, 0.0) * (n_up_i * n_dn_i));

        Nup += n_up_i;
        Ndown += n_dn_i;
    }

    Hubbard.fill_Weisse_table();

    Hubbard.enumerate_basis_repr({Q}, {Nup, Ndown}, {Nup_total, Ndn_total});
    Hubbard.generate_Ham_sparse_repr();
    std::cout << std::endl;

    Hubbard.locate_E0_lanczos(qbasis::which_sym::repr);
    std::cout << std::endl;

    std::vector<int> q_list;
    for (int momentum = 0; momentum < 1; momentum++) {
        q_list.push_back(momentum);
    }

    std::string output_name_p =
        "L" + std::to_string(L) + "Nup" + std::to_string(int(Nup_total)) +
        "Ndn" + std::to_string(int(Ndn_total)) + "f" + std::to_string(k) +
        "Q" + std::to_string(Q) + "_Gp.dat";
    std::ofstream fout_p(output_name_p, std::ios::out);
    fout_p << std::scientific << std::right;
    fout_p << std::setprecision(18);
    fout_p << "L\t" << L << std::endl;
    fout_p << "U\t" << U << std::endl;
    fout_p << "E0\t" << Hubbard.energy_min() << std::endl;
    fout_p << "Gap\t" << Hubbard.energy_gap() << std::endl << std::endl;

    std::string output_name_h =
        "L" + std::to_string(L) + "Nup" + std::to_string(int(Nup_total)) +
        "Ndn" + std::to_string(int(Ndn_total)) + "f" + std::to_string(k) +
        "Q" + std::to_string(Q) + "_Gh.dat";
    std::ofstream fout_h(output_name_h, std::ios::out);
    fout_h << std::scientific << std::right;
    fout_h << std::setprecision(18);
    fout_h << "L\t" << L << std::endl;
    fout_h << "U\t" << U << std::endl;
    fout_h << "E0\t" << Hubbard.energy_min() << std::endl;
    fout_h << "Gap\t" << Hubbard.energy_gap() << std::endl << std::endl;

    const MKL_INT maxit = 1000;
    MKL_INT step_used_p, step_used_h;
    double norm_p, norm_h;

    for (int q : q_list) {
        int Qp = (Q + q) % L;
        int Qh = (Q - q) % L;
        if (Qh < 0) {
            Qh += L;
        }
        std::cout << "k\t" << q << std::endl;
        std::vector<double> Hessenberg_h(2 * maxit), Hessenberg_p(2 * maxit);

        qbasis::mopr<std::complex<double>> c_up_q, c_up_dg_q;
        for (int x = 0; x < L; x++) {
            double qdotr = 2.0 * PI * q * x / static_cast<double>(L);
            auto coeff = 1.0 / std::sqrt(static_cast<double>(L)) *
                         std::exp(std::complex<double>{0.0, qdotr});
            auto coeffm = 1.0 / std::sqrt(static_cast<double>(L)) *
                          std::exp(std::complex<double>{0.0, -qdotr});
            uint32_t site;
            std::vector<int> work(lattice.dim);
            lattice.coor2site({x}, 0, site, work);
            auto c_up_i = qbasis::opr<std::complex<double>>(site, 0, true, c_up);
            auto c_up_dg_i = c_up_i;
            c_up_dg_i.dagger();
            c_up_q += coeff * c_up_i;
            c_up_dg_q += coeffm * c_up_dg_i;
        }

        Hubbard.enumerate_basis_repr({Qp}, {Nup, Ndown}, {Nup_total + 1, Ndn_total}, 1);
        Hubbard.generate_Ham_sparse_repr(1);
        Hubbard.switch_sec_mat(1);

        Hubbard.measure_repr_dynamic(c_up_dg_q, 0, 1, maxit, step_used_p, norm_p, Hessenberg_p.data());

        std::cout << "norm_p = " << norm_p << std::endl;
        fout_p << "Q\t" << q << std::endl;
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

        Hubbard.enumerate_basis_repr({Qh}, {Nup, Ndown}, {Nup_total - 1, Ndn_total}, 2);
        Hubbard.generate_Ham_sparse_repr(2);
        Hubbard.switch_sec_mat(2);

        Hubbard.measure_repr_dynamic(c_up_q, 0, 2, maxit, step_used_h, norm_h, Hessenberg_h.data());

        std::cout << "norm_h = " << norm_h << std::endl;
        fout_h << "Q\t" << q << std::endl;
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
