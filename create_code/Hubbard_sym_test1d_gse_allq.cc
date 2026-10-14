#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "qbasis.h"

#define PI 3.1415926535897932

// Fermi-Hubbard model on a 1D chain
int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " N_up N_dn <k>" << std::endl;
        return 1;
    }

    int Nupinput = std::atoi(argv[1]);
    int Ndninput = std::atoi(argv[2]);
    int k = std::atoi(argv[3]);

    qbasis::initialize();
    std::cout << std::setprecision(10);

    double t = 1.0;
    double U = 6.0;
    int L = 16;
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

    std::string output_name =
        "L" + std::to_string(L) + "Nup" + std::to_string(int(Nup_total)) +
        "Ndn" + std::to_string(int(Ndn_total)) + "f" + std::to_string(k) +
        "_gse.dat";
    std::ofstream fout(output_name, std::ios::out);
    fout << std::scientific << std::right;
    fout << std::setprecision(18);
    fout << "L\t" << L << std::endl;
    fout << "U\t" << U << std::endl;

    for (int Q = 0; Q < L; Q += L / 2) {
        Hubbard.enumerate_basis_repr({Q}, {Nup, Ndown}, {Nup_total, Ndn_total});
        Hubbard.generate_Ham_sparse_repr();
        std::cout << std::endl;

        Hubbard.locate_E0_lanczos(qbasis::which_sym::repr, 1, 0);
        std::cout << std::endl;

        fout << "Q\t" << Q << std::endl;
        fout << "E0\t" << Hubbard.energy_min() << std::endl;
        fout << "Gap\t" << Hubbard.energy_gap() << std::endl << std::endl;
    }

    fout.close();
}
