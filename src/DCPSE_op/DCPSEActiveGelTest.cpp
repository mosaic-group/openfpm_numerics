//
// Created by Abhinav Singh on 23.04.20.
//
#include "config.h"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "DCPSE_op.hpp"
#include "DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "EqnsStruct.hpp"

template<typename particle_type, typename particle_type2>
void indexUpdate(
        particle_type &Part,
        particle_type2 &Part_subset,
        openfpm::vector<aggregate<int>> &up_p, openfpm::vector<aggregate<int>> &dw_p,
        openfpm::vector<aggregate<int>> &l_p, openfpm::vector<aggregate<int>> &r_p,
        openfpm::vector<aggregate<int>> &bulk, Box<2, double> &up, Box<2, double> &down, Box<2, double> &left,
        Box<2, double> &right) {
    up_p.clear();
    dw_p.clear();
    l_p.clear();
    r_p.clear();
    bulk.clear();
    Part_subset.clear();

    auto it2 = Part.getDomainIterator();
    while (it2.isNext()) {
        auto p = it2.get();
        Point<2, double> xp = Part.getPos(p);
        if (up.isInside(xp) == true) {
            up_p.add();
            up_p.last().get<0>() = p.getKey();
        } else if (down.isInside(xp) == true) {
            dw_p.add();
            dw_p.last().get<0>() = p.getKey();
        } else if (left.isInside(xp) == true) {
            l_p.add();
            l_p.last().get<0>() = p.getKey();
        } else if (right.isInside(xp) == true) {
            r_p.add();
            r_p.last().get<0>() = p.getKey();
        } else {
            bulk.add();
            bulk.last().get<0>() = p.getKey();
        }
        ++it2;
    }

    for (int i = 0; i < bulk.size(); i++) {
        Part_subset.add();
        Part_subset.getLastPos()[0] = Part.getPos(bulk.template get<0>(i))[0];
        Part_subset.getLastPos()[1] = Part.getPos(bulk.template get<0>(i))[1];
    }

}

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)

    BOOST_AUTO_TEST_CASE(Active2DPetsc) {
        timer tt2;
        tt2.start();
        const size_t sz[2] = {17, 17};
        Box<2, double> box({0, 0}, {10, 10});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.6 * spacing;
        double rCut2 = 3.6 * spacing;
        int ord = 2;
        double sampling_factor = 2.5;
        int ord2 = 2;
        double sampling_factor2 = 2.5;
        Ghost<2, double> ghost(rCut);
/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div   H_t                                                                                      delmu */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2], double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double,double,double>> Particles(
                0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkP;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        constexpr int x = 0;
        constexpr int y = 1;

        constexpr int Polarization = 0;
        constexpr int Velocity = 1;
        constexpr int Vorticity = 2;
        constexpr int ExtForce = 3;
        constexpr int Pressure = 4;
        constexpr int Strain_rate = 5;
        constexpr int Stress = 6;
        constexpr int MolField = 7;
        auto Pos = getV<PROP_POS>(Particles);

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dPol = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);
        auto Df1 = getV<21>(Particles);
        auto Df2 = getV<22>(Particles);
        auto Df3 = getV<23>(Particles);
        auto Df4 = getV<24>(Particles);
        auto Df5 = getV<25>(Particles);
        auto Df6 = getV<26>(Particles);
        auto delmu = getV<27>(Particles);
        auto k1 = getV<28>(Particles);
        auto k2 = getV<29>(Particles);
        auto k3 = getV<30>(Particles);
        auto k4 = getV<31>(Particles);
        auto Hperp_boundary = getV<32>(Particles);
        auto FranckEnergyDensity = getV<33>(Particles);
        auto r = getV<34>(Particles);


        double eta = 1.0;
        double nu = -0.5;
        double gama = 0.1;
        double zeta = 0.07;
        double Ks = 1.0;
        double Kb = 1.0;
        double lambda = 0.1;
        //double delmu = -1.0;
        g = 0;
        delmu = -1.0;
        P = 0;
        V = 0;
        // Here fill up the boxes for particle boundary detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        /*Box<2, double> mid({box.getHigh(0) / 2.0 - spacing, box.getHigh(1) / 2.0 - spacing / 2.0},
                           {box.getHigh(0) / 2.0, box.getHigh(1) / 2.0 + spacing / 2.0});
*/
        Box<2, double> mid({box.getLow(0) + 3.1 * spacing, box.getLow(1) + 3.1 * spacing},
                           {box.getHigh(0) - 3.1 * spacing, box.getHigh(1) - 3.1 * spacing});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            } else {
                if (mid.isInside(xp) == true) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                    Particles.getProp<4>(p) = 0;
                } else {
                    bulkP.add();
                    bulkP.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }


            ++it2;
        }

        Derivative_x Dx(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_y Dy(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);

        Hperp_boundary = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                          Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));


        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        //solverPetsc.setPreconditioner(PCJACOBI);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setSolver(KSPGMRES);
        //solverPetsc2.setRestart(250);
        /*solverPetsc2.setPreconditioner(PCJACOBI);
        solverPetsc2.setRelTol(1e-6);
        solverPetsc2.setAbsTol(1e-5);
        solverPetsc2.setDivTol(1e4);*/

        eq_id vx, vy;
        timer tt;
        timer tt3;
        vx.setId(0);
        vy.setId(1);
        double V_err_eps = 1e-4;
        double V_err = 1;
        int n = 0;
        int ctr = 0;
        double dt = 1e-3;
        double tim = 0;
        double tf = 1e-1;
        while (tim <= tf) {
            tt.start();
            Particles.ghost_get<Polarization>();
            sigma[x][x] =
                    -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
            sigma[x][y] =
                    -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
            sigma[y][x] =
                    -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
            sigma[y][y] =
                    -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
            Particles.ghost_get<Stress>();


            h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));
/*
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<7>(p)[y] = Particles.getProp<32>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<7>(p)[y] = Particles.getProp<32>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<7>(p)[y] = Particles.getProp<32>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<7>(p)[y] = Particles.getProp<32>(p);
            }*/

            Particles.ghost_get<MolField>();

            f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Particles.ghost_get<11, 12, 13, 14, 15, 16>();
            Df1[x] = Dx(f1);
            Df2[x] = Dx(f2);
            Df3[x] = Dx(f3);
            Df4[x] = Dx(f4);
            Df5[x] = Dx(f5);
            Df6[x] = Dx(f6);

            Df1[y] = Dy(f1);
            Df2[y] = Dy(f2);
            Df3[y] = Dy(f3);
            Df4[y] = Dy(f4);
            Df5[y] = Dy(f5);
            Df6[y] = Dy(f6);
            Particles.ghost_get<21, 22, 23, 24, 25, 26>();


            dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                    zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) -
                    g[x]
                    - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));


            dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                    zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) -
                    g[y]
                    - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
            Particles.ghost_get<9>();


            Particles.write("PolarI");
            //Velocity Solution n iterations
            double sum, sum1;

            auto Stokes1 = eta * (Dxx(V[x]) + Dyy(V[x]))
                           + 0.5 * nu * (Df1[x] * Dx(V[x]) + f1 * Dxx(V[x]))
                           + 0.5 * nu * (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y]))
                           + 0.5 * nu * (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x]))
                           + 0.5 * nu * (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           + 0.5 * nu * (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
            auto Stokes2 = eta * (Dxx(V[y]) + Dyy(V[y]))
                           - 0.5 * nu * (Df1[y] * Dx(V[x]) + f1 * Dxy(V[x]))
                           - 0.5 * nu * (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           - 0.5 * nu * (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y]))
                           + 0.5 * nu * (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x]))
                           + 0.5 * nu * (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));
            tt.stop();
            std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
            tt.start();
            V_err = 1;
            n = 0;
            while (V_err >= V_err_eps && n <= 30) {
                RHS[x] = Dx(P) + dV[x];
                RHS[y] = Dy(P) + dV[y];
                Particles.ghost_get<10>();
                DCPSE_scheme<equations2d2, decltype(Particles)> Solver(Particles);
                Solver.impose(Stokes1, bulk, RHS[0], vx);
                Solver.impose(Stokes2, bulk, RHS[1], vy);
                Solver.impose(V[x], up_p, 0, vx);
                Solver.impose(V[y], up_p, 0, vy);
                Solver.impose(V[x], dw_p, 0, vx);
                Solver.impose(V[y], dw_p, 0, vy);
                Solver.impose(V[x], l_p, 0, vx);
                Solver.impose(V[y], l_p, 0, vy);
                Solver.impose(V[x], r_p, 0, vx);
                Solver.impose(V[y], r_p, 0, vy);
                Solver.solve_with_solver(solverPetsc, V[x], V[y]);

                //std::cout << "Stokes Solved in " << tt.getwct() << " seconds." << std::endl;
                Particles.ghost_get<Velocity>();
                div = -(Dx(V[x]) + Dy(V[y]));
                Particles.ghost_get<19>();
                auto Helmholtz = Dxx(H) + Dyy(H);
                DCPSE_scheme<equations2d1, decltype(Particles)> SolverH(Particles);
                SolverH.impose(Helmholtz, bulk, prop_id<19>());
                SolverH.impose(H, up_p, 0);
                SolverH.impose(H, dw_p, 0);
                SolverH.impose(H, l_p, 0);
                SolverH.impose(H, r_p, 0);
                SolverH.solve_with_solver(solverPetsc2, H);
                //std::cout << "Helmholtz Solved in " << tt.getwct() << " seconds." << std::endl;
                Particles.ghost_get<17>();
                Particles.ghost_get<Velocity>();
                V[x] = V[x] + Dx(H);
                V[y] = V[y] + Dy(H);
                P = P + (Dxx(H) + Dyy(H));
                for (int j = 0; j < up_p.size(); j++) {
                    auto p = up_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    //Particles.getProp<4>(p) = 0;

                }
                for (int j = 0; j < dw_p.size(); j++) {
                    auto p = dw_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    // Particles.getProp<4>(p) = 0;
                }
                for (int j = 0; j < l_p.size(); j++) {
                    auto p = l_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    //Particles.getProp<4>(p) = 0;
                }
                for (int j = 0; j < r_p.size(); j++) {
                    auto p = r_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    // Particles.getProp<4>(p) = 0;
                }
                Particles.ghost_get<Velocity>();
                Particles.ghost_get<Pressure>();
                sum = 0;
                sum1 = 0;
                for (int j = 0; j < bulk.size(); j++) {
                    auto p = bulk.get<0>(j);
                    sum += (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) *
                           (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) +
                           (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]) *
                           (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]);
                    sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                            Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1];
                }
                sum = sqrt(sum);
                sum1 = sqrt(sum1);
                V_t = V;
                V_err = sum / sum1;
                n++;
            }
            tt.stop();
            std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds with " << n
                      << " iterations."
                      << std::endl;
            u[x][x] = Dx(V[x]);
            u[x][y] = 0.5 * (Dx(V[y]) + Dy(V[x]));
            u[y][x] = 0.5 * (Dy(V[x]) + Dx(V[y]));
            u[y][y] = Dy(V[y]);

            W[x][x] = 0;
            W[x][y] = 0.5 * (Dy(V[x]) - Dx(V[y]));
            W[y][x] = 0.5 * (Dx(V[y]) - Dy(V[x]));
            W[y][y] = 0;


            Particles.write_frame("Polar_Petsc", ctr);
            ctr++;

            k1[x] = ((h[x] * Pol[x] - h[y] * Pol[y]) / gama + lambda * delmu * Pol[x] -
                     nu * (u[x][x] * Pol[x] + u[x][y] * Pol[y]) + W[x][x] * Pol[x] +
                     W[x][y] * Pol[y]);// - V[x] * Dx(Pol[x]) - V[y] * Dy(Pol[x]));
            k1[y] = ((h[x] * Pol[y] - h[y] * Pol[x]) / gama + lambda * delmu * Pol[y] -
                     nu * (u[y][x] * Pol[x] + u[y][y] * Pol[y]) + W[y][x] * Pol[x] +
                     W[y][y] * Pol[y]);// - V[x] * Dx(Pol[y]) - V[y] * Dy(Pol[y]));

            //dPol = Pol + (dt / 6.0) * k1;
            dPol = Pol + (0.5 * dt) * k1;

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum = sqrt(Particles.getProp<8>(p)[x] * Particles.getProp<8>(p)[x] +
                           Particles.getProp<8>(p)[y] * Particles.getProp<8>(p)[y]);
                Particles.getProp<0>(p)[x] /= sum;
                Particles.getProp<0>(p)[y] /= sum;
            }

            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu
                            - nu * u[x][x] * dPol[x] * dPol[x] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
                            - nu * u[y][y] * dPol[y] * dPol[y] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
                            - (2 * nu) * u[x][y] * dPol[x] * dPol[y] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
            );


            k2[x] = ((h[x] * (0.5 * dt * k1[x] + Pol[x]) - h[y] * (0.5 * dt * k1[y] + Pol[y])) / gama +
                     lambda * delmu * (0.5 * dt * k1[x] + Pol[x]) -
                     nu * (u[x][x] * (0.5 * dt * k1[x] + Pol[x]) + u[x][y] * (0.5 * dt * k1[y] + Pol[y])) +
                     W[x][x] * (0.5 * dt * k1[x] + Pol[x]) + W[x][y] * (0.5 * dt * k1[y] + Pol[y]));
            //-V[x] * Dx((0.5 * dt * k1[x] + Pol[x])) - V[y] * Dy((0.5 * dt * k1[x] + Pol[x])));
            k2[y] = ((h[x] * (0.5 * dt * k1[y] + Pol[y]) - h[y] * (0.5 * dt * k1[x] + Pol[x])) / gama +
                     lambda * delmu * (0.5 * dt * k1[y] + Pol[y]) -
                     nu * (u[y][x] * (0.5 * dt * k1[x] + Pol[x]) + u[y][y] * (0.5 * dt * k1[y] + Pol[y])) +
                     W[y][x] * (0.5 * dt * k1[x] + Pol[x]) + W[y][y] * (0.5 * dt * k1[y] + Pol[y]));
            //-V[x] * Dx((0.5 * dt * k1[y] + Pol[y])) - V[y] * Dy((0.5 * dt * k1[y] + Pol[y])));

            dPol = Pol + (0.5 * dt) * k2;

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum = sqrt(Particles.getProp<8>(p)[x] * Particles.getProp<8>(p)[x] +
                           Particles.getProp<8>(p)[y] * Particles.getProp<8>(p)[y]);
                Particles.getProp<0>(p)[x] /= sum;
                Particles.getProp<0>(p)[y] /= sum;
            }


            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu
                            - nu * u[x][x] * dPol[x] * dPol[x] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
                            - nu * u[y][y] * dPol[y] * dPol[y] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
                            - (2 * nu) * u[x][y] * dPol[x] * dPol[y] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
            );


            k3[x] = ((h[x] * (0.5 * dt * k2[x] + Pol[x]) - h[y] * (0.5 * dt * k2[y] + Pol[y])) / gama +
                     lambda * delmu * (0.5 * dt * k2[x] + Pol[x]) -
                     nu * (u[x][x] * (0.5 * dt * k2[x] + Pol[x]) + u[x][y] * (0.5 * dt * k2[y] + Pol[y])) +
                     W[x][x] * (0.5 * dt * k2[x] + Pol[x]) + W[x][y] * (0.5 * dt * k2[y] + Pol[y]));
            // -V[x] * Dx((0.5 * dt * k2[x] + Pol[x])) - V[y] * Dy((0.5 * dt * k2[x] + Pol[x])));
            k3[y] = ((h[x] * (0.5 * dt * k2[y] + Pol[y]) - h[y] * (0.5 * dt * k2[x] + Pol[x])) / gama +
                     lambda * delmu * (0.5 * dt * k2[y] + Pol[y]) -
                     nu * (u[y][x] * (0.5 * dt * k2[x] + Pol[x]) + u[y][y] * (0.5 * dt * k2[y] + Pol[y])) +
                     W[y][x] * (0.5 * dt * k2[x] + Pol[x]) + W[y][y] * (0.5 * dt * k2[y] + Pol[y]));
            // -V[x] * Dx((0.5 * dt * k2[y] + Pol[y])) - V[y] * Dy((0.5 * dt * k2[y] + Pol[y])));

            dPol = Pol + (dt * k3);

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum = sqrt(Particles.getProp<8>(p)[x] * Particles.getProp<8>(p)[x] +
                           Particles.getProp<8>(p)[y] * Particles.getProp<8>(p)[y]);
                Particles.getProp<0>(p)[x] /= sum;
                Particles.getProp<0>(p)[y] /= sum;
            }

            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<8>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu
                            - nu * u[x][x] * dPol[x] * dPol[x] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
                            - nu * u[y][y] * dPol[y] * dPol[y] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
                            - (2 * nu) * u[x][y] * dPol[x] * dPol[y] / (dPol[x] * dPol[x] + dPol[y] * dPol[y])
            );


            k4[x] = ((h[x] * (dt * k3[x] + Pol[x]) - h[y] * (dt * k3[y] + Pol[y])) / gama +
                     lambda * delmu * (dt * k3[x] + Pol[x]) -
                     nu * (u[x][x] * (dt * k3[x] + Pol[x]) + u[x][y] * (dt * k3[y] + Pol[y])) +
                     W[x][x] * (dt * k3[x] + Pol[x]) +
                     W[x][y] *
                     (dt * k3[y] + Pol[y]));//   -V[x]*Dx( (dt * k3[x]+Pol[x])) -V[y]*Dy( (dt * k3[x]+Pol[x])));
            k4[y] = ((h[x] * (dt * k3[y] + Pol[y]) - h[y] * (dt * k3[x] + Pol[x])) / gama +
                     lambda * delmu * (dt * k3[y] + Pol[y]) -
                     nu * (u[y][x] * (dt * k3[x] + Pol[x]) + u[y][y] * (dt * k3[y] + Pol[y])) +
                     W[y][x] * (dt * k3[x] + Pol[x]) +
                     W[y][y] *
                     (dt * k3[y] + Pol[y]));//  -V[x]*Dx( (dt * k3[y]+Pol*[y])) -V[y]*Dy( (dt * k3[y]+Pol[y])));

            Pol = Pol + (dt / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }

            Pos = Pos + dt * V;
            Particles.map();
            Particles.ghost_get<0>();
            //indexUpdate(Particles, Particles_subset, up_p, dw_p, l_p, r_p, bulk, up, down, left, right);
            //Particles_subset.map();
            //Particles_subset.ghost_get<0>();

            //Particles_subset.write("debug");
            Particles.ghost_get<0, Vorticity, MolField>();
            //Particles_subset.ghost_get<0>();

            tt.start();
            Dx.update(Particles);
            Dy.update(Particles);
            Dxy.update(Particles);
            auto Dyx = Dxy;
            Dxx.update(Particles);
            Dyy.update(Particles);

            //Bulk_Dx.update(Particles_subset);
            //Bulk_Dy.update(Particles_subset);

            tt.stop();
            std::cout << "Updation of operators took " << tt.getwct() << " seconds." << std::endl;
            //Pol = Pol + (dt / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);

            /* for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum = sqrt(Particles.getProp<0>(p)[x] * Particles.getProp<0>(p)[x] +
                           Particles.getProp<0>(p)[y] * Particles.getProp<0>(p)[y]);
                Particles.getProp<0>(p)[x] /= sum;
                Particles.getProp<0>(p)[y] /= sum;
            }*/
            Particles.ghost_get<0>();


            std::cout << "Time step " << ctr << " : " << tim << " over." << std::endl;
            tim += dt;
            std::cout << "----------------------------------------------------------" << std::endl;
        }
        Particles.deleteGhost();
        tt2.stop();
        std::cout << "The simulation took " << tt2.getwct() << "Seconds.";
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    BOOST_AUTO_TEST_CASE(Active2DEigen) {
        timer tt2;
        tt2.start();
        const size_t sz[2] = {81, 81};
        Box<2, double> box({0, 0}, {10, 10});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        double rCut2 = 3.1 * spacing;
        int ord = 2;
        int ord2 = 2;
        double sampling_factor = 1.9;
        double sampling_factor2 = 1.9;
        double alpha_V = 1.0;
        double alpha_P = 1.0;
        Ghost<2, double> ghost(rCut);

        auto &v_cl = create_vcluster();

/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div    H_t                                                                                      delmu */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2], double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double,double,double>> Particles(
                0, box, bc, ghost);
        vector_dist<2, double, aggregate<double, double, VectorS<2, double>>> Particles_subset(
                Particles.getDecomposition(), 0);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        constexpr int x = 0;
        constexpr int y = 1;

        constexpr int Polarization = 0;
        constexpr int Velocity = 1;
        constexpr int Vorticity = 2;
        constexpr int ExtForce = 3;
        constexpr int Pressure = 4;
        constexpr int Strain_rate = 5;
        constexpr int Stress = 6;
        constexpr int MolField = 7;
        auto Pos = getV<PROP_POS>(Particles);

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto P_bulk = getV<0>(Particles_subset); //Pressure only on inside
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dPol = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto H_bulk = getV<1>(Particles_subset); //Pressure only on inside
        auto Grad_bulk = getV<2>(Particles_subset);

        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);
        auto Df1 = getV<21>(Particles);
        auto Df2 = getV<22>(Particles);
        auto Df3 = getV<23>(Particles);
        auto Df4 = getV<24>(Particles);
        auto Df5 = getV<25>(Particles);
        auto Df6 = getV<26>(Particles);
        auto delmu = getV<27>(Particles);
        auto k1 = getV<28>(Particles);
        auto k2 = getV<29>(Particles);
        auto k3 = getV<30>(Particles);
        auto k4 = getV<31>(Particles);
        auto H_p_b = getV<32>(Particles);
        auto FranckEnergyDensity = getV<33>(Particles);
        auto r = getV<34>(Particles);


        double eta = 1.0;
        double nu = -0.5;
        double gama = 0.1;
        double zeta = 0.07;
        double Ks = 1.0;
        double Kb = 1.0;
        double lambda = 0.1;
        //double delmu = -1.0;
        g = 0;
        delmu = -1.0;
        P = 0;
        P_bulk = 0;
        V = 0;
        // Here fill up the boxes for particle boundary detection.
        Particles.ghost_get<ExtForce,28>();


        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        /*Box<2, double> mid({box.getHigh(0) / 2.0 - spacing, box.getHigh(1) / 2.0 - spacing / 2.0},
                           {box.getHigh(0) / 2.0, box.getHigh(1) / 2.0 + spacing / 2.0});
*/
        Box<2, double> mid({box.getLow(0) + 3.1 * spacing, box.getLow(1) + 3.1 * spacing},
                           {box.getHigh(0) - 3.1 * spacing, box.getHigh(1) - 3.1 * spacing});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }


        for (int i = 0; i < bulk.size(); i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
        }


        Particles_subset.map();
        Particles_subset.ghost_get<0>();

        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dx(Particles_subset, ord2,
                                                                                                 rCut2,
                                                                                                 sampling_factor2,
                                                                                                 support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dy(Particles_subset, ord2,
                                                                                                 rCut2,
                                                                                                 sampling_factor2,
                                                                                                 support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut, sampling_factor,
                          support_options::RADIUS);//,Bulk_Dxx(Particles_subset, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut, sampling_factor,
                          support_options::RADIUS);//,Bulk_Dyy(Particles_subset, ord2, rCut2, sampling_factor2, support_options::RADIUS);

        /*Derivative_x Dx(Particles, ord, rCut, sampling_factor), Bulk_Dx(Particles_subset, ord, rCut, sampling_factor2);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor), Bulk_Dy(Particles_subset, ord, rCut, sampling_factor2);
        Derivative_xy Dxy(Particles, ord2, rCut2, sampling_factor2);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord2, rCut2, sampling_factor2),Bulk_Dxx(Particles_subset, ord2, rCut2, sampling_factor2);
        Derivative_yy Dyy(Particles, ord2, rCut2, sampling_factor2),Bulk_Dyy(Particles_subset, ord2, rCut2, sampling_factor2);*/

        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc.setPreconditioner(PCJACOBI);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setSolver(KSPGMRES);
        solverPetsc2.setPreconditioner(PCJACOBI);


        eq_id vx, vy;
        timer tt;
        timer tt3;
        vx.setId(0);
        vy.setId(1);
        double V_err_eps = 5e-4;
        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 300;
        int ctr = 0, errctr;
        double dt = 3e-3;
        double tim = 0;
        double tf = 1e-1;
        div = 0;
        double sum, sum1, sum_k;
        while (tim <= tf) {
            tt.start();
            Particles.ghost_get<Polarization>();
            sigma[x][x] =
                    -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
            sigma[x][y] =
                    -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
            sigma[y][x] =
                    -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
            sigma[y][y] =
                    -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
            Particles.ghost_get<Stress>();


            h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));

            Particles.ghost_get<MolField>();

            FranckEnergyDensity = (Ks / 2.0) * ((Dx(Pol[x]) * Dx(Pol[x])) + (Dy(Pol[x]) * Dy(Pol[x])) + (Dx(Pol[y]) * Dx(Pol[y])) + (Dy(Pol[y]) * Dy(Pol[y]))) + ((Kb - Ks) / 2.0) * ((Dx(Pol[y]) - Dy(Pol[x])) * (Dx(Pol[y]) - Dy(Pol[x])));
            Particles.ghost_get<33>();


            f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Particles.ghost_get<11, 12, 13, 14, 15, 16>();
            Df1[x] = Dx(f1);
            Df2[x] = Dx(f2);
            Df3[x] = Dx(f3);
            Df4[x] = Dx(f4);
            Df5[x] = Dx(f5);
            Df6[x] = Dx(f6);

            Df1[y] = Dy(f1);
            Df2[y] = Dy(f2);
            Df3[y] = Dy(f3);
            Df4[y] = Dy(f4);
            Df5[y] = Dy(f5);
            Df6[y] = Dy(f6);
            Particles.ghost_get<21, 22, 23, 24, 25, 26>();


            dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                    zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) -
                    g[x]
                    - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));


            dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                    zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) -
                    g[y]
                    - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
            Particles.ghost_get<9>();


            //Particles.write("PolarI");
            //Velocity Solution n iterations


            auto Stokes1 = eta * (Dxx(V[x]) + Dyy(V[x]))
                           + 0.5 * nu * (Df1[x] * Dx(V[x]) + f1 * Dxx(V[x]))
                           + 0.5 * nu * (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y]))
                           + 0.5 * nu * (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x]))
                           + 0.5 * nu * (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           + 0.5 * nu * (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
            auto Stokes2 = eta * (Dxx(V[y]) + Dyy(V[y]))
                           - 0.5 * nu * (Df1[y] * Dx(V[x]) + f1 * Dxy(V[x]))
                           - 0.5 * nu * (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           - 0.5 * nu * (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y]))
                           + 0.5 * nu * (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x]))
                           + 0.5 * nu * (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));
            tt.stop();
            std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
            tt.start();
            V_err = 1;
            n = 0;
            errctr = 0;
            while (V_err >= V_err_eps && n <= nmax) {

                petsc_solver<double> solverPetsc;
                solverPetsc.setSolver(KSPGMRES);
                //solverPetsc.setRestart(250);
                solverPetsc.setPreconditioner(PCJACOBI);
                petsc_solver<double> solverPetsc2;
                solverPetsc2.setSolver(KSPGMRES);
                solverPetsc2.setPreconditioner(PCJACOBI);

                RHS[x] = dV[x];
                RHS[y] = dV[y];
                Particles_subset.ghost_get<0>();
                Grad_bulk[x] = Bulk_Dx(P_bulk);
                Grad_bulk[y] = Bulk_Dy(P_bulk);
                Particles_subset.ghost_get<2>();
                for (int i = 0; i < bulk.size(); i++) {
                    Particles.template getProp<10>(bulk.template get<0>(i))[x] += Particles_subset.getProp<2>(i)[x];
                    Particles.template getProp<10>(bulk.template get<0>(i))[y] += Particles_subset.getProp<2>(i)[y];
                }
                Particles.ghost_get<10>();
                DCPSE_scheme<equations2d2, decltype(Particles)> Solver(Particles);
                Solver.impose(Stokes1, bulk, RHS[0], vx);
                Solver.impose(Stokes2, bulk, RHS[1], vy);
                Solver.impose(V[x], up_p, 0, vx);
                Solver.impose(V[y], up_p, 0, vy);
                Solver.impose(V[x], dw_p, 0, vx);
                Solver.impose(V[y], dw_p, 0, vy);
                Solver.impose(V[x], l_p, 0, vx);
                Solver.impose(V[y], l_p, 0, vy);
                Solver.impose(V[x], r_p, 0, vx);
                Solver.impose(V[y], r_p, 0, vy);
                //Solver.solve(V[x], V[y]);
                Solver.solve_with_solver(solverPetsc, V[x], V[y]);
                Particles.ghost_get<Velocity>();
                div = -(Dx(V[x]) + Dy(V[y]));
                Particles.ghost_get<19>();
                auto Helmholtz = Dxx(H) + Dyy(H);
                DCPSE_scheme<equations2d1, decltype(Particles)> SolverH(Particles);
                SolverH.impose(Helmholtz, bulk, prop_id<19>());
                SolverH.impose(H, up_p, 0);
                SolverH.impose(H, dw_p, 0);
                SolverH.impose(H, l_p, 0);
                SolverH.impose(H, r_p, 0);
                //SolverH.solve(H);
                SolverH.solve_with_solver(solverPetsc2, H);
                //std::cout << "Helmholtz Solved in " << tt.getwct() << " seconds." << std::endl;
                Particles.ghost_get<17>();
                Particles.ghost_get<Velocity>();
                P = P + div;
                Particles.ghost_get<Pressure>();
                for (int i = 0; i < bulk.size(); i++) {
                    Particles_subset.getProp<0>(i) = Particles.template getProp<4>(bulk.template get<0>(i));
                }
                V[x] = V[x] + Dx(H);
                V[y] = V[y] + Dy(H);
                for (int j = 0; j < up_p.size(); j++) {
                    auto p = up_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    Particles.getProp<4>(p) = 0;

                }
                for (int j = 0; j < dw_p.size(); j++) {
                    auto p = dw_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    Particles.getProp<4>(p) = 0;
                }
                for (int j = 0; j < l_p.size(); j++) {
                    auto p = l_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    Particles.getProp<4>(p) = 0;
                }
                for (int j = 0; j < r_p.size(); j++) {
                    auto p = r_p.get<0>(j);
                    Particles.getProp<1>(p)[0] = 0;
                    Particles.getProp<1>(p)[1] = 0;
                    Particles.getProp<4>(p) = 0;
                }
                Particles.ghost_get<Velocity>();
                Particles.ghost_get<Pressure>();
                sum = 0;
                sum1 = 0;
                for (int j = 0; j < bulk.size(); j++) {
                    auto p = bulk.get<0>(j);
                    sum += (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) *
                           (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) +
                           (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]) *
                           (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]);
                    sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                            Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1];
                }
                sum = sqrt(sum);
                sum1 = sqrt(sum1);

                v_cl.sum(sum);
                v_cl.sum(sum1);
                v_cl.execute();
                V_t = V;
                Particles.ghost_get<18>();
                V_err_old = V_err;
                V_err = sum / sum1;
                if (V_err > V_err_old) {
                    errctr++;
                    //alpha_P -= 0.1;
                } else {
                    errctr = 0;
                }
                if (errctr > 5)
                    break;
                n++;
                if (v_cl.rank() == 0) { std::cout << "Rel l2 cgs err in V = " << V_err << std::endl; }
            }
            tt.stop();
            std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds with " << n
                      << " iterations."
                      << std::endl;

            u[x][x] = Dx(V[x]);
            u[x][y] = 0.5 * (Dx(V[y]) + Dy(V[x]));
            u[y][x] = 0.5 * (Dy(V[x]) + Dx(V[y]));
            u[y][y] = Dy(V[y]);

            W[x][x] = 0;
            W[x][y] = 0.5 * (Dy(V[x]) - Dx(V[y]));
            W[y][x] = 0.5 * (Dx(V[y]) - Dy(V[x]));
            W[y][y] = 0;


            h[x] = -gama * (lambda * delmu - nu * (u[x][x] * Pol[x] * Pol[x] +  u[y][y] * Pol[y] * Pol[y] + 2* u[x][y] * Pol[x] * Pol[y]) / (Pol[x] * Pol[x] + Pol[y] * Pol[y]));


            Particles.ghost_get<MolField, Strain_rate, Vorticity>();
            Particles.deleteGhost();
            Particles.write_frame("Polar_3e-3", ctr);
            Particles.ghost_get<MolField, Strain_rate, Vorticity>();
            ctr++;
            H_p_b = Pol[x] * Pol[x] + Pol[y] * Pol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<32>(p) =(Particles.getProp<32>(p)== 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<32>(p) =(Particles.getProp<32>(p)== 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<32>(p) =(Particles.getProp<32>(p)== 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<32>(p) =(Particles.getProp<32>(p)== 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<32>(p) =(Particles.getProp<32>(p)== 0) ? 1 : Particles.getProp<32>(p);
            }
            H_p_b=sqrt(H_p_b);

            k1[x] = ((h[x] * Pol[x] - h[y] * Pol[y]) / gama + lambda * delmu * Pol[x] -
                     nu * (u[x][x] * Pol[x] + u[x][y] * Pol[y]) + W[x][x] * Pol[x] +
                     W[x][y] * Pol[y]);// - V[x] * Dx(Pol[x]) - V[y] * Dy(Pol[x]));
            k1[y] = ((h[x] * Pol[y] + h[y] * Pol[x]) / gama + lambda * delmu * Pol[y] -
                     nu * (u[y][x] * Pol[x] + u[y][y] * Pol[y]) + W[y][x] * Pol[x] +
                     W[y][y] * Pol[y]);// - V[x] * Dx(Pol[y]) - V[y] * Dy(Pol[y]));
            Particles.ghost_get<28>();

            H_t = H_p_b;//+0.5*dt*(k1[x]*k1[x]+k1[y]*k1[y]);
            dPol = Pol + (0.5 * dt) * k1;
            dPol = dPol / H_t;
            Particles.ghost_get<8>();
            r=dPol[x]*dPol[x]+dPol[y]*dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] +  u[y][y] * dPol[y] * dPol[y] + 2* u[x][y] * dPol[x] * dPol[y])/(r)));
            Particles.ghost_get<7>();


            k2[x] = ((h[x] * (dPol[x]) - h[y] * (dPol[y])) / gama +
                     lambda * delmu * (dPol[x]) -
                     nu * (u[x][x] * (dPol[x]) + u[x][y] * (dPol[y])) +
                     W[x][x] * (dPol[x]) + W[x][y] * (dPol[y])); //-V[x] * Dx((dPol[x])) - V[y] * Dy((dPol[x])));
            k2[y] = ((h[x] * (dPol[y]) + h[y] * (dPol[x])) / gama +
                     lambda * delmu * (dPol[y]) -
                     nu * (u[y][x] * (dPol[x]) + u[y][y] * (dPol[y])) +
                     W[y][x] * (dPol[x]) + W[y][y] * (dPol[y])); //-V[x] * Dx((dPol[y])) - V[y] * Dy((dPol[y])));
            Particles.ghost_get<29>();
            H_t = H_p_b;//+0.5*dt*(k2[x]*k2[x]+k2[y]*k2[y]);
            dPol = Pol + (0.5 * dt) * k2;
            dPol = dPol / H_t;
            Particles.ghost_get<8>();
            r=dPol[x]*dPol[x]+dPol[y]*dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] +  u[y][y] * dPol[y] * dPol[y] + 2* u[x][y] * dPol[x] * dPol[y])/(r)));

            k3[x] = ((h[x] * (dPol[x]) - h[y] * (dPol[y])) / gama +
                     lambda * delmu * (dPol[x]) -
                     nu * (u[x][x] * (dPol[x]) + u[x][y] * (dPol[y])) +
                     W[x][x] * (dPol[x]) + W[x][y] * (dPol[y]));
            // -V[x] * Dx((dPol[x])) - V[y] * Dy((dPol[x])));
            k3[y] = ((h[x] * (dPol[y]) + h[y] * (dPol[x])) / gama +
                     lambda * delmu * (dPol[y]) -
                     nu * (u[y][x] * (dPol[x]) + u[y][y] * (dPol[y])) +
                     W[y][x] * (dPol[x]) + W[y][y] * (dPol[y]));
            // -V[x] * Dx((dPol[y])) - V[y] * Dy((dPol[y])));
            Particles.ghost_get<30>();
            H_t =  H_p_b;//+dt*(k3[x]*k3[x]+k3[y]*k3[y]);
            dPol = Pol + (dt * k3);
            dPol = dPol / H_t;
            Particles.ghost_get<8>();
            r=dPol[x]*dPol[x]+dPol[y]*dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) =(Particles.getProp<34>(p)== 0) ? 1 : Particles.getProp<34>(p);
            }

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] +  u[y][y] * dPol[y] * dPol[y] + 2* u[x][y] * dPol[x] * dPol[y])/(r)));
            Particles.ghost_get<7>();


            k4[x] = ((h[x] * (dPol[x]) - h[y] * (dPol[y])) / gama +
                     lambda * delmu * (dPol[x]) -
                     nu * (u[x][x] * (dPol[x]) + u[x][y] * (dPol[y])) +
                     W[x][x] * (dPol[x]) +
                     W[x][y] * (dPol[y]));//   -V[x]*Dx( (dt * k3[x]+Pol[x])) -V[y]*Dy( (dt * k3[x]+Pol[x])));
            k4[y] = ((h[x] * (dPol[y]) + h[y] * (dPol[x])) / gama +
                     lambda * delmu * (dPol[y]) -
                     nu * (u[y][x] * (dPol[x]) + u[y][y] * (dPol[y])) +
                     W[y][x] * (dPol[x]) +
                     W[y][y] * (dPol[y]));//  -V[x]*Dx( (dt * k3[y]+Pol*[y])) -V[y]*Dy( (dt * k3[y]+Pol[y])));
            Particles.ghost_get<31>();

            Pol = Pol + (dt / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
            Pol = Pol / H_p_b;

            H_p_b = sqrt(Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Pol = Pol / H_p_b;

            Pos = Pos + dt * V;
            Particles.map();
            Particles.ghost_get<0>();
            indexUpdate(Particles, Particles_subset, up_p, dw_p, l_p, r_p, bulk, up, down, left, right);
            Particles_subset.map();
            Particles_subset.ghost_get<0>();

            //Particles_subset.write("debug");
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }

            Particles.ghost_get<0, Vorticity, MolField>();
            Particles_subset.ghost_get<0,1,2>();

            tt.start();
            Dx.update(Particles);
            Dy.update(Particles);
            Dxy.update(Particles);
            auto Dyx = Dxy;
            Dxx.update(Particles);
            Dyy.update(Particles);

            Bulk_Dx.update(Particles_subset);
            Bulk_Dy.update(Particles_subset);

            tt.stop();
            std::cout << "Updation of operators took " << tt.getwct() << " seconds." << std::endl;

            std::cout << "Time step " << ctr << " : " << tim << " over." << std::endl;
            tim += dt;
            std::cout << "----------------------------------------------------------" << std::endl;
        }
        Particles.deleteGhost();
        tt2.stop();
        std::cout << "The simulation took " << tt2.getwct() << "Seconds.";
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(Active2DEigenP) {
        const size_t sz[2] = {21, 21};
        Box<2, double> box({0, 0}, {10, 10});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        int ord = 2;
        double rCut = 3.1 * spacing;
        double sampling_factor = 1.9;


        int ord2 = 2;
        double rCut2 = 3.1 * spacing;
        double sampling_factor2 = 1.9;

        Ghost<2, double> ghost(rCut);

/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div   H_t   */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2]>> Particles(
                0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y * 0.99999;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkP;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        constexpr int x = 0;
        constexpr int y = 1;

        constexpr int Polarization = 0;
        constexpr int Velocity = 1;
        constexpr int Vorticity = 2;
        constexpr int ExtForce = 3;
        constexpr int Pressure = 4;
        constexpr int Strain_rate = 5;
        constexpr int Stress = 6;
        constexpr int MolField = 7;

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        V.setVarId(0);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        P.setVarId(0);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dP = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);
        auto Df1 = getV<21>(Particles);
        auto Df2 = getV<22>(Particles);
        auto Df3 = getV<23>(Particles);
        auto Df4 = getV<24>(Particles);
        auto Df5 = getV<25>(Particles);
        auto Df6 = getV<26>(Particles);


        double eta = 1.0;
        double nu = -0.5;
        double gama = 0.1;
        double zeta = 0.07;
        double Ks = 1.0;
        double Kb = 1.0;
        double lambda = 0.1;
        double delmu = -1.0;
        g = 0;
        V = 0;
        P = 0;

        // Here fill up the boxes for particle boundary detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        Box<2, double> mid({box.getHigh(0) / 2.0 - spacing, box.getHigh(1) / 2.0 - spacing / 2.0},
                           {box.getHigh(0) / 2.0, box.getHigh(1) / 2.0 + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            } else {
                if (mid.isInside(xp) == true) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                    Particles.getProp<4>(p) = 0;
                } else {
                    bulkP.add();
                    bulkP.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }
            ++it2;
        }


        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Dx2(Particles, ord2, rCut2,
                                                                                             sampling_factor2,
                                                                                             support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Dy2(Particles, ord2, rCut2,
                                                                                             sampling_factor2,
                                                                                             support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Gradient Grad(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Laplacian Lap(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Advection Adv(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Divergence Div(Particles, ord, rCut, sampling_factor, support_options::RADIUS);

        Particles.ghost_get<Polarization>();
        sigma[x][x] =
                -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
        sigma[x][y] =
                -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
        sigma[y][x] =
                -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
        sigma[y][y] =
                -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
        Particles.ghost_get<Stress>();


        h[y] = Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
               Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y]));
        Particles.ghost_get<MolField>();


        f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
        f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
             (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
        f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
        f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
        f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
        f6 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
        Particles.ghost_get<11, 12, 13, 14, 15, 16>();
        Df1[x] = Dx(f1);
        Df2[x] = Dx(f2);
        Df3[x] = Dx(f3);
        Df4[x] = Dx(f4);
        Df5[x] = Dx(f5);
        Df6[x] = Dx(f6);
        Df1[y] = Dy(f1);
        Df2[y] = Dy(f2);
        Df3[y] = Dy(f3);
        Df4[y] = Dy(f4);
        Df5[y] = Dy(f5);
        Df6[y] = Dy(f6);
        Particles.ghost_get<21, 22, 23, 24, 25, 26>();


        dV[x] = -0.5 * Dy(h[y])
                + zeta * Dx(delmu * Pol[x] * Pol[x])
                + zeta * Dy(delmu * Pol[x] * Pol[y])
                - zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y]))
                - 0.5 * nu * Dx(-2 * h[y] * Pol[x] * Pol[y])
                - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) - g[x]
                - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) -
                0.5 * Dy(-2 * gama * lambda * delmu * (Pol[x] * Pol[y]));


        dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                0.5 * nu * Dy(-2 * h[y] * Pol[x] * Pol[y])
                - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) - g[y]
                - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) -
                0.5 * Dx(-2 * gama * lambda * delmu * (Pol[x] * Pol[y]));
        Particles.ghost_get<9>();


        Particles.write_frame("Polar", 0);
        //Velocity Solution n iterations
        eq_id vx, vy;
        timer tt;
        vx.setId(0);
        vy.setId(1);
        double sum = 0, sum1 = 0;
        int n = 100;
        /* auto Stokes1 = nu * Lap(V[x]) + 0.5 * nu * (f1 * Dxx(V[x]) + Dx(f1) * Dx(V[x])) +
                        (Dx(f2) * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x]))) +
                        (Dx(f3) * Dy(V[y]) + f3 * Dyx(V[y])) + (Dy(f4) * Dx(V[x]) + f4 * Dxy(V[x])) +
                        (Dy(f5) * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x]))) +
                        (Dy(f6) * Dy(V[y]) + f6 * Dyy(V[y]));
         auto Stokes2 = nu * Lap(V[y]) + 0.5 * nu * (f1 * Dxy(V[x]) + Dy(f1) * Dx(V[x])) +
                        (Dy(f2) * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x]))) +
                        (Dy(f3) * Dy(V[y]) + f3 * Dyy(V[y])) + (Dx(f4) * Dx(V[x]) + f4 * Dxx(V[x])) +
                        (Dx(f5) * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x]))) +
                        (Dx(f6) * Dy(V[y]) + f6 * Dyx(V[y]));
 */
        auto Stokes1 = eta * Lap(V[x]) + 0.5 * nu * (f1 * Dxx(V[x]) + Df1[x] * Dx(V[x])) +
                       (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x]))) +
                       (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y])) + (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x])) +
                       (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x]))) +
                       (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
        auto Stokes2 = eta * Lap(V[y]) + 0.5 * nu * (f1 * Dxy(V[x]) + Df1[y] * Dx(V[x])) +
                       (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x]))) +
                       (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y])) + (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x])) +
                       (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x]))) +
                       (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));

        auto Helmholtz = Lap(H);
        auto D_y = Dy2(H);
        auto D_x = Dx2(H);

        for (int i = 1; i <= n; i++) {
            RHS[x] = Dx(P) + dV[x];
            RHS[y] = Dy(P) + dV[y];
            Particles.ghost_get<10>();
            DCPSE_scheme<equations2d2pE, decltype(Particles)> Solver(Particles);
            Solver.impose(Stokes1, bulk, RHS[0], vx);
            Solver.impose(Stokes2, bulk, RHS[1], vy);
            Solver.impose(V[x], up_p, 0, vx);
            Solver.impose(V[y], up_p, 0, vy);
            Solver.impose(V[x], dw_p, 0, vx);
            Solver.impose(V[y], dw_p, 0, vy);
            tt.start();
            Solver.solve(V[x], V[y]);
            tt.stop();
            std::cout << "Stokes Solved in " << tt.getwct() << " seconds." << std::endl;
            Particles.ghost_get<Velocity>();
            div = -Div(V);
            Particles.ghost_get<19>();
            DCPSE_scheme<equations2d1pE, decltype(Particles)> SolverH(Particles, options_solver::LAGRANGE_MULTIPLIER);
            SolverH.impose(Helmholtz, bulk, prop_id<19>());
            SolverH.impose(H, up_p, 0);
            SolverH.impose(H, dw_p, 0);
            tt.start();
            SolverH.solve(H);
            tt.stop();
            std::cout << "Helmholtz Solved in " << tt.getwct() << " seconds." << std::endl;
            Particles.ghost_get<17>();
            V = V + Grad(H);
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            }
            P = P + Lap(H);
            Particles.ghost_get<Velocity>();
            Particles.ghost_get<Pressure>();
            sum = 0;
            sum1 = 0;
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) *
                       (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) +
                       (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]) *
                       (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]);
                sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                        Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1];
            }
            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_t = V;
            std::cout << "Rel l2 cgs err in V at " << i << "= " << sum / sum1 << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            if (i % 10 == 0)
                Particles.write_frame("Polar", i);
        }
        Particles.deleteGhost();
        Particles.write_frame("Polar", n + 1);

    }


    BOOST_AUTO_TEST_CASE(Active2DEigen_decouple) {
        timer tt2;
        tt2.start();
        const size_t sz[2] = {41, 41};
        Box<2, double> box({0, 0}, {10, 10});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        double rCut2 = 3.1 * spacing;
        int ord = 2;
        double sampling_factor = 1.9;
        int ord2 = 2;
        double sampling_factor2 = 1.9;
        Ghost<2, double> ghost(rCut);
/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div   H_t                                                                                      delmu */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2], double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>>> Particles(
                0, box, bc, ghost);

        vector_dist<2, double, aggregate<double, double, VectorS<2, double>>> Particles_subset(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkP;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        constexpr int x = 0;
        constexpr int y = 1;

        constexpr int Polarization = 0;
        constexpr int Velocity = 1;
        constexpr int Vorticity = 2;
        constexpr int ExtForce = 3;
        constexpr int Pressure = 4;
        constexpr int Strain_rate = 5;
        constexpr int Stress = 6;
        constexpr int MolField = 7;
        auto Pos = getV<PROP_POS>(Particles);

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dP = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);
        auto Df1 = getV<21>(Particles);
        auto Df2 = getV<22>(Particles);
        auto Df3 = getV<23>(Particles);
        auto Df4 = getV<24>(Particles);
        auto Df5 = getV<25>(Particles);
        auto Df6 = getV<26>(Particles);
        auto delmu = getV<27>(Particles);
        auto k1 = getV<28>(Particles);
        auto k2 = getV<29>(Particles);
        auto k3 = getV<30>(Particles);
        auto k4 = getV<31>(Particles);

        auto P_bulk = getV<0>(Particles_subset); //Pressure only on inside
        auto H_bulk = getV<1>(Particles_subset); //Pressure only on inside
        auto Grad_bulk = getV<2>(Particles_subset);


        double eta = 1.0;
        double nu = -0.5;
        double gama = 0.1;
        double zeta = 0.07;
        double Ks = 1.0;
        double Kb = 1.0;
        double lambda = 0.1;
        //double delmu = -1.0;
        g = 0;
        delmu = -1.0;
        P = 0;
        V = 0;
        P_bulk = 0;
        // Here fill up the boxes for particle boundary detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        Box<2, double> mid({box.getHigh(0) / 2.0 - 0.75 * spacing, box.getHigh(1) / 2.0 - 0.75 * spacing},
                           {box.getHigh(0) / 2.0 + 0.75 * spacing, box.getHigh(1) / 2.0 + 0.75 * spacing});

/*
        Box<2, double> mid({box.getLow(0)/ + 3.1 * spacing, box.getLow(1) + 3.1 * spacing},
                           {box.getHigh(0) - 3.1 * spacing, box.getHigh(1) - 3.1 * spacing});
*/


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            } else {
                if (mid.isInside(xp) == true) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                    Particles.getProp<4>(p) = 0;
                } else {
                    bulkP.add();
                    bulkP.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }


            ++it2;
        }
        for (int i = 0; i < bulk.size(); i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
        }

        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dx(Particles_subset, ord2,
                                                                                                 rCut2,
                                                                                                 sampling_factor2,
                                                                                                 support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dy(Particles_subset, ord2,
                                                                                                 rCut2,
                                                                                                 sampling_factor2,
                                                                                                 support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dxx(Particles_subset,
                                                                                                    ord2, rCut2,
                                                                                                    sampling_factor2,
                                                                                                    support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dyy(Particles_subset,
                                                                                                    ord2, rCut2,
                                                                                                    sampling_factor2,
                                                                                                    support_options::RADIUS);

        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        //solverPetsc.setPreconditioner(PCJACOBI);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setSolver(KSPGMRES);
        //solverPetsc2.setRestart(250);
        /*solverPetsc2.setPreconditioner(PCJACOBI);
        solverPetsc2.setRelTol(1e-6);
        solverPetsc2.setAbsTol(1e-5);
        solverPetsc2.setDivTol(1e4);*/
        eq_id vx, vy, ic;
        timer tt;
        timer tt3;
        vx.setId(0);
        vy.setId(1);
        double V_err_eps = 1e-4;
        double V_err = 1;
        int n = 0;
        int ctr = 0;
        double dt = 3e-3;
        double tim = 0;
        double tf = 3e-3;
        while (tim < tf) {
            tt.start();
            Particles.ghost_get<Polarization>();
            sigma[x][x] =
                    -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
            sigma[x][y] =
                    -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
            sigma[y][x] =
                    -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
            sigma[y][y] =
                    -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
            Particles.ghost_get<Stress>();


            h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));
            Particles.ghost_get<MolField>();

            f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Particles.ghost_get<11, 12, 13, 14, 15, 16>();
            Df1[x] = Dx(f1);
            Df2[x] = Dx(f2);
            Df3[x] = Dx(f3);
            Df4[x] = Dx(f4);
            Df5[x] = Dx(f5);
            Df6[x] = Dx(f6);

            Df1[y] = Dy(f1);
            Df2[y] = Dy(f2);
            Df3[y] = Dy(f3);
            Df4[y] = Dy(f4);
            Df5[y] = Dy(f5);
            Df6[y] = Dy(f6);
            Particles.ghost_get<21, 22, 23, 24, 25, 26>();


            dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                    zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) -
                    g[x]
                    - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));


            dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                    zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) -
                    g[y]
                    - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
            Particles.ghost_get<9>();
            //Velocity Solution
            auto Stokes1 = eta * (Dxx(V[x]) + Dyy(V[x]))
                           + 0.5 * nu * (Df1[x] * Dx(V[x]) + f1 * Dxx(V[x]))
                           + 0.5 * nu * (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y]))
                           + 0.5 * nu * (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x]))
                           + 0.5 * nu * (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           + 0.5 * nu * (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
            auto Stokes2 = eta * (Dxx(V[y]) + Dyy(V[y]))
                           - 0.5 * nu * (Df1[y] * Dx(V[x]) + f1 * Dxy(V[x]))
                           - 0.5 * nu * (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           - 0.5 * nu * (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y]))
                           + 0.5 * nu * (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x]))
                           + 0.5 * nu * (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));
            tt.stop();
            std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
            tt.start();
            //auto Pressure_Poisson = Bulk_Dxx(P_bulk) + Bulk_Dyy(P_bulk);
            auto Pressure_Poisson = Dxx(P) + Dyy(P);
            div = -(Dx(dV[x]) + Dy(dV[y]));
            /*for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_subset.getProp<1>(i) = Particles.template getProp<19>(bulk.template get<0>(i));
            }*/
            Particles.ghost_get<19>();
            DCPSE_scheme<equations2d1E, decltype(Particles)> SolverH(Particles);
            SolverH.impose(Pressure_Poisson, bulkP, prop_id<19>());
            SolverH.impose(Dy(P), up_p, 0);
            SolverH.impose(Dy(P), dw_p, 0);
            SolverH.impose(Dx(P), l_p, 0);
            SolverH.impose(Dx(P), r_p, 0);
            SolverH.impose(P, ref_p, 0);
            SolverH.solve(P);
            std::cout << "Pressure Poisson Solved in " << tt.getwct() << " seconds." << std::endl;
            Particles.ghost_get<Pressure>();
            for (int i = 0; i < bulk.size(); i++) {
                Particles_subset.getProp<0>(i) = Particles.template getProp<4>(bulk.template get<0>(i));
            }
            RHS[x] = dV[x];
            RHS[y] = dV[y];
            Grad_bulk[x] = Bulk_Dx(P_bulk);
            Grad_bulk[y] = Bulk_Dy(P_bulk);
            for (int i = 0; i < bulk.size(); i++) {
                Particles.template getProp<10>(bulk.template get<0>(i))[x] += Particles_subset.getProp<2>(i)[x];
                Particles.template getProp<10>(bulk.template get<0>(i))[y] += Particles_subset.getProp<2>(i)[y];
            }
            Particles.ghost_get<10>();
            DCPSE_scheme<equations2d2E, decltype(Particles)> Solver(Particles);
            Solver.impose(Stokes1, bulk, RHS[0], vx);
            Solver.impose(Stokes2, bulk, RHS[1], vy);
            Solver.impose(V[x], up_p, 0, vx);
            Solver.impose(V[y], up_p, 0, vy);
            Solver.impose(V[x], dw_p, 0, vx);
            Solver.impose(V[y], dw_p, 0, vy);
            Solver.impose(V[x], l_p, 0, vx);
            Solver.impose(V[y], l_p, 0, vy);
            Solver.impose(V[x], r_p, 0, vx);
            Solver.impose(V[y], r_p, 0, vy);
            tt.start();
            Solver.solve(V[x], V[y]);
            tt.stop();
            for (int i = 0; i < bulk.size(); i++) {
                Particles.template getProp<4>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            }

            std::cout << "Stokes Solved in " << tt.getwct() << " seconds." << std::endl;
            std::cout << "----------------------------------------------------------" << std::endl;
            Particles.write("Polar_decouple");
            tim += dt;
        }

    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(Active2DEigen_saddle) {
        timer tt2;
        tt2.start();
        const size_t sz[2] = {41, 41};
        Box<2, double> box({0, 0}, {10, 10});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.6 * spacing;
        double rCut2 = 3.6 * spacing;
        int ord = 2;
        double sampling_factor = 2.5;
        int ord2 = 2;
        double sampling_factor2 = 2.5;
        Ghost<2, double> ghost(rCut);
/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div   H_t                                                                                      delmu */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2], double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>>> Particles(
                0, box, bc, ghost);


        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkP;
        openfpm::vector<aggregate<int>> bulkF;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        constexpr int x = 0;
        constexpr int y = 1;

        constexpr int Polarization = 0;
        constexpr int Velocity = 1;
        constexpr int Vorticity = 2;
        constexpr int ExtForce = 3;
        constexpr int Pressure = 4;
        constexpr int Strain_rate = 5;
        constexpr int Stress = 6;
        constexpr int MolField = 7;
        auto Pos = getV<PROP_POS>(Particles);

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dP = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);
        auto Df1 = getV<21>(Particles);
        auto Df2 = getV<22>(Particles);
        auto Df3 = getV<23>(Particles);
        auto Df4 = getV<24>(Particles);
        auto Df5 = getV<25>(Particles);
        auto Df6 = getV<26>(Particles);
        auto delmu = getV<27>(Particles);
        auto k1 = getV<28>(Particles);
        auto k2 = getV<29>(Particles);
        auto k3 = getV<30>(Particles);
        auto k4 = getV<31>(Particles);


        double eta = 1.0;
        double nu = -0.5;
        double gama = 0.1;
        double zeta = 0.07;
        double Ks = 1.0;
        double Kb = 1.0;
        double lambda = 0.1;
        //double delmu = -1.0;
        g = 0;
        delmu = -1.0;
        P = 0;
        V = 0;
        // Here fill up the boxes for particle boundary detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        Box<2, double> mid({box.getHigh(0) / 2.0 - spacing, box.getHigh(1) / 2.0 - spacing / 2.0},
                           {box.getHigh(0) / 2.0, box.getHigh(1) / 2.0 + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();
            } else {
                if (mid.isInside(xp) == true) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                    Particles.getProp<4>(p) = 0;
                } else {
                    bulkP.add();
                    bulkP.last().get<0>() = p.getKey();
                    bulkF.add();
                    bulkF.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }


            ++it2;
        }
        Derivative_x Dx(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_y Dy(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        V.setVarId(0);
        P.setVarId(2);
        eq_id vx, vy, ic;
        timer tt;
        timer tt3;
        vx.setId(0);
        vy.setId(1);
        ic.setId(2);
        double V_err_eps = 1e-4;
        double V_err = 1;
        int n = 0;
        int ctr = 0;
        double dt = 3e-3;
        double tim = 0;
        double tf = 3e-3;
        while (tim < tf) {
            tt.start();
            Particles.ghost_get<Polarization>();
            sigma[x][x] =
                    -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
            sigma[x][y] =
                    -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
            sigma[y][x] =
                    -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
            sigma[y][y] =
                    -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
            Particles.ghost_get<Stress>();


            h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));
            Particles.ghost_get<MolField>();

            f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) /
                 (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Particles.ghost_get<11, 12, 13, 14, 15, 16>();
            Df1[x] = Dx(f1);
            Df2[x] = Dx(f2);
            Df3[x] = Dx(f3);
            Df4[x] = Dx(f4);
            Df5[x] = Dx(f5);
            Df6[x] = Dx(f6);

            Df1[y] = Dy(f1);
            Df2[y] = Dy(f2);
            Df3[y] = Dy(f3);
            Df4[y] = Dy(f4);
            Df5[y] = Dy(f5);
            Df6[y] = Dy(f6);
            Particles.ghost_get<21, 22, 23, 24, 25, 26>();


            dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                    zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) -
                    g[x]
                    - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));


            dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                    zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) -
                    g[y]
                    - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
            Particles.ghost_get<9>();
            auto Stokes1 = eta * (Dxx(V[x]) + Dyy(V[x]))
                           + 0.5 * nu * (Df1[x] * Dx(V[x]) + f1 * Dxx(V[x]))
                           + 0.5 * nu * (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y]))
                           + 0.5 * nu * (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x]))
                           + 0.5 * nu * (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           + 0.5 * nu * (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
            auto Stokes2 = eta * (Dxx(V[y]) + Dyy(V[y]))
                           - 0.5 * nu * (Df1[y] * Dx(V[x]) + f1 * Dxy(V[x]))
                           - 0.5 * nu * (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           - 0.5 * nu * (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y]))
                           + 0.5 * nu * (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x]))
                           + 0.5 * nu * (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));

            auto continuity = Dx(V[x]) + Dy(V[y]) + 1e-2 * (Dxx(P) + Dyy(P));
            tt.stop();
            std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
            tt.start();
            Particles.ghost_get<Pressure>();
            RHS[x] = dV[x];
            RHS[y] = dV[y];
            Particles.ghost_get<10>();
            DCPSE_scheme<equations2d3E, decltype(Particles)> Solver(Particles);
            Solver.impose(Stokes1, bulk, RHS[0], vx);
            Solver.impose(Stokes2, bulk, RHS[1], vy);
            Solver.impose(continuity, bulkF, 0, ic);
            Solver.impose(V[x], up_p, 0, vx);
            Solver.impose(V[y], up_p, 0, vy);
            Solver.impose(V[x], dw_p, 0, vx);
            Solver.impose(V[y], dw_p, 0, vy);
            Solver.impose(V[x], l_p, 0, vx);
            Solver.impose(V[y], l_p, 0, vy);
            Solver.impose(V[x], r_p, 0, vx);
            Solver.impose(V[y], r_p, 0, vy);
            Solver.impose(P, ref_p, 0, ic);
            tt.start();
            BOOST_TEST_MESSAGE("Solver Imposed...");
            Solver.solve(V[x], V[y], P);
            BOOST_TEST_MESSAGE("Solver Solved...");
            //Solver.solve(V[x], V[y],P);
            tt.stop();
            std::cout << "Stokes Solved in " << tt.getwct() << " seconds." << std::endl;
            Particles.ghost_get<Velocity>();
            std::cout << "----------------------------------------------------------" << std::endl;
            Particles.write_frame("Polar_saddle", 0);
            Particles.deleteGhost();
            Particles.write_frame("Polar_saddle", 1);
            tt2.stop();
            std::cout << "Simulation took " << tt2.getwct() << " seconds." << std::endl;
            tim += dt;
        }
    }


BOOST_AUTO_TEST_SUITE_END()
