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
        particle_type &Particles,
        particle_type2 &Particles_subset,
        openfpm::vector<aggregate<int>> &up_p, openfpm::vector<aggregate<int>> &dw_p,
        openfpm::vector<aggregate<int>> &l_p, openfpm::vector<aggregate<int>> &r_p,
        openfpm::vector<aggregate<int>> &up_p1, openfpm::vector<aggregate<int>> &dw_p1,
        openfpm::vector<aggregate<int>> &l_p1, openfpm::vector<aggregate<int>> &r_p1,
        openfpm::vector<aggregate<int>> &corner_ul, openfpm::vector<aggregate<int>> &corner_ur,
        openfpm::vector<aggregate<int>> &corner_dl, openfpm::vector<aggregate<int>> &corner_dr,
        openfpm::vector<aggregate<int>> &bulk, Box<2, double> &up, Box<2, double> &down, Box<2, double> &left,
        Box<2, double> &right) {
    up_p.clear();
    dw_p.clear();
    l_p.clear();
    r_p.clear();
    up_p1.clear();
    dw_p1.clear();
    l_p1.clear();
    r_p1.clear();
    corner_ul.clear();
    corner_ur.clear();
    corner_dl.clear();
    corner_dr.clear();
    bulk.clear();
    Particles_subset.clear();

    auto it2 = Particles.getDomainIterator();

    while (it2.isNext()) {
        auto p = it2.get();
        Point<2, double> xp = Particles.getPos(p);
        if (up.isInside(xp) == true) {
            if (left.isInside(xp) == true) {
                corner_ul.add();
                corner_ul.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                corner_ur.add();
                corner_ur.last().get<0>() = p.getKey();
            } else {
                up_p1.add();
                up_p1.last().get<0>() = p.getKey();
            }
            up_p.add();
            up_p.last().get<0>() = p.getKey();

        } else if (down.isInside(xp) == true) {
            if (left.isInside(xp) == true) {
                corner_dl.add();
                corner_dl.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                corner_dr.add();
                corner_dr.last().get<0>() = p.getKey();
            } else {
                dw_p1.add();
                dw_p1.last().get<0>() = p.getKey();
            }
            dw_p.add();
            dw_p.last().get<0>() = p.getKey();
        } else if (left.isInside(xp) == true) {
            if (up.isInside(xp) == false && down.isInside(xp) == false) {
                l_p1.add();
                l_p1.last().get<0>() = p.getKey();
            }
            l_p.add();
            l_p.last().get<0>() = p.getKey();
        } else if (right.isInside(xp) == true) {
            if (up.isInside(xp) == false && down.isInside(xp) == false) {
                r_p1.add();
                r_p1.last().get<0>() = p.getKey();
            }
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

}

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)

    BOOST_AUTO_TEST_CASE(Active2DEigen) {
        timer tt2;
        tt2.start();
        double boxsize = 10;
        const size_t sz[2] = {41, 41};
        Box<2, double> box({0, 0}, {boxsize, boxsize});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.9 * spacing;
        double rCut2 = 3.9 * spacing;
        int ord = 2;
        int ord2 = 2;
        double sampling_factor = 4.0;
        double sampling_factor2 = 2.4;
        double alpha_V = 1.0;
        double alpha_P = 1.0;
        Ghost<2, double> ghost(rCut);

        auto &v_cl = create_vcluster();

/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div    H_t                                                                                      delmu */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2], double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double>> Particles(
                0, box, bc, ghost);
        vector_dist<2, double, aggregate<double, double, VectorS<2, double>>> Particles_subset(
                Particles.getDecomposition(), 0);
        double x0, y0, x1, y1;
        x0 = box.getLow(0);
        y0 = box.getLow(1);
        x1 = box.getHigh(0);
        y1 = box.getHigh(1);

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

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Particles.map();
        Particles.ghost_get<0>();

        //Particles.write("Par");

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;
        openfpm::vector<aggregate<int>> corner_ul;
        openfpm::vector<aggregate<int>> corner_ur;
        openfpm::vector<aggregate<int>> corner_dl;
        openfpm::vector<aggregate<int>> corner_dr;


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
        Particles.ghost_get<ExtForce, 27>(SKIP_LABELLING);


        Box<2, double> up({x0 - spacing / 2.0, y1 - spacing / 2.0},
                          {x1 + spacing / 2.0, y1 + spacing / 2.0});

        Box<2, double> down({x0 - spacing / 2.0, y0 - spacing / 2.0},
                            {x1 + spacing / 2.0, y0 + spacing / 2.0});

        Box<2, double> left({x0 - spacing / 2.0, y0 - spacing / 2.0},
                            {x0 + spacing / 2.0, y1 + spacing / 2.0});

        Box<2, double> right({x1 - spacing / 2.0, y0 - spacing / 2.0},
                             {x1 + spacing / 2.0, y1 + spacing / 2.0});
        /*Box<2, double> mid({box.getHigh(0) / 2.0 - spacing, box.getHigh(1) / 2.0 - spacing / 2.0},
                           {box.getHigh(0) / 2.0, box.getHigh(1) / 2.0 + spacing / 2.0});
*/


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);


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
                if (left.isInside(xp) == true) {
                    corner_ul.add();
                    corner_ul.last().get<0>() = p.getKey();
                } else if (right.isInside(xp) == true) {
                    corner_ur.add();
                    corner_ur.last().get<0>() = p.getKey();
                } else {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                }
                up_p.add();
                up_p.last().get<0>() = p.getKey();

            } else if (down.isInside(xp) == true) {
                if (left.isInside(xp) == true) {
                    corner_dl.add();
                    corner_dl.last().get<0>() = p.getKey();
                } else if (right.isInside(xp) == true) {
                    corner_dr.add();
                    corner_dr.last().get<0>() = p.getKey();
                } else {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                }
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                if (up.isInside(xp) == false && down.isInside(xp) == false) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                }
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                if (up.isInside(xp) == false && down.isInside(xp) == false) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                }
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

        //Particles_subset.write("Pars");
        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dx(Particles_subset, ord,
                                                                                                 rCut,
                                                                                                 sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dy(Particles_subset, ord,
                                                                                                 rCut,
                                                                                                 sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);//, Dxx2(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);//, Dyy2(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);

/*        Derivative_x Dx(Particles, ord, rCut, sampling_factor), Bulk_Dx(Particles_subset, ord, rCut, sampling_factor);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor), Bulk_Dy(Particles_subset, ord, rCut, sampling_factor);
        Derivative_xy Dxy(Particles, ord2, rCut, sampling_factor);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord2, rCut2, sampling_factor2),Bulk_Dxx(Particles_subset, ord2, rCut2, sampling_factor2);
        Derivative_yy Dyy(Particles, ord2, rCut2, sampling_factor2),Bulk_Dyy(Particles_subset, ord2, rCut2, sampling_factor2);*/




        eq_id vx, vy;
        timer tt;
        timer tt3;
        vx.setId(0);
        vy.setId(1);
        double V_err_eps = 5 * 1e-4;
        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 300;
        int ctr = 0, errctr, Vreset = 0;
        double dt = 3 * 1e-3;

        double tim = 0;
        double tf = 1.25;
        div = 0;
        double sum, sum1, sum_k;
        while (tim <= tf) {
            tt.start();
            petsc_solver<double> solverPetsc;
            solverPetsc.setSolver(KSPGMRES);
            //solverPetsc.setRestart(250);
            solverPetsc.setPreconditioner(PCJACOBI);
            petsc_solver<double> solverPetsc2;
            solverPetsc2.setSolver(KSPGMRES);
            solverPetsc2.setPreconditioner(PCJACOBI);

            Particles.ghost_get<Polarization>(SKIP_LABELLING);
            sigma[x][x] =
                    -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
            sigma[x][y] =
                    -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
            sigma[y][x] =
                    -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
            sigma[y][y] =
                    -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
            Particles.ghost_get<Stress>(SKIP_LABELLING);


            r = Pol[x] * Pol[x] + Pol[y] * Pol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }

            h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));

            Particles.ghost_get<MolField>(SKIP_LABELLING);

            FranckEnergyDensity = (Ks / 2.0) *
                                  ((Dx(Pol[x]) * Dx(Pol[x])) + (Dy(Pol[x]) * Dy(Pol[x])) +
                                   (Dx(Pol[y]) * Dx(Pol[y])) +
                                   (Dy(Pol[y]) * Dy(Pol[y]))) +
                                  ((Kb - Ks) / 2.0) * ((Dx(Pol[y]) - Dy(Pol[x])) * (Dx(Pol[y]) - Dy(Pol[x])));
            Particles.ghost_get<33>(SKIP_LABELLING);


            f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
            f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
            f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
            f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (r);
            f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (r);
            f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (r);
            Particles.ghost_get<11, 12, 13, 14, 15, 16>(SKIP_LABELLING);
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
            Particles.ghost_get<21, 22, 23, 24, 25, 26>(SKIP_LABELLING);


            dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                    zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) -
                    Dy(sigma[x][y]) -
                    g[x]
                    - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));


            dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                    zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) -
                    Dy(sigma[y][y]) -
                    g[y]
                    - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
            Particles.ghost_get<9>(SKIP_LABELLING);


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
            if (Vreset == 1) {
                P_bulk = 0;
                P = 0;
                Vreset = 0;
            }
            P = 0;
            P_bulk = 0;
            while (V_err >= V_err_eps && n <= nmax) {
                RHS[x] = dV[x];
                RHS[y] = dV[y];
                Particles_subset.ghost_get<0>(SKIP_LABELLING);
                Grad_bulk[x] = Bulk_Dx(P_bulk);
                Grad_bulk[y] = Bulk_Dy(P_bulk);
                for (int i = 0; i < bulk.size(); i++) {
                    Particles.template getProp<10>(bulk.template get<0>(i))[x] += Particles_subset.getProp<2>(i)[x];
                    Particles.template getProp<10>(bulk.template get<0>(i))[y] += Particles_subset.getProp<2>(i)[y];
                }
                Particles.ghost_get<10>(SKIP_LABELLING);
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
                //Solver.solve(V[x], V[y]);
                Particles.ghost_get<Velocity>(SKIP_LABELLING);
                div = -(Dx(V[x]) + Dy(V[y]));
                auto Helmholtz = Dxx(H) + Dyy(H);
                DCPSE_scheme<equations2d1, decltype(Particles)> SolverH(Particles);
                SolverH.impose(Helmholtz, bulk, prop_id<19>());
                SolverH.impose(H, up_p1, 0);
                SolverH.impose(H, dw_p1, 0);
                SolverH.impose(H, l_p1, 0);
                SolverH.impose(H, r_p1, 0);
                SolverH.impose(-Dx(H) + Dy(H), corner_ul, 0);
                SolverH.impose(Dx(H) + Dy(H), corner_ur, 0);
                SolverH.impose(-Dx(H) - Dy(H), corner_dl, 0);
                SolverH.impose(Dx(H) - Dy(H), corner_dr, 0);
                SolverH.solve_with_solver(solverPetsc2, H);
                //SolverH.solve(H);
                P = P + div;
                for (int i = 0; i < bulk.size(); i++) {
                    Particles_subset.getProp<0>(i) = Particles.template getProp<4>(bulk.template get<0>(i));
                }
                for (int j = 0; j < up_p.size(); j++) {
                    auto p = up_p.get<0>(j);
                    Particles.getProp<4>(p) = 0;

                }
                for (int j = 0; j < dw_p.size(); j++) {
                    auto p = dw_p.get<0>(j);
                    Particles.getProp<4>(p) = 0;
                }
                for (int j = 0; j < l_p.size(); j++) {
                    auto p = l_p.get<0>(j);
                    Particles.getProp<4>(p) = 0;
                }
                for (int j = 0; j < r_p.size(); j++) {
                    auto p = r_p.get<0>(j);
                    Particles.getProp<4>(p) = 0;
                }
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
                Particles.ghost_get<18>(SKIP_LABELLING);
                V_err_old = V_err;
                V_err = sum / sum1;
                if (V_err > V_err_old) {
                    errctr++;
                    //alpha_P -= 0.1;
                } else {
                    errctr = 0;
                }
                if (n > 3) {
                    if (errctr > 5 || abs(V_err_old - V_err) < 1e-6) {
                        std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/Oscillation IN ERROR" << std::endl;
                        Vreset = 1;
                        break;
                    } else {
                        Vreset = 0;
                    }
                }
                n++;
                //Particles.write_frame("V_debug", n);
                if (v_cl.rank() == 0) {
                    std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;
                }
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

            H_p_b = Pol[x] * Pol[x] + Pol[y] * Pol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }

            h[x] = -gama * (lambda * delmu - nu * (u[x][x] * Pol[x] * Pol[x] + u[y][y] * Pol[y] * Pol[y] +
                                                   2 * u[x][y] * Pol[x] * Pol[y]) / (H_p_b));


            //Particles.ghost_get<MolField, Strain_rate, Vorticity>(SKIP_LABELLING);
            //Particles.write_frame("Polar_withGhost_3e-3", ctr);
            //Particles.deleteGhost();
            Particles.write_frame("Polar_3e-3", ctr);
            Particles.ghost_get<0>(SKIP_LABELLING);

            ctr++;

            H_p_b = sqrt(H_p_b);

            k1[x] = ((h[x] * Pol[x] - h[y] * Pol[y]) / gama + lambda * delmu * Pol[x] -
                     nu * (u[x][x] * Pol[x] + u[x][y] * Pol[y]) + W[x][x] * Pol[x] +
                     W[x][y] * Pol[y]);// - V[x] * Dx(Pol[x]) - V[y] * Dy(Pol[x]));
            k1[y] = ((h[x] * Pol[y] + h[y] * Pol[x]) / gama + lambda * delmu * Pol[y] -
                     nu * (u[y][x] * Pol[x] + u[y][y] * Pol[y]) + W[y][x] * Pol[x] +
                     W[y][y] * Pol[y]);// - V[x] * Dx(Pol[y]) - V[y] * Dy(Pol[y]));

            H_t = H_p_b;//+0.5*dt*(k1[x]*k1[x]+k1[y]*k1[y]);
            dPol = Pol + (0.5 * dt) * k1;
            dPol = dPol / H_t;
            r = dPol[x] * dPol[x] + dPol[y] * dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            Particles.ghost_get<8>(SKIP_LABELLING);


            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] + u[y][y] * dPol[y] * dPol[y] +
                                                    2 * u[x][y] * dPol[x] * dPol[y]) / (r)));

            k2[x] = ((h[x] * (dPol[x]) - h[y] * (dPol[y])) / gama +
                     lambda * delmu * (dPol[x]) -
                     nu * (u[x][x] * (dPol[x]) + u[x][y] * (dPol[y])) +
                     W[x][x] * (dPol[x]) + W[x][y] * (dPol[y])); //-V[x] * Dx((dPol[x])) - V[y] * Dy((dPol[x])));
            k2[y] = ((h[x] * (dPol[y]) + h[y] * (dPol[x])) / gama +
                     lambda * delmu * (dPol[y]) -
                     nu * (u[y][x] * (dPol[x]) + u[y][y] * (dPol[y])) +
                     W[y][x] * (dPol[x]) + W[y][y] * (dPol[y])); //-V[x] * Dx((dPol[y])) - V[y] * Dy((dPol[y])));

            H_t = H_p_b;//+0.5*dt*(k2[x]*k2[x]+k2[y]*k2[y]);
            dPol = Pol + (0.5 * dt) * k2;
            dPol = dPol / H_t;
            r = dPol[x] * dPol[x] + dPol[y] * dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            Particles.ghost_get<8>(SKIP_LABELLING);

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] + u[y][y] * dPol[y] * dPol[y] +
                                                    2 * u[x][y] * dPol[x] * dPol[y]) / (r)));

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
            H_t = H_p_b;//+dt*(k3[x]*k3[x]+k3[y]*k3[y]);
            dPol = Pol + (dt * k3);
            dPol = dPol / H_t;
            Particles.ghost_get<8>(SKIP_LABELLING);
            r = dPol[x] * dPol[x] + dPol[y] * dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));

            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            Particles.ghost_get<8>(SKIP_LABELLING);

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] + u[y][y] * dPol[y] * dPol[y] +
                                                    2 * u[x][y] * dPol[x] * dPol[y]) / (r)));

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

            Pol = Pol + (dt / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
            Pol = Pol / H_p_b;

            H_p_b = sqrt(Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Pol = Pol / H_p_b;
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

            k1 = V;
            k2 = 0.5 * dt * k1 + V;
            k3 = 0.5 * dt * k2 + V;
            k4 = dt * k3 + V;
            //Pos = Pos + dt * V;
            Pos = Pos + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
            Particles.map();
            Particles.ghost_get<0, ExtForce, 27>();
            indexUpdate(Particles, Particles_subset, up_p, dw_p, l_p, r_p, up_p1, dw_p1, l_p1, r_p1, corner_ul,
                        corner_ur, corner_dl, corner_dr, bulk, up, down, left, right);
            Particles_subset.map();
            Particles_subset.ghost_get<0>();

            //Particles_subset.write("debug");

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

            std::cout << "Time step " << ctr - 1 << " : " << tim << " over." << std::endl;
            tim += dt;
            std::cout << "----------------------------------------------------------" << std::endl;
        }
        Particles.deleteGhost();
        tt2.stop();
        std::cout << "The simulation took " << tt2.getwct() << "Seconds.";
    }


    BOOST_AUTO_TEST_CASE(Active2DExp) {
        timer tt2;
        tt2.start();
        double boxsize = 10;
        const size_t sz[2] = {41, 41};
        Box<2, double> box({0, 0}, {boxsize, boxsize});
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.9 * spacing;
        double rCut2 = 3.9 * spacing;
        double rCut3 = 5 * spacing;

        int ord = 2;
        int ord2 = 2;
        double sampling_factor = 4.0;
        double sampling_factor2 = 2.4;
        double sampling_factor3 = 1.6;

        double alpha_V = 1.0;
        double alpha_P = 1.0;
        Ghost<2, double> ghost(rCut3);

        auto &v_cl = create_vcluster();

/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div    H_t                                                                                      delmu */
        vector_dist<2, double, aggregate<VectorS<2, double>, VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double, double, double, double, double, VectorS<2, double>, double, double, double[2], double[2], double[2], double[2], double[2], double[2], double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, double, double, double>> Particles(
                0, box, bc, ghost);
        vector_dist<2, double, aggregate<double, double, VectorS<2, double>>> Particles_subset(
                Particles.getDecomposition(), 0);
        double x0, y0, x1, y1;
        x0 = box.getLow(0);
        y0 = box.getLow(1);
        x1 = box.getHigh(0);
        y1 = box.getHigh(1);

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

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Particles.map();
        Particles.ghost_get<0>();

        Particles.write("Par");

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulk_withoutmid;
        openfpm::vector<aggregate<int>> mid_ref;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;
        openfpm::vector<aggregate<int>> corner_ul;
        openfpm::vector<aggregate<int>> corner_ur;
        openfpm::vector<aggregate<int>> corner_dl;
        openfpm::vector<aggregate<int>> corner_dr;


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
        Particles.ghost_get<ExtForce, 27>(SKIP_LABELLING);


        Box<2, double> up({x0 - spacing / 2.0, y1 - spacing / 2.0},
                          {x1 + spacing / 2.0, y1 + spacing / 2.0});

        Box<2, double> down({x0 - spacing / 2.0, y0 - spacing / 2.0},
                            {x1 + spacing / 2.0, y0 + spacing / 2.0});

        Box<2, double> left({x0 - spacing / 2.0, y0 - spacing / 2.0},
                            {x0 + spacing / 2.0, y1 + spacing / 2.0});

        Box<2, double> right({x1 - spacing / 2.0, y0 - spacing / 2.0},
                             {x1 + spacing / 2.0, y1 + spacing / 2.0});
        Box<2, double> mid({box.getHigh(0) / 2.0 - 0.75*spacing, box.getHigh(1) / 2.0 - 0.75*spacing },
                           {box.getHigh(0) / 2.0+0.75*spacing, box.getHigh(1) / 2.0 + 0.75*spacing});


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
                if (left.isInside(xp) == true) {
                    corner_ul.add();
                    corner_ul.last().get<0>() = p.getKey();
                } else if (right.isInside(xp) == true) {
                    corner_ur.add();
                    corner_ur.last().get<0>() = p.getKey();
                } else {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                }
                up_p.add();
                up_p.last().get<0>() = p.getKey();

            } else if (down.isInside(xp) == true) {
                if (left.isInside(xp) == true) {
                    corner_dl.add();
                    corner_dl.last().get<0>() = p.getKey();
                } else if (right.isInside(xp) == true) {
                    corner_dr.add();
                    corner_dr.last().get<0>() = p.getKey();
                } else {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                }
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                if (up.isInside(xp) == false && down.isInside(xp) == false) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                }
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                if (up.isInside(xp) == false && down.isInside(xp) == false) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                }
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            } else {
                if (mid.isInside(xp) == true) {
                    mid_ref.add();
                    mid_ref.last().get<0>() = p.getKey();
                } else{
                    bulk_withoutmid.add();
                    bulk_withoutmid.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }

        double sigma2 = spacing * spacing / (4);
        std::mt19937 rng{7};
        std::normal_distribution<> gaussian{0, sigma2};


        for(int j=0;j<bulk.size();j++) {
            auto p = bulk.get<0>(j);
            Particles.getPos(p)[0]=Particles.getPos(p)[0]+gaussian(rng);
            Particles.getPos(p)[1]=Particles.getPos(p)[1]+gaussian(rng);
        }


        for (int i = 0; i < bulk.size(); i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
        }


        Particles_subset.map();
        Particles_subset.ghost_get<0>();

        //Particles_subset.write("Pars");
        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dx(Particles_subset, ord,
                                                                                                 rCut,
                                                                                                 sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dy(Particles_subset, ord,
                                                                                                 rCut,
                                                                                                 sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);//, Dxx2(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);//, Dyy2(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);


        Derivative_xxx Dxxx(Particles, ord, rCut3, sampling_factor3,
                          support_options::RADIUS);
        Derivative_xxy Dxxy(Particles, ord, rCut3, sampling_factor3,
                            support_options::RADIUS);

        Derivative_yyx Dyyx(Particles, ord, rCut3, sampling_factor3,
                          support_options::RADIUS);

        Derivative_yyy Dyyy(Particles, ord, rCut3, sampling_factor3,
                            support_options::RADIUS);

/*        Derivative_x Dx(Particles, ord, rCut, sampling_factor), Bulk_Dx(Particles_subset, ord, rCut, sampling_factor);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor), Bulk_Dy(Particles_subset, ord, rCut, sampling_factor);
        Derivative_xy Dxy(Particles, ord2, rCut, sampling_factor);
        auto Dyx = Dxy;
        Derivative_xx Dxx(Particles, ord2, rCut2, sampling_factor2),Bulk_Dxx(Particles_subset, ord2, rCut2, sampling_factor2);
        Derivative_yy Dyy(Particles, ord2, rCut2, sampling_factor2),Bulk_Dyy(Particles_subset, ord2, rCut2, sampling_factor2);*/


        V.setVarId(0);
        P.setVarId(2);
        H.setVarId(3);
        eq_id vx, vy, ic,helm;
        timer tt;
        timer tt3;
        vx.setId(0);
        vy.setId(1);
        ic.setId(2);
        helm.setId(3);
        double V_err_eps = 5 * 1e-2;
        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 300;
        int ctr = 0, errctr, Vreset = 0;
        double dt = 3 * 1e-3;

        double tim = 0;
        double tf = 1.25;
        div = 0;
        double sum, sum1, sum_k;
        while (tim <= tf) {
            tt.start();
            /*petsc_solver<double> solverPetsc;
            solverPetsc.setSolver(KSPGMRES);
            //solverPetsc.setRestart(250);
            solverPetsc.setPreconditioner(PCJACOBI);*/

            Particles.ghost_get<Polarization>(SKIP_LABELLING);
            sigma[x][x] =
                    -Ks * Dx(Pol[x]) * Dx(Pol[x]) - Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
            sigma[x][y] =
                    -Ks * Dy(Pol[y]) * Dx(Pol[y]) - Kb * Dy(Pol[x]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dx(Pol[x]);
            sigma[y][x] =
                    -Ks * Dx(Pol[x]) * Dy(Pol[x]) - Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
            sigma[y][y] =
                    -Ks * Dy(Pol[y]) * Dy(Pol[y]) - Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dx(Pol[y]) * Dy(Pol[x]);
            Particles.ghost_get<Stress>(SKIP_LABELLING);


            r = Pol[x] * Pol[x] + Pol[y] * Pol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }

            h[y] = (Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) -
                    Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y])));

            Particles.ghost_get<MolField>(SKIP_LABELLING);

            FranckEnergyDensity = (Ks / 2.0) *
                                  ((Dx(Pol[x]) * Dx(Pol[x])) + (Dy(Pol[x]) * Dy(Pol[x])) +
                                   (Dx(Pol[y]) * Dx(Pol[y])) +
                                   (Dy(Pol[y]) * Dy(Pol[y]))) +
                                  ((Kb - Ks) / 2.0) * ((Dx(Pol[y]) - Dy(Pol[x])) * (Dx(Pol[y]) - Dy(Pol[x])));
            Particles.ghost_get<33>(SKIP_LABELLING);


            f1 = gama * nu * Pol[x] * Pol[x] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
            f2 = 2.0 * gama * nu * Pol[x] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
            f3 = gama * nu * Pol[y] * Pol[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y]) / (r);
            f4 = 2.0 * gama * nu * Pol[x] * Pol[x] * Pol[x] * Pol[y] / (r);
            f5 = 4.0 * gama * nu * Pol[x] * Pol[x] * Pol[y] * Pol[y] / (r);
            f6 = 2.0 * gama * nu * Pol[x] * Pol[y] * Pol[y] * Pol[y] / (r);
            Particles.ghost_get<11, 12, 13, 14, 15, 16>(SKIP_LABELLING);
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
            Particles.ghost_get<21, 22, 23, 24, 25, 26>(SKIP_LABELLING);


            dV[x] = -0.5 * Dy(h[y]) + zeta * Dx(delmu * Pol[x] * Pol[x]) + zeta * Dy(delmu * Pol[x] * Pol[y]) -
                    zeta * Dx(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dx(-2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dy(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[x][x]) -
                    Dy(sigma[x][y]) -
                    g[x]
                    - 0.5 * nu * Dx(-gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dy(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));


            dV[y] = -0.5 * Dx(-h[y]) + zeta * Dy(delmu * Pol[y] * Pol[y]) + zeta * Dx(delmu * Pol[x] * Pol[y]) -
                    zeta * Dy(0.5 * delmu * (Pol[x] * Pol[x] + Pol[y] * Pol[y])) -
                    0.5 * nu * Dy(2.0 * h[y] * Pol[x] * Pol[y])
                    - 0.5 * nu * Dx(h[y] * (Pol[x] * Pol[x] - Pol[y] * Pol[y])) - Dx(sigma[y][x]) -
                    Dy(sigma[y][y]) -
                    g[y]
                    - 0.5 * nu * Dy(gama * lambda * delmu * (Pol[x] * Pol[x] - Pol[y] * Pol[y]))
                    - 0.5 * Dx(-2.0 * gama * lambda * delmu * (Pol[x] * Pol[y]));
            Particles.ghost_get<9>(SKIP_LABELLING);


            //Particles.write("PolarI");
            //Velocity Solution n iterations


/*            auto Stokes1 = -(Dxxx(H) + Dyyx(H)) + eta * (Dxx(V[x]) + Dyy(V[x]))
                           + 0.5 * nu * (Df1[x] * Dx(V[x]) + f1 * Dxx(V[x]))
                           + 0.5 * nu * (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y]))
                           + 0.5 * nu * (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x]))
                           + 0.5 * nu * (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           + 0.5 * nu * (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
            auto Stokes2 = -(Dxxy(H) + Dyyy(H)) + eta * (Dxx(V[y]) + Dyy(V[y]))
                           - 0.5 * nu * (Df1[y] * Dx(V[x]) + f1 * Dxy(V[x]))
                           - 0.5 * nu * (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           - 0.5 * nu * (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y]))
                           + 0.5 * nu * (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x]))
                           + 0.5 * nu * (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));

            auto Helmholtz = Dxx(H) + Dyy(H)+Dx(Stokes1) + Dy(Stokes2);*/

            auto Stokes1 = -(Bulk_Dx(P)) + eta * (Dxx(V[x]) + Dyy(V[x]))
                           + 0.5 * nu * (Df1[x] * Dx(V[x]) + f1 * Dxx(V[x]))
                           + 0.5 * nu * (Df2[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df3[x] * Dy(V[y]) + f3 * Dyx(V[y]))
                           + 0.5 * nu * (Df4[y] * Dx(V[x]) + f4 * Dxy(V[x]))
                           + 0.5 * nu * (Df5[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           + 0.5 * nu * (Df6[y] * Dy(V[y]) + f6 * Dyy(V[y]));
            auto Stokes2 = -(Bulk_Dx(P)) + eta * (Dxx(V[y]) + Dyy(V[y]))
                           - 0.5 * nu * (Df1[y] * Dx(V[x]) + f1 * Dxy(V[x]))
                           - 0.5 * nu * (Df2[y] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f2 * 0.5 * (Dxy(V[y]) + Dyy(V[x])))
                           - 0.5 * nu * (Df3[y] * Dy(V[y]) + f3 * Dyy(V[y]))
                           + 0.5 * nu * (Df4[x] * Dx(V[x]) + f4 * Dxx(V[x]))
                           + 0.5 * nu * (Df5[x] * 0.5 * (Dx(V[y]) + Dy(V[x])) + f5 * 0.5 * (Dxx(V[y]) + Dyx(V[x])))
                           + 0.5 * nu * (Df6[x] * Dy(V[y]) + f6 * Dyx(V[y]));

           auto incompressibility=(Dx(V[x]) + Dy(V[y]));

            tt.stop();
            std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
            tt.start();
            RHS[x] = dV[x];
            RHS[y] = dV[y];
            div=-(Dx(dV[x])+Dy(dV[y]));
            Particles.ghost_get<10>(SKIP_LABELLING);
            DCPSE_scheme<equations2d3E, decltype(Particles)> Solver(Particles, 41+41+39+39);
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
            /*Solver.impose(Helmholtz, bulk, prop_id<19>(), helm);
            Solver.impose(H, up_p1, 0, helm);
            Solver.impose(H, dw_p1, 0, helm);
            Solver.impose(H, l_p1, 0, helm);
            Solver.impose(H, r_p1, 0, helm);
            Solver.impose(-Dx(H) + Dy(H), corner_ul, 0, helm);
            Solver.impose(Dx(H) + Dy(H), corner_ur, 0, helm);
            Solver.impose(-Dx(H) - Dy(H), corner_dl, 0, helm);
            Solver.impose(Dx(H) - Dy(H), corner_dr, 0, helm);
             Solver.impose(P+Dxx(H) + Dyy(H), mid_ref, 0, ic);*/
            Solver.impose(incompressibility, bulk_withoutmid, 0, ic);
           /* Solver.impose(incompressibility, up_p, 0, ic);
            Solver.impose(incompressibility, dw_p, 0, ic);
            Solver.impose(incompressibility, l_p, 0, ic);
            Solver.impose(incompressibility, r_p, 0, ic);*/
            Solver.impose(P, mid_ref, 0, ic);
            Solver.solve(V[x], V[y], P);
            //Solver.solve(V[x], V[y]);
            Particles.ghost_get<Velocity>(SKIP_LABELLING);
            //SolverH.solve(H);
            //P = (Dxx(H) + Dyy(H));
            for (int i = 0; i < bulk.size(); i++) {
                Particles_subset.getProp<0>(i) = Particles.template getProp<4>(bulk.template get<0>(i));
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<4>(p) = 0;

            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<4>(p) = 0;
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<4>(p) = 0;
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<4>(p) = 0;
            }
            tt.stop();
            std::cout << "Velocity Solved" << std::endl;
            Particles.write("V_DEBUG");
            return;
            u[x][x] = Dx(V[x]);
            u[x][y] = 0.5 * (Dx(V[y]) + Dy(V[x]));
            u[y][x] = 0.5 * (Dy(V[x]) + Dx(V[y]));
            u[y][y] = Dy(V[y]);

            W[x][x] = 0;
            W[x][y] = 0.5 * (Dy(V[x]) - Dx(V[y]));
            W[y][x] = 0.5 * (Dx(V[y]) - Dy(V[x]));
            W[y][y] = 0;

            H_p_b = Pol[x] * Pol[x] + Pol[y] * Pol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<32>(p) = (Particles.getProp<32>(p) == 0) ? 1 : Particles.getProp<32>(p);
            }

            h[x] = -gama * (lambda * delmu - nu * (u[x][x] * Pol[x] * Pol[x] + u[y][y] * Pol[y] * Pol[y] +
                                                   2 * u[x][y] * Pol[x] * Pol[y]) / (H_p_b));


            //Particles.ghost_get<MolField, Strain_rate, Vorticity>(SKIP_LABELLING);
            //Particles.write_frame("Polar_withGhost_3e-3", ctr);
            Particles.deleteGhost();
            Particles.write_frame("Polar_3e-3", ctr);
            Particles.ghost_get<0>();

            ctr++;

            H_p_b = sqrt(H_p_b);

            k1[x] = ((h[x] * Pol[x] - h[y] * Pol[y]) / gama + lambda * delmu * Pol[x] -
                     nu * (u[x][x] * Pol[x] + u[x][y] * Pol[y]) + W[x][x] * Pol[x] +
                     W[x][y] * Pol[y]);// - V[x] * Dx(Pol[x]) - V[y] * Dy(Pol[x]));
            k1[y] = ((h[x] * Pol[y] + h[y] * Pol[x]) / gama + lambda * delmu * Pol[y] -
                     nu * (u[y][x] * Pol[x] + u[y][y] * Pol[y]) + W[y][x] * Pol[x] +
                     W[y][y] * Pol[y]);// - V[x] * Dx(Pol[y]) - V[y] * Dy(Pol[y]));

            H_t = H_p_b;//+0.5*dt*(k1[x]*k1[x]+k1[y]*k1[y]);
            dPol = Pol + (0.5 * dt) * k1;
            dPol = dPol / H_t;
            r = dPol[x] * dPol[x] + dPol[y] * dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            Particles.ghost_get<8>(SKIP_LABELLING);


            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] + u[y][y] * dPol[y] * dPol[y] +
                                                    2 * u[x][y] * dPol[x] * dPol[y]) / (r)));

            k2[x] = ((h[x] * (dPol[x]) - h[y] * (dPol[y])) / gama +
                     lambda * delmu * (dPol[x]) -
                     nu * (u[x][x] * (dPol[x]) + u[x][y] * (dPol[y])) +
                     W[x][x] * (dPol[x]) + W[x][y] * (dPol[y])); //-V[x] * Dx((dPol[x])) - V[y] * Dy((dPol[x])));
            k2[y] = ((h[x] * (dPol[y]) + h[y] * (dPol[x])) / gama +
                     lambda * delmu * (dPol[y]) -
                     nu * (u[y][x] * (dPol[x]) + u[y][y] * (dPol[y])) +
                     W[y][x] * (dPol[x]) + W[y][y] * (dPol[y])); //-V[x] * Dx((dPol[y])) - V[y] * Dy((dPol[y])));

            H_t = H_p_b;//+0.5*dt*(k2[x]*k2[x]+k2[y]*k2[y]);
            dPol = Pol + (0.5 * dt) * k2;
            dPol = dPol / H_t;
            r = dPol[x] * dPol[x] + dPol[y] * dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            Particles.ghost_get<8>(SKIP_LABELLING);

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] + u[y][y] * dPol[y] * dPol[y] +
                                                    2 * u[x][y] * dPol[x] * dPol[y]) / (r)));

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
            H_t = H_p_b;//+dt*(k3[x]*k3[x]+k3[y]*k3[y]);
            dPol = Pol + (dt * k3);
            dPol = dPol / H_t;
            Particles.ghost_get<8>(SKIP_LABELLING);
            r = dPol[x] * dPol[x] + dPol[y] * dPol[y];
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
            }
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));

            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<34>(p) = (Particles.getProp<34>(p) == 0) ? 1 : Particles.getProp<34>(p);
                Particles.getProp<8>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
            }
            Particles.ghost_get<8>(SKIP_LABELLING);

            h[y] = (dPol[x] * (Ks * Dyy(dPol[y]) + Kb * Dxx(dPol[y]) + (Ks - Kb) * Dxy(dPol[x])) -
                    dPol[y] * (Ks * Dxx(dPol[x]) + Kb * Dyy(dPol[x]) + (Ks - Kb) * Dxy(dPol[y])));

            h[x] = -gama * (lambda * delmu - nu * ((u[x][x] * dPol[x] * dPol[x] + u[y][y] * dPol[y] * dPol[y] +
                                                    2 * u[x][y] * dPol[x] * dPol[y]) / (r)));

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

            Pol = Pol + (dt / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
            Pol = Pol / H_p_b;

            H_p_b = sqrt(Pol[x] * Pol[x] + Pol[y] * Pol[y]);
            Pol = Pol / H_p_b;
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

            k1 = V;
            k2 = 0.5 * dt * k1 + V;
            k3 = 0.5 * dt * k2 + V;
            k4 = dt * k3 + V;
            //Pos = Pos + dt * V;
            Pos = Pos + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
            Particles.map();
            Particles.ghost_get<0, ExtForce, 27>();
            indexUpdate(Particles, Particles_subset, up_p, dw_p, l_p, r_p, up_p1, dw_p1, l_p1, r_p1, corner_ul,
                        corner_ur, corner_dl, corner_dr, bulk, up, down, left, right);
            Particles_subset.map();
            Particles_subset.ghost_get<0>();

            //Particles_subset.write("debug");

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

            std::cout << "Time step " << ctr - 1 << " : " << tim << " over." << std::endl;
            tim += dt;
            std::cout << "----------------------------------------------------------" << std::endl;
        }
        Particles.deleteGhost();
        tt2.stop();
        std::cout << "The simulation took " << tt2.getwct() << "Seconds.";
    }

BOOST_AUTO_TEST_SUITE_END()
