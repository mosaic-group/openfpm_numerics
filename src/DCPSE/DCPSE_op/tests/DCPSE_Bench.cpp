//
// Created by Abhinav Singh on 27.10.20.
//

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

#include "config.h"

#define BOOST_TEST_DYN_LINK


#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <iostream>
#include "../DCPSE_op.hpp"
#include "../DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "../EqnsStruct.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests_paper)

    BOOST_AUTO_TEST_CASE(Lid_cavity_paper) {
        /*Data from Ghia et al.
         * y,u
        1.0000,1.00000
        0.9766,0.84123
        0.9688,0.78871
        0.9609,0.73722
        0.9531,0.68717
        0.8516,0.23151
        0.7344,0.00332
        0.6172,-0.13641
        0.5000,-0.20581
        0.4531,-0.21090
        0.2813,-0.15662
        0.1719,-0.10150
        0.1016,-0.06434
        0.0703,-0.04775
        0.0625,-0.04192
        0.0547,-0.03717
        0.0000,0.00000*/
        constexpr int x = 0;
        constexpr int y = 1;
        size_t gd_sz = 256;
        const size_t sz[2] = {gd_sz,gd_sz};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing;
        spacing = 1.0 / (sz[0] - 1);
        double rCut = 4.1 * spacing;
        int ord = 2;
        double sampling_factor = 4.0;
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        auto &v_cl = create_vcluster();
        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,double>> Particles(0, box,
                                                                                                                                                                           bc,
                                                                                                                                                                           ghost);
        //Init_DCPSE(Particles)
        BOOST_TEST_MESSAGE("Init Particles...");

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


        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing;
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing;
            Particles.getLastPos()[1] = yp0;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync Particles across processors...");
        Particles.map();
        Particles.ghost_get<0>();

        double x0=box.getLow(0),x1=box.getHigh(0),y0=box.getLow(1),y1=box.getHigh(1);

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

        vector_dist_subset<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>,double,VectorS<2, double>,double,double>> Particles_bulk(Particles,bulk);

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto RHS = getV<2>(Particles);
        auto dV = getV<3>(Particles);
        auto div = getV<4>(Particles);
        auto V_star = getV<5>(Particles);
        auto H = getV<6>(Particles);
        auto dP = getV<7>(Particles);



        auto P_bulk = getV<0>(Particles_bulk);
        auto V_bulk = getV<1>(Particles_bulk);
        auto RHS_bulk =getV<2>(Particles_bulk);
        auto V_star_bulk = getV<5>(Particles_bulk);
        auto dP_bulk = getV<7>(Particles_bulk);


        P_bulk = 0;
        Advection Adv(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Dx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        //Laplacian Lap(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_xx Dxx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);

        Derivative_y Dy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Bulk_Dx(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Bulk_Dy(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);
        Laplacian Bulk_Lap(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);

        int n = 0, nmax = 50, ctr = 0, errctr=1, Vreset = 0;
        double V_err=1;
        if (Vreset == 1) {
            P_bulk = 0;
            P = 0;
            Vreset = 0;
        }
        P=0;
        eq_id vx,vy;
        vx.setId(0);
        vy.setId(1);
        double Re=100;
        double sum, sum1, sum_k,V_err_eps=1e-3,V_err_old;
        auto StokesX=V[x]*Dx(V_star[x])-(1.0/Re)*(Dxx(V_star[x])+Dyy(V_star[x]));
        auto StokesY=V[y]*Dy(V_star[y])-(1.0/Re)*(Dxx(V_star[y])+Dyy(V_star[y]));
/*
        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc.setPreconditioner(PCJACOBI);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc2.setPreconditioner(PCJACOBI);*/
        RHS[x] = V[x];
        RHS[y] = V[y];
        dV=0;
        timer tt;
        while (V_err >= V_err_eps && n <= nmax) {
            Particles.write_frame("LID",n);
            Particles.ghost_get<0>(SKIP_LABELLING);
            RHS_bulk[x] = -Bulk_Dx(P);
            RHS_bulk[y] = -Bulk_Dy(P);
            DCPSE_scheme<equations2d2E, decltype(Particles)> Solver(Particles);
            Solver.impose(StokesX, bulk, RHS[0], vx);
            Solver.impose(StokesY, bulk, RHS[1], vy);
            Solver.impose(V_star[x], up_p1, 1.0, vx);
            Solver.impose(V_star[y], up_p1, 0.0, vy);
            Solver.impose(V_star[x], l_p, 0.0, vx);
            Solver.impose(V_star[y], l_p, 0.0, vy);
            Solver.impose(V_star[x], r_p, 0.0, vx);
            Solver.impose(V_star[y], r_p, 0.0, vy);
            Solver.impose(V_star[x], dw_p, 0.0, vx);
            Solver.impose(V_star[y], dw_p, 0.0, vy);
            Solver.impose(V_star[x], corner_ul, 0.0, vx);
            Solver.impose(V_star[y], corner_ul, 0.0, vy);
            Solver.impose(V_star[x], corner_ur, 0.0, vx);
            Solver.impose(V_star[y], corner_ur, 0.0, vy);
            //Solver.solve_with_solver(solverPetsc, V_star[x], V_star[y]);
            Solver.solve(V_star[x], V_star[y]);
            //Particles.write_frame("LIDV",n);
            Particles.ghost_get<5>(SKIP_LABELLING);
            div = (Dx(V_star[x]) + Dy(V_star[y]));
            DCPSE_scheme<equations2d1E,decltype(Particles)> SolverH(Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Dxx(H)+Dyy(H);
            SolverH.impose(Helmholtz,bulk,prop_id<4>());
/*            SolverH.impose(H, up_p,0);
            SolverH.impose(H, l_p,0);
            SolverH.impose(H, r_p,0);
            SolverH.impose(H, dw_p,0);*/
            SolverH.impose(Dy(H), up_p1,0);
             SolverH.impose(Dx(H), l_p1,0);
             SolverH.impose(Dx(H), r_p1,0);
             SolverH.impose(Dy(H), dw_p1,0);
             SolverH.impose(-Dx(H) + Dy(H), corner_ul, 0);
             SolverH.impose(Dx(H) + Dy(H), corner_ur, 0);
             SolverH.impose(-Dx(H) - Dy(H), corner_dl, 0);
             SolverH.impose(Dx(H) - Dy(H), corner_dr, 0);
            //SolverH.solve_with_solver(solverPetsc2,H);
            SolverH.solve(H);
            Particles.ghost_get<6>(SKIP_LABELLING);
            //dP_bulk=Bulk_Lap(H);
            P_bulk = P - 0.0125*div;
            //dV[x]=Dx(H);
            //dV[y]=Dy(H);
            V_star_bulk[0] = V_star[0] - Bulk_Dx(H);
            V_star_bulk[1] = V_star[1] - Bulk_Dy(H);
            sum = 0;
            sum1 = 0;

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) *
                       (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) +
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]) *
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]);
                sum1 += Particles.getProp<5>(p)[0] * Particles.getProp<5>(p)[0] +
                        Particles.getProp<5>(p)[1] * Particles.getProp<5>(p)[1];
            }

            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V = V_star;
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            V_err_old = V_err;
            V_err = sum / sum1;
/*            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
                errctr++;
                //alpha_P -= 0.1;
            } else {
                errctr = 0;
            }
            if (n > 3) {
                if (errctr > 3) {
                    std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN ERROR" << std::endl;
                    Vreset = 1;
                    break;
                } else {
                    Vreset = 0;
                }
            }*/
            n++;
            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;
            }
        }
        double worst1 = 0.0;
        double worst2 = 0.0;

      /*  for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(Particles.getProp<6>(p)[0] - Particles.getProp<1>(p)[0]) >= worst1) {
                worst1 = fabs(Particles.getProp<6>(p)[0] - Particles.getProp<1>(p)[0]);
            }
        }
        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(Particles.getProp<6>(p)[1] - Particles.getProp<1>(p)[1]) >= worst2) {
                worst2 = fabs(Particles.getProp<6>(p)[1] - Particles.getProp<1>(p)[1]);
            }
        }*/
        //Particles.deleteGhost();
        //Particles.write("PC_subset_lid");
/*        std::cout << "Maximum Analytic Error in Vx: " << worst1 << std::endl;
        std::cout << "Maximum Analytic Error in Vy: " << worst2 << std::endl;*/
        Particles.deleteGhost();
        Particles.write("LID");

    }
    BOOST_AUTO_TEST_CASE(Active2d_paper) {

    }
    BOOST_AUTO_TEST_CASE(Sphere3d_paper) {

    }



BOOST_AUTO_TEST_SUITE_END()

