/*
 * DCPSE_op_test.cpp
 *
 *  Created on: May 15, 2020
 *      Author: Abhinav Singh
 *
 */
#include "config.h"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../DCPSE_op.hpp"
#include "../DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "../EqnsStruct.hpp"

//template<typename T>
//struct Debug;



BOOST_AUTO_TEST_SUITE(dcpse_op_subset_suite_tests)

    BOOST_AUTO_TEST_CASE(dcpse_op_subset_tests) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {1.0,1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        int ord = 2;
        double sampling_factor = 4.0;
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, double, double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double>> Particles(0, box,
                                                                                                                 bc,
                                                                                                                 ghost);

        //Init_DCPSE(Particles)
        BOOST_TEST_MESSAGE("Init Particles...");
        std::mt19937 rng{6666666};

        std::normal_distribution<> gaussian{0, sigma2};

        openfpm::vector<aggregate<int>> bulk;

        auto it = Particles.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing[0];
            Particles.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            Particles.getLastPos()[1] = y;//+gaussian(rng);
            // Here fill the function value
            Particles.template getLastProp<0>() = sin(Particles.getLastPos()[0]) + sin(Particles.getLastPos()[1]);

            if (k0 != 0 && k1 != 0 && k0 != sz[0] -1 && k1 != sz[1] - 1)
            {
            	bulk.add();
            	bulk.template get<0>(bulk.size()-1) = Particles.size_local() - 1;
            }

            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync Particles across processors...");

        Particles.map();
        Particles.ghost_get<0>();

        vector_dist_subset<2, double, aggregate<double, double, double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double>> Particles_bulk(Particles,bulk);

        // move particles

        auto P = getV<0>(Particles);
        auto Out = getV<1>(Particles);
        auto Pb = getV<2>(Particles);
        auto Out_V = getV<3>(Particles);

        auto P_bulk = getV<2>(Particles_bulk);
        auto Out_bulk = getV<1>(Particles_bulk);
	    auto Out_V_bulk = getV<3>(Particles_bulk);
        Out=10;
        P_bulk = 5;

        P_bulk=Pb+Out;
//        Particles.write("Test_output_subset");

        // Create the subset

       /* Derivative_x Dx(Particles, 2, rCut);
        Derivative_y Dy(Particles, 2, rCut);
        Derivative_x Dx_bulk(Particles_bulk, 2, rCut);
*/

        Derivative_x Dx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Dx_bulk(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Dy_bulk(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);

        Out_bulk = Dx_bulk(P);
	    Out_V_bulk[0] = P + Dx_bulk(P);
        Out_V_bulk[1] = Out_V[0] +Dy_bulk(P);

	// Check
	auto it2 = Particles_bulk.getDomainIterator();
        while (it2.isNext())
        {
		auto p = it2.get();

		BOOST_REQUIRE_EQUAL(Particles_bulk.getProp<2>(p),15.0);
		BOOST_REQUIRE(fabs(Particles_bulk.getProp<1>(p) - cos(Particles_bulk.getPos(p)[0])) < 0.005 );
		BOOST_REQUIRE(fabs(Particles_bulk.getProp<3>(p)[0] - Particles_bulk.getProp<0>(p) - cos(Particles_bulk.getPos(p)[0])) < 0.001 );
		BOOST_REQUIRE(fabs(Particles_bulk.getProp<3>(p)[1] - Particles_bulk.getProp<3>(p)[0] - cos(Particles_bulk.getPos(p)[1])) < 0.001 );

            ++it2;
	}

//        P_bulk = Dx_bulk(P_bulk);  <------------ Incorrect produce error message
//        P = Dx_bulk(P);   <------- Incorrect produce overflow

        Particles.write("Out");
    }


    BOOST_AUTO_TEST_CASE(dcpse_op_subset_PC_lid) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        constexpr int x = 0;
        constexpr int y = 1;
        size_t edgeSemiSize = 20;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {1.0,1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        int ord = 2;
        double sampling_factor = 4.0;
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);
        auto &v_cl = create_vcluster();


        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,VectorS<2, double>,double>> Particles(0, box,
                                                                                                                 bc,
                                                                                                                 ghost);

        //Init_DCPSE(Particles)
        BOOST_TEST_MESSAGE("Init Particles...");

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> boundary;

        auto it = Particles.getGridIterator(sz);
        size_t pointId = 0;
        double minNormOne = 999;
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync Particles across processors...");
        Particles.map();
        Particles.ghost_get<0>();

        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (xp[0] != 0 && xp[1] != 0 && xp[0] != 1.0 && xp[1] != 1.0) {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<3>(p)[x] = 3.0;
                Particles.getProp<3>(p)[y] = 3.0;
            } else {
                boundary.add();
                boundary.last().get<0>() = p.getKey();
                Particles.getProp<3>(p)[x] = xp[0]*xp[0]+xp[1]*xp[1];
                Particles.getProp<3>(p)[y] = xp[0]*xp[0]-2*xp[0]*xp[1];
            }
            Particles.getProp<6>(p)[x] = xp[0]*xp[0]+xp[1]*xp[1];
            Particles.getProp<6>(p)[y] = xp[0]*xp[0]-2*xp[0]*xp[1];
            Particles.getProp<7>(p) = xp[0]+xp[1]-1.0;

            ++it2;
        }

        vector_dist_subset<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>,double,VectorS<2, double>,VectorS<2, double>,double>> Particles_bulk(Particles,bulk);

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto RHS = getV<2>(Particles);
        auto dV = getV<3>(Particles);
        auto div = getV<4>(Particles);
        auto V_star = getV<5>(Particles);


        auto P_bulk = getV<0>(Particles_bulk);
        auto RHS_bulk =getV<2>(Particles_bulk);

        P_bulk = 0;

        Derivative_x Dx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_xx Dxx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Bulk_Dx(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Bulk_Dy(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);

        int n = 0, nmax = 5, ctr = 0, errctr=1, Vreset = 0;
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
        double sum, sum1, sum_k,V_err_eps=1e-3,V_err_old;
        auto Stokes1=Dxx(V[x])+Dyy(V[x]);
        auto Stokes2=Dxx(V[y])+Dyy(V[y]);

        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc.setPreconditioner(PCJACOBI);
        V_star=0;
        RHS[x] = dV[x];
        RHS[y] = dV[y];
        while (V_err >= V_err_eps && n <= nmax) {
            RHS_bulk[x] = dV[x] + Bulk_Dx(P);
            RHS_bulk[y] = dV[y] + Bulk_Dy(P);

            DCPSE_scheme<equations2d2, decltype(Particles)> Solver(Particles);
            Solver.impose(Stokes1, bulk, RHS[0], vx);
            Solver.impose(Stokes2, bulk, RHS[1], vy);
            Solver.impose(V[x], boundary, RHS[0], vx);
            Solver.impose(V[y], boundary, RHS[1], vy);
            if (n == 0)
            {
            	auto & A=Solver.getA(options_solver::STANDARD);
            	A.write("Mat_lid");
            	auto & b = Solver.getB();
            	b.write("vect_lid");
            }
            Solver.solve_with_solver(solverPetsc, V[x], V[y]);

            if (n == 0)
            {
            	Particles.write_frame("PC_subset_lid",n);
            }

            Particles.ghost_get<0>(SKIP_LABELLING);
            div = -(Dx(V[x]) + Dy(V[y]));
            P_bulk = P + div;
            //Particles.write_frame("PC_subset_lid",n);
            //P_bulk = P_bulk + div_bulk;
            sum = 0;
            sum1 = 0;

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) *
                       (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) +
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]) *
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]);
                sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                        Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1];
            }

            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_star = V;
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            V_err_old = V_err;
            V_err = sum / sum1;
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
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
            }
            n++;
            //if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;
            //}
        }

        Particles.write("PC_subset_lid");
    }


    BOOST_AUTO_TEST_CASE(dcpse_op_subset_PC_lid2) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        constexpr int x = 0;
        constexpr int y = 1;
        size_t edgeSemiSize = 20;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {1.0,1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        int ord = 2;
        double sampling_factor = 4.0;
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        auto &v_cl = create_vcluster();

        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,VectorS<2, double>,double>> Particles(0, box,
                                                                                                                                                 bc,
                                                                                                                                                 ghost);
        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,VectorS<2, double>,double>> Particles_subset(Particles.getDecomposition(), 0);


        //Init_DCPSE(Particles)
        BOOST_TEST_MESSAGE("Init Particles...");

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> boundary;

        auto it = Particles.getGridIterator(sz);
        size_t pointId = 0;
        double minNormOne = 999;
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync Particles across processors...");
        Particles.map();
        Particles.ghost_get<0>();
        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (xp[0] != 0 && xp[1] != 0 && xp[0] != 1.0 && xp[1] != 1.0) {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<3>(p)[x] = 3.0;
                Particles.getProp<3>(p)[y] = 3.0;
            } else {
                boundary.add();
                boundary.last().get<0>() = p.getKey();
                Particles.getProp<3>(p)[x] = xp[0]*xp[0]+xp[1]*xp[1];
                Particles.getProp<3>(p)[y] = xp[0]*xp[0]-2*xp[0]*xp[1];
            }
            Particles.getProp<6>(p)[x] = xp[0]*xp[0]+xp[1]*xp[1];
            Particles.getProp<6>(p)[y] = xp[0]*xp[0]-2*xp[0]*xp[1];
            Particles.getProp<7>(p) = xp[0]+xp[1]-1.0;

            ++it2;
        }

        for (int i = 0; i < bulk.size(); i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
        }
        Particles_subset.map();
        Particles_subset.ghost_get<0>();



        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto RHS = getV<2>(Particles);
        auto dV = getV<3>(Particles);
        auto div = getV<4>(Particles);
        auto V_star = getV<5>(Particles);

        auto P_bulk = getV<0>(Particles_subset);
        auto Grad_bulk= getV<2>(Particles_subset);

        P_bulk = 0;

        Derivative_x Dx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Bulk_Dx(Particles_subset, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_xx Dxx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,sampling_factor, support_options::RADIUS),Bulk_Dy(Particles_subset, 2, rCut,sampling_factor, support_options::RADIUS);;

        int n = 0, nmax = 5, ctr = 0, errctr=0, Vreset = 0;
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
        double sum, sum1, sum_k,V_err_eps=1e-3,V_err_old;
        auto Stokes1=Dxx(V[x])+Dyy(V[x]);
        auto Stokes2=Dxx(V[y])+Dyy(V[y]);

        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc.setPreconditioner(PCJACOBI);
        V_star=0;
        while (V_err >= V_err_eps && n <= nmax) {
            RHS[x] = dV[x];
            RHS[y] = dV[y];
            Particles_subset.ghost_get<0>(SKIP_LABELLING);

            Grad_bulk[x] = Bulk_Dx(P_bulk);
            Grad_bulk[y] = Bulk_Dy(P_bulk);
            for (int i = 0; i < bulk.size(); i++) {
                Particles.template getProp<2>(bulk.template get<0>(i))[x] += Particles_subset.getProp<2>(i)[x];
                Particles.template getProp<2>(bulk.template get<0>(i))[y] += Particles_subset.getProp<2>(i)[y];
            }

            DCPSE_scheme<equations2d2, decltype(Particles)> Solver(Particles);
            Solver.impose(Stokes1, bulk, RHS[0], vx);
            Solver.impose(Stokes2, bulk, RHS[1], vy);
            Solver.impose(V[x], boundary, RHS[0], vx);
            Solver.impose(V[y], boundary, RHS[1], vy);
            Solver.solve_with_solver(solverPetsc, V[x], V[y]);

            Particles.ghost_get<0>(SKIP_LABELLING);
            div = -(Dx(V[x]) + Dy(V[y]));
            P = P + div;
            for (int i = 0; i < bulk.size(); i++) {
                Particles_subset.getProp<0>(i) = Particles.template getProp<0>(bulk.template get<0>(i));
            }
            sum = 0;
            sum1 = 0;

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) *
                       (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) +
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]) *
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]);
                sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                        Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1];
            }

            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_star=V;
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            V_err_old = V_err;
            V_err = sum / sum1;
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
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
            }
            n++;
            //if (v_cl.rank() == 0) {
            std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;

            //}
        }

        Particles.write("PC_subset_lid2");
    }


BOOST_AUTO_TEST_SUITE_END()
