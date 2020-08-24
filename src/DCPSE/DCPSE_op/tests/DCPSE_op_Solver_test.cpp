/*
 * DCPSE_op_Solver_test.cpp
 *
 *  Created on: Jan 7, 2020
 *      Author: Abhinav Singh, Pietro Incardona
 *
 */
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 30
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



BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)

    BOOST_AUTO_TEST_CASE(dcpse_op_solver) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31, 31};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Laplacian Lap(domain, 2, rCut, 2);

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) - spacing / 2.0, box.getHigh(1) - 8*spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - 6*spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) - spacing / 2.0, box.getLow(1) + 3*spacing / 2.0},
                              {box.getHigh(0) + spacing / 2.0, box.getLow(1) + 4*spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                               {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(up_d);
        boxes.add(down);
        boxes.add(down_u);
        boxes.add(left);
        boxes.add(left_r);
        boxes.add(right);
        boxes.add(right_l);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            domain.getProp<2>(p)=1+xp[0]*xp[0]+2*xp[1]*xp[1];
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  3 + xp.get(0)*xp.get(0);
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  1 + xp.get(0)*xp.get(0);
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  1 + 2*xp.get(1)*xp.get(1);
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  2 + 2*xp.get(1)*xp.get(1);
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }


            ++it2;
        }


        auto eq1 = Lap(v);

        Solver.impose(eq1, bulk, 6);
        Solver.impose(v, up_p, prop_id<1>());
        Solver.impose(v, dw_p, prop_id<1>());
        Solver.impose(v, l_p, prop_id<1>());
        Solver.impose(v, r_p, prop_id<1>());
        Solver.solve(v);
        anasol=Lap(v);

        double worst1 = 0.0;

        it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(domain.getProp<0>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<0>(p) - domain.getProp<2>(p));
            }

            domain.getProp<1>(p) = fabs(domain.getProp<0>(p) - domain.getProp<2>(p));

            ++it2;
        }

        domain.write("particles");
        BOOST_REQUIRE(worst1 < 0.03);
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_poisson_Robin_anal) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {81,81};
        Box<2, double> box({0, 0}, {0.5, 0.5});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3.1);
        double rCut = 3.1 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,double>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            ++it;
        }

        // Add multi res patch 1

        {
        const size_t sz2[2] = {40,40};
        Box<2,double> bx({0.25 + it.getSpacing(0)/4.0,0.25 + it.getSpacing(0)/4.0},{sz2[0]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0, sz2[1]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<2,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[1] = y;

            ++it2;
        }
        }

        // Add multi res patch 2

        {
        const size_t sz2[2] = {40,40};
        Box<2,double> bx({0.25 + 21.0*spacing/8.0,0.25 + 21.0*spacing/8.0},{sz2[0]*spacing/4.0 + 0.25 + 21.0*spacing/8.0, sz2[1]*spacing/4.0 + 0.25 + 21*spacing/8.0});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<2,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[1] = y;

            ++it2;
        }
        }

        ///////////////////////

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut / 3.0 ,1.9/*,support_options::RADIUS*/);
        Derivative_y Dy(domain, 2, rCut / 3.0 ,1.9/*,support_options::RADIUS*/);
        Laplacian Lap(domain, 2, rCut / 3.0 ,1.9/*,support_options::RADIUS*/);


        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            }
            ++it2;
        }

        domain.ghost_get<1,3>();

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);
        auto Poisson = Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(D_y, up_p, 0);
        Solver.impose(D_x, r_p, 0);
        Solver.impose(v, dw_p, 0);
        Solver.impose(v, l_p, 0);

        petsc_solver<double> solver;



        solver.setRestart(500);

        Solver.solve_with_solver(solver,sol);

        domain.ghost_get<2>();

        DCPSE_sol=Lap(sol);
        domain.ghost_get<5>();

        double worst1 = 0.0;

        v=abs(DCPSE_sol-RHS);

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));

        }
        std::cout << "Maximum Analytic Error: " << worst1 << std::endl;

        domain.ghost_get<4>();
        domain.write("Robin_anasol");
        BOOST_REQUIRE(worst1 < 0.03);

    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //81 0.00131586
    //161 0.000328664
    //320 8.30297e-05
    //520 3.12398e-05
    //1024 8.08087e-06
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_poisson_Dirichlet_anal) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {81,81};
        Box<2, double> box({0, 0}, {1, 1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3.1);
        double rCut = 3.1 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,double>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut,1.9,support_options::RADIUS);
        Derivative_y Dy(domain, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(domain, 2, rCut, 1.9,support_options::RADIUS);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkF;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                bulkF.add();
                bulkF.last().get<0>() = p.getKey();
            }
            ++it2;
        }
        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);
        auto Poisson = Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(v, up_p, prop_id<1>());
        Solver.impose(v, r_p, prop_id<1>());
        Solver.impose(v, dw_p, prop_id<1>());
        Solver.impose(v, l_p, prop_id<1>());
        Solver.solve(sol);
        DCPSE_sol=Lap(sol);


        for (int j = 0; j < up_p.size(); j++) {
            auto p = up_p.get<0>(j);
            domain.getProp<5>(p) = 0;
        }
        for (int j = 0; j < dw_p.size(); j++) {
            auto p = dw_p.get<0>(j);
            domain.getProp<5>(p) = 0;
        }
        for (int j = 0; j < l_p.size(); j++) {
            auto p = l_p.get<0>(j);
            domain.getProp<5>(p) = 0;
        }
        for (int j = 0; j < r_p.size(); j++) {
            auto p = r_p.get<0>(j);
            domain.getProp<5>(p) = 0;
        }

        double worst1 = 0.0;

        v=abs(DCPSE_sol-RHS);

        for(int j=0;j<bulkF.size();j++)
        {   auto p=bulkF.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));

        }
        BOOST_REQUIRE(worst1 < 0.03);

        domain.write("Dirichlet_anasol");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_poisson_Periodic) {
        //https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/periodic/python/documentation.html
        //  int rank;
        //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3.1);
        double rCut = 3.1 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> domain(0, box, bc, ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1)*0.99999;
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();


        Laplacian Lap(domain, 2, rCut, 1.9, support_options::RADIUS);

        DCPSE_scheme<equations2d1p,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto u = getV<5>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            //domain.getProp<3>(p)=1+xp[0]*xp[0]+2*xp[1]*xp[1];
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  0;
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  0;
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = xp.get(0)*sin(5*M_PI*xp.get(1))+exp(-((xp.get(0)-0.5)*(xp.get(0)-0.5)+(xp.get(1)-0.5)*(xp.get(1)-0.5))/0.02);
            }

            ++it2;
        }

        domain.ghost_get<1>();
        auto Poisson = -Lap(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(v, up_p, 0);
        Solver.impose(v, dw_p, 0);
        Solver.solve(v);

        domain.ghost_get<0>();
        anasol=-Lap(v);
        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<1>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<1>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<3>(p));

        }
        //Auto Error
        BOOST_REQUIRE(worst1 < 1.0);

        domain.write("Poisson_Periodic");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_poisson_Robin) {
        //http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation#Full_Neumann_boundary_conditions
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_y Dy(domain, 2, rCut,2);
        Laplacian Lap(domain, 2, rCut, 3);

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto u = getV<5>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            //domain.getProp<3>(p)=1+xp[0]*xp[0]+2*xp[1]*xp[1];
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -10.0*exp(-((xp.get(0)-0.5)*(xp.get(0)-0.5)+(xp.get(1)-0.5)*(xp.get(1)-0.5))/0.02);
            }

            ++it2;
        }


        auto Poisson = Lap(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(D_y, up_p, prop_id<1>());
        Solver.impose(-D_y, dw_p, prop_id<1>());
        Solver.impose(v, l_p, 0);
        Solver.impose(v, r_p, 0);
        Solver.solve(sol);
        anasol=Lap(sol);
        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<1>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<1>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<3>(p));

        }
        //Auto Error
        BOOST_REQUIRE(worst1 < 1.0);

        domain.write("Mixed");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_poisson_Neumann) {
    //https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/neumann-poisson/python/documentation.html
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(spacing * 3.1);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {

            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut,1.9);
        Derivative_y Dy(domain, 2, rCut,1.9);
        Laplacian Lap(domain, 2, rCut, 1.9);
        petsc_solver<double> solver;
        solver.setRestart(500);
        solver.setSolver(KSPGMRES);
        solver.setPreconditioner(PCSVD);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        //auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            //domain.getProp<3>(p)=1+xp[0]*xp[0]+2*xp[1]*xp[1];
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -10*exp(-((xp.get(0)-0.5)*(xp.get(0)-0.5)+(xp.get(1)-0.5)*(xp.get(1)-0.5))/0.02);
            }

            ++it2;
        }

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain,options_solver::LAGRANGE_MULTIPLIER);
        auto Poisson = -Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(D_y, up_p, prop_id<1>());
        Solver.impose(-D_y, dw_p, prop_id<1>());
        Solver.impose(-D_x, l_p, prop_id<1>());
        Solver.impose(D_x, r_p, prop_id<1>());
        Solver.solve_with_solver(solver,sol);

//       Solver.solve(sol);
        anasol=-Lap(sol);
        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<1>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<1>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<3>(p));

        }
        //Auto Error
        BOOST_REQUIRE(worst1 < 1.0);

        domain.write("Neumann");
    }


    BOOST_AUTO_TEST_CASE(dcpse_slice_solver) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1, 1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut,1.9,support_options::RADIUS);
        Derivative_y Dy(domain, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(domain, 2, rCut,1.9,support_options::RADIUS);


        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));


            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            }
            ++it2;
        }

        eq_id vx,vy;

        vx.setId(0);
        vy.setId(1);

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);
        auto Poisson0 = Lap(v[0]);
        auto Poisson1 = Lap(v[1]);
        //auto D_x = Dx(v[1]);
        //auto D_y = Dy(v[1]);
        Solver.impose(Poisson0, bulk, RHS[0],vx);
        Solver.impose(Poisson1, bulk, RHS[1],vy);
        Solver.impose(v[0], up_p, RHS[0],vx);
        Solver.impose(v[1], up_p, RHS[1],vy);
        Solver.impose(v[0], r_p,  RHS[0],vx);
        Solver.impose(v[1], r_p,  RHS[1],vy);
        Solver.impose(v[0], dw_p, RHS[0],vx);
        Solver.impose(v[1], dw_p, RHS[1],vy);
        Solver.impose(v[0], l_p,  RHS[0],vx);
        Solver.impose(v[1], l_p,  RHS[1],vy);
        Solver.solve(sol[0],sol[1]);
        DCPSE_sol=Lap(sol);
        double worst1 = 0.0;
        double worst2 = 0.0;


        v=sol-RHS;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p)[0] - domain.getProp<2>(p)[0]) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p)[0] - domain.getProp<2>(p)[0]);
            }
            domain.getProp<4>(p)[0] = fabs(domain.getProp<3>(p)[0] - domain.getProp<2>(p)[0]);

        }
        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p)[1] - domain.getProp<2>(p)[1]) >= worst2) {
                worst2 = fabs(domain.getProp<3>(p)[1] - domain.getProp<2>(p)[1]);
            }
            domain.getProp<4>(p)[1] = fabs(domain.getProp<3>(p)[1] - domain.getProp<2>(p)[1]);

        }
        //std::cout << "Maximum Analytic Error in slice x: " << worst1 << std::endl;
        //std::cout << "Maximum Analytic Error in slice y: " << worst2 << std::endl;
        domain.write("Slice_anasol");
        BOOST_REQUIRE(worst1 < 0.03);
        BOOST_REQUIRE(worst2 < 0.03);

    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_Lid) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        //                                  P        V                 Dv              RHS    Vtemp
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,double>> Particles(0, box, bc, ghost);
        vector_dist<2, double, aggregate<double,VectorS<2, double>,double> > Particles_subset(0, box, bc, ghost);

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
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> ref_bulk;
        openfpm::vector<aggregate<int>> ref_bulk2;
        openfpm::vector<aggregate<int>> bulk_1;


        auto P = getV<0>(Particles);
        auto P_bulk = getV<2>(Particles_subset);
        auto V = getV<1>(Particles);
        auto V_bulk = getV<1>(Particles_subset);
        auto divV_bulk = getV<0>(Particles_subset);
        auto dV = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Vtemp = getV<4>(Particles);
        auto k1 = getV<5>(Particles);
        auto k2 = getV<6>(Particles);
        auto k3 = getV<7>(Particles);
        auto k4 = getV<8>(Particles);
        auto divV = getV<9>(Particles);


        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("Re1000-1e-3-vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                if (xp[0]==0 && xp[1]==0) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                else{
                    dw_p.add();
                    dw_p.last().get<0>() = p.getKey();
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) = 0.0;
            }

            ++it2;
        }


        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
            Particles_subset.getLastProp<1>() = Particles.template getProp<1>(bulk.template get<0>(i));
        }



        Particles_subset.map();
        Particles_subset.ghost_get<0>();

        auto it3 = Particles_subset.getDomainIterator();

        while (it3.isNext()) {
            auto p = it3.get();
            Point<2, double> xp = Particles.getPos(p);

            if (xp[0]==0.5 && xp[1]==0.5) {
                ref_bulk.add();
                ref_bulk.last().get<0>() = p.getKey();
            }
            else if(xp[0]==0.5+spacing && xp[1]==0.5+spacing)
            {
                ref_bulk2.add();
                ref_bulk2.last().get<0>() = p.getKey();
            }
            else
            {
                bulk_1.add();
                bulk_1.last().get<0>() = p.getKey();
            }



            ++it3;
        }

        Derivative_x Dx(Particles, 2, rCut,2);
        Derivative_y Dy(Particles, 2, rCut,2);
        Gradient Grad_sub(Particles_subset, 2, rCut, 2),Grad(Particles, 2, rCut,2);
        Laplacian Lap_sub(Particles_subset, 2, rCut, 2),Lap(Particles, 2, rCut, 2);
        Advection Adv(Particles, 2, rCut, 2);
        Divergence Div(Particles, 2, rCut, 2),Div_bulk(Particles_subset, 2, rCut, 2);


        double dt=0.005;
        int n=50;
        double nu=1e-3;
        divV_bulk=Div_bulk(V_bulk);
        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<3>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            Particles.template getProp<9>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
        }
        //subset_create<0,1,2,4>(Particles,Particles_subset,bulk);
        DCPSE_scheme<equations2d1,decltype(Particles_subset)> Solver(Particles_subset);
        auto Pressure_Poisson = Lap_sub(P_bulk);
        Solver.impose(Pressure_Poisson, bulk_1,prop_id<0>());
        Solver.impose(P_bulk, ref_bulk, 1);
        Solver.impose(P_bulk, ref_bulk2, 1);
        Solver.solve(P_bulk);

        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<0>(bulk.template get<0>(i)) = Particles_subset.getProp<2>(i);
        }

        std::cout<<"Poisson Solved"<<std::endl;
        V_bulk=V_bulk-Grad_sub(P_bulk);
        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<1>(bulk.template get<0>(i)) = Particles_subset.getProp<1>(i);
        }

        divV_bulk=Div_bulk(V_bulk);
        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<9>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
        }

        Particles.write_frame("Re1000-1e-4-Lid",0);
        std::cout<<"Init Done"<<std::endl;
        return;
        for(int i =1; i<=n ;i++)
        {   dV=V+dt*(nu*Lap(V)-Adv(V,V));
            std::cout<<"dV Done"<<std::endl;
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_subset.getProp<1>(i) = Particles.template getProp<2>(bulk.template get<0>(i));
            }
            divV_bulk=1.0/dt*Div_bulk(V_bulk);
            for (int i = 0 ; i < bulk.size() ; i++)
            {
                Particles.template getProp<3>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            }
            std::cout<<"RHS Done"<<std::endl;
            DCPSE_scheme<equations2d1,decltype(Particles)> Solver(Particles);
            auto Pressure_Poisson = Lap(P);
            auto D_x = Dx(P);
            auto D_y = Dy(P);
            Solver.impose(Pressure_Poisson, bulk, prop_id<3>());
            Solver.impose(D_y, up_p, 0);
            Solver.impose(-D_y, dw_p,0);
            Solver.impose(-D_x, l_p,0);
            Solver.impose(D_x, r_p, 0);
            Solver.impose(P, ref_p, 0);
            Solver.solve(P);
            std::cout<<"Poisson Solved"<<std::endl;
            Vtemp=dV-dt*Grad(P);
            V=Vtemp;
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            Particles.getProp<1>(0)[0] =  0;
            Particles.getProp<1>(0)[1] =  0;
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_subset.getProp<1>(i) = Particles.template getProp<1>(bulk.template get<0>(i));
            }
            divV_bulk=Div_bulk(V_bulk);
            for (int i = 0 ; i < bulk.size() ; i++)
            {
                Particles.template getProp<9>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            }
            //divV=Div(V);
            std::cout<<"Velocity updated"<<std::endl;
            Particles.write_frame("Re1000-1e-4-Lid",i);
            std::cout<<i<<std::endl;
        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_Lid_sf) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {10,10});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut =spacing*(3.1);
        double ord = 3;
        double sampling = 2.3;
        double rCut2 = 3.1*spacing;
        double ord2 = 2;
        double sampling2 = 1.9;

        Ghost<2, double> ghost(rCut);
        std::cout<<spacing<<std::endl;
        //                                  sf    W   DW        RHS    Wnew  V
        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> Particles(0, box, bc, ghost);

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

        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;

        openfpm::vector<aggregate<int>> BULK;



        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) -  spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) +  spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                            {box.getHigh(0) -  spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> down({box.getLow(0) +  spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) -  spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) +   spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) -  spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) -  spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) +  spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) + spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                              {box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) -  spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) +  spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                               {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});
        Box<2, double> Bulk_box({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                                {box.getHigh(0) + spacing / 2.0, box.getHigh(1)  +spacing / 2.0});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(up_d);
        boxes.add(down);
        boxes.add(down_u);
        boxes.add(left);
        boxes.add(left_r);
        boxes.add(right);
        boxes.add(right_l);
        boxes.add(Bulk_box);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("box_sf.vtk");
        // Particles.write_frame("Re1000-1e-3-Lid_sf",0);

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;//-M_PI*M_PI*sin(M_PI*xp[0])*sin(M_PI*xp[1]);
            Particles.getProp<1>(p) =0;//sin(M_PI*xp[0])*sin(M_PI*xp[1]);
            Particles.getProp<2>(p) =0.0;
            Particles.getProp<3>(p) =0.0;
            Particles.getProp<4>(p) =0.0;
            Particles.getProp<5>(p)[0] =0.0;
            Particles.getProp<5>(p)[1] =0.0;

            BULK.add();
            BULK.last().get<0>() = p.getKey();

            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            }
            else if (Bulk_box.isInside(xp) == true)
            {
                if (up_d.isInside(xp) == true) {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                }
                else if (down_u.isInside(xp) == true) {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                }else if (left_r.isInside(xp) == true) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                } else if (right_l.isInside(xp) == true) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }

        for(int j=0;j<up_p.size();j++) {
            auto p = up_p.get<0>(j);
            Particles.getProp<1>(p) =  -12.0/spacing;
        }

        Laplacian Lap(Particles, ord2, rCut2,sampling,support_options::RADIUS);
        //Gradient Grad(Particles, ord, rCut,sampling,support_options::RADIUS);
        Derivative_xy Dxy(Particles,ord, rCut, sampling,support_options::RADIUS);
        Derivative_xy Dxy2(Particles,ord2, rCut2, sampling2,support_options::RADIUS);
        Derivative_x Dx(Particles,ord, rCut, sampling,support_options::RADIUS);
        Derivative_y Dy(Particles,ord, rCut, sampling,support_options::RADIUS);


        auto its = Particles.getDomainIterator();
        int ctr=0;
        while (its.isNext()) {
            auto p = its.get();
            Dx.DrawKernel<0>(Particles, p.getKey());
            Dy.DrawKernel<1>(Particles, p.getKey());
            Dxy.DrawKernel<2>(Particles, p.getKey());
            Dxy2.DrawKernel<3>(Particles, p.getKey());
            Lap.DrawKernel<4>(Particles, p.getKey());
            //Grad.DrawKernel<4>(Particles, p.getKey());

            Particles.write_frame("LapKer",ctr);
            for(int j=0;j<BULK.size();j++) {
                auto p1 = BULK.get<0>(j);
                Particles.getProp<0>(p1) =  0;
                Particles.getProp<1>(p1) =  0;
                Particles.getProp<2>(p1) =  0;
                Particles.getProp<3>(p1) =  0;
                Particles.getProp<4>(p1) =  0;
                Particles.getProp<5>(p1)[0] =  0;
                Particles.getProp<5>(p1)[1] =  0;

            }
            ++its;
            ctr++;
        }


    }

BOOST_AUTO_TEST_SUITE_END()


