/*
 * DCPSE_op_test.cpp
 *
 *  Created on: April 9, 2020
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



BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests3)
    BOOST_AUTO_TEST_CASE(dcpse_op_vec3d) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 21;
        const size_t sz[3] = {edgeSemiSize,  edgeSemiSize,edgeSemiSize};
        Box<3, double> box({0, 0,0}, {1,1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<3, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing * spacing/ (2 * 4);

        vector_dist<3, double, aggregate<double, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>,double,double>> domain(
                0, box, bc, ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext()) {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing;
            domain.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing;
            domain.getLastPos()[1] = y;//+gaussian(rng);
            mem_id k2 = key.get(2);
            double z = k2 * spacing;
            domain.getLastPos()[2] = z;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<0>()    = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]) + sin(domain.getLastPos()[2]) ;
            domain.template getLastProp<1>()[0] = cos(domain.getLastPos()[0]);
            domain.template getLastProp<1>()[1] = cos(domain.getLastPos()[1]) ;
            domain.template getLastProp<1>()[2] = cos(domain.getLastPos()[2]);
            // Here fill the validation value for Df/Dx
            domain.template getLastProp<2>()[0] = 0;//cos(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<2>()[1] = 0;//-sin(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<3>()[0] = 0;//cos(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<3>()[1] = 0;//-sin(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<3>()[2] = 0;

            domain.template getLastProp<4>()[0] = -cos(domain.getLastPos()[0]) * sin(domain.getLastPos()[0]);
            domain.template getLastProp<4>()[1] = -cos(domain.getLastPos()[1]) * sin(domain.getLastPos()[1]);
            domain.template getLastProp<4>()[2] = -cos(domain.getLastPos()[2]) * sin(domain.getLastPos()[2]);


            /*  domain.template getLastProp<4>()[0] = cos(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) +
                                                    cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
              domain.template getLastProp<4>()[1] = -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                                                    sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
              domain.template getLastProp<4>()[2] = -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                                                    sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));*/
            domain.template getLastProp<5>()    = cos(domain.getLastPos()[0]) * cos(domain.getLastPos()[0])+cos(domain.getLastPos()[1]) * cos(domain.getLastPos()[1])+cos(domain.getLastPos()[2]) * cos(domain.getLastPos()[2]) ;
            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Advection Adv(domain, 2, rCut, 1.9,support_options::RADIUS);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);
        auto dv = getV<3>(domain);
        auto dP = getV<6>(domain);


//        typedef boost::mpl::int_<std::is_fundamental<point_expression_op<Point<2U, double>, point_expression<double>, Point<2U, double>, 3>>::value>::blabla blabla;
//        std::is_fundamental<decltype(o1.value(key))>

        domain.ghost_get<1>();
        dv = Adv(v, v);
        auto it2 = domain.getDomainIterator();

        double worst1 = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();

            if (fabs(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1]) > worst1) {
                worst1 = fabs(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1]);

            }

            ++it2;
        }
        //std::cout << "Maximum Error in component 2: " << worst1 << std::endl;
        BOOST_REQUIRE(worst1 < 0.03);

        //Adv.checkMomenta(domain);
        //Adv.DrawKernel<2>(domain,0);

        //domain.deleteGhost();

        dP = Adv(v, P);//+Dy(P);
        auto it3 = domain.getDomainIterator();

        double worst2 = 0.0;

        while (it3.isNext()) {
            auto p = it3.get();
            if (fabs(domain.getProp<6>(p) - domain.getProp<5>(p)) > worst2) {
                worst2 = fabs(domain.getProp<6>(p) - domain.getProp<5>(p));

            }

            ++it3;
        }
        domain.deleteGhost();
        BOOST_REQUIRE(worst2 < 0.03);


    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_poisson_dirichlet_anal3d) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t grd_sz=21;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({0, 0,0}, {1.0, 1.0,1.0});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<3, double> ghost(spacing * 3.1);
        double rCut = 3.1 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<3, double, aggregate<double,double,double,double,double,double>> domain(0, box, bc, ghost);


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
            double z = key.get(2) * it.getSpacing(2);
            domain.getLastPos()[2] = z;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut,3.1,support_options::RADIUS);
        Derivative_y Dy(domain, 2, rCut,3.1,support_options::RADIUS);
        Laplacian Lap(domain, 2, rCut,1.9,support_options::RADIUS);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> front_p;
        openfpm::vector<aggregate<int>> back_p;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> down_p;
        openfpm::vector<aggregate<int>> left_p;
        openfpm::vector<aggregate<int>> right_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        Box<3, double> up(
                {box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> down(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> left(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> right(
                {box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> front(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getLow(2) + spacing / 2.0});

        Box<3, double> back(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getHigh(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(front);
        boxes.add(back);
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        //vtk_box.write("boxes_3d.vtk");
        auto Particles=domain;
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            domain.getProp<1>(p) = -3.0*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1))*sin(M_PI*xp.get(2));
            domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1))*sin(M_PI*xp.get(2));
            if (front.isInside(xp) == true) {
                front_p.add();
                front_p.last().get<0>() = p.getKey();
            } else if (back.isInside(xp) == true) {
                back_p.add();
                back_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                left_p.add();
                left_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                right_p.add();
                right_p.last().get<0>() = p.getKey();
            } else if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                down_p.add();
                down_p.last().get<0>() = p.getKey();
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }


        DCPSE_scheme<equations3d1,decltype(domain)> Solver( domain);
        auto Poisson = Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(v, up_p, prop_id<1>());
        Solver.impose(v, right_p, prop_id<1>());
        Solver.impose(v, down_p, prop_id<1>());
        Solver.impose(v, left_p, prop_id<1>());
        Solver.impose(v, front_p, prop_id<1>());
        Solver.impose(v, back_p, prop_id<1>());
        Solver.solve(sol);
        DCPSE_sol=Lap(sol);

        double worst1 = 0.0;

        v=abs(DCPSE_sol-RHS);

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));

        }
        //std::cout << "Maximum Analytic Error: " << worst1 << std::endl;

        BOOST_REQUIRE(worst1 < 0.03);

        //domain.write("Dirichlet_anasol_3d");




    }

BOOST_AUTO_TEST_SUITE_END()


