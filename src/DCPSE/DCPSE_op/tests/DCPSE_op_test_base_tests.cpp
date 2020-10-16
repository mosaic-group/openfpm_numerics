/*
 * DCPSE_op_test_base_tests.cpp
 *
 *  Created on: April 9, 2020
 *      Author: Abhinav Singh
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

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)
BOOST_AUTO_TEST_CASE(dcpse_op_tests) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 2.0 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, double, double, VectorS<2, double>, VectorS<2, double>>> domain(0, box,
                                                                                                                 bc,
                                                                                                                 ghost);

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
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<0>() = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            domain.template getLastProp<2>() = cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]);
            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut);
        Derivative_y Dy(domain, 2, rCut);
        Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);

        v = Dx(P) + Dy(P);
        auto it2 = domain.getDomainIterator();

        double worst = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();

            if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) > worst) {
                worst = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
            }

            ++it2;
        }

        domain.deleteGhost();
        BOOST_REQUIRE(worst < 0.03);

    }


    BOOST_AUTO_TEST_CASE(dcpse_op_test_lap) {
        size_t edgeSemiSize = 81;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 2.0 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, double, double, double, VectorS<2, double>>> domain(0, box,
                                                                                                     bc,
                                                                                                     ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");
        std::normal_distribution<> gaussian{0, sigma2};

        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext()) {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<0>() = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            // Here fill the validation value for Lap
            domain.template getLastProp<1>() = - sin(domain.getLastPos()[0]) - sin(domain.getLastPos()[1]);
            domain.template getLastProp<4>()[0] = -sin(domain.getLastPos()[0]);
            domain.template getLastProp<4>()[1] = -sin(domain.getLastPos()[1]);

            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Laplacian Lap(domain, 2, rCut,1.9);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);
        auto vv = getV<2>(domain);
        auto errv = getV<3>(domain);

        vv = Lap(P);
        auto it2 = domain.getDomainIterator();

        double worst = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();

            if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) > worst) {
                worst = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
            }

            ++it2;
        }
        errv=v-vv;
        domain.deleteGhost();
        BOOST_REQUIRE(worst < 0.3);
    }

    BOOST_AUTO_TEST_CASE(dcpse_op_div) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 31;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {10.0, 10.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = box.getHigh(0) / (sz[0] - 1);
        spacing[1] = box.getHigh(1) / (sz[1] - 1);
        double rCut = 3.1* spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>,double,double>> domain(
                0, box, bc, ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");
//            std::random_device rd{};
//            std::mt19937 rng{rd()};
        std::mt19937 rng{6666666};

        std::normal_distribution<> gaussian{0, sigma2};

        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext()) {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<1>()[0] = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            domain.template getLastProp<1>()[1] = cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]);


            // Here fill the validation value for Divergence
            domain.template getLastProp<0>()= cos(domain.getLastPos()[0]) - sin(domain.getLastPos()[1]);
           /* domain.template getLastProp<4>()[0] =
                    cos(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) +
                    cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
            domain.template getLastProp<4>()[1] =
                    -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                    sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));*/


            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Divergence Div(domain, 2, rCut, 1.9,support_options::RADIUS);
        Derivative_x Dx(domain, 2, rCut, 1.9,support_options::RADIUS);
        Derivative_y Dy(domain, 2, rCut, 1.9,support_options::RADIUS);

        auto v = getV<1>(domain);
        auto anasol = getV<0>(domain);
        auto div = getV<5>(domain);
        auto div2 = getV<6>(domain);

        div = Div(v);
        div2=Dx(v[0])+Dy(v[1]);
        auto it2 = domain.getDomainIterator();

        double worst1 = 0.0;
        double worst2 = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(domain.getProp<0>(p) - domain.getProp<5>(p)) > worst1) {
                worst1 = fabs(domain.getProp<0>(p) - domain.getProp<5>(p));
            }
            if (fabs(domain.getProp<0>(p) - domain.getProp<6>(p)) > worst2) {
                worst2 = fabs(domain.getProp<0>(p) - domain.getProp<6>(p));
            }
            ++it2;
        }

        domain.deleteGhost();
/*
        domain.write("DIV");

        std::cout<<worst1<<":"<<worst2;
*/

        BOOST_REQUIRE(worst1 < 0.05);
        BOOST_REQUIRE(worst2 < 0.05);

    }

    BOOST_AUTO_TEST_CASE(dcpse_op_vec) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 160;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {10,10});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = box.getHigh(0)  / (sz[0] - 1);
        spacing[1] = box.getHigh(1) / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3.1);
        double rCut = 3.1 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>,double,double>> domain(
                0, box, bc, ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");
//            std::random_device rd{};
//            std::mt19937 rng{rd()};
        std::mt19937 rng{6666666};

        std::normal_distribution<> gaussian{0, sigma2};

        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext()) {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<0>()    = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            domain.template getLastProp<1>()[0] = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            domain.template getLastProp<1>()[1] = cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]);

            // Here fill the validation value for Df/Dx
            domain.template getLastProp<2>()[0] = 0;
            domain.template getLastProp<2>()[1] = 0;
            domain.template getLastProp<4>()[0] =
                    cos(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) +
                    cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
            domain.template getLastProp<4>()[1] =
                    -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                    sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));

            domain.template getLastProp<5>()    = cos(domain.getLastPos()[0])*(sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]))+cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));


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
        //domain.write("v1");
        BOOST_REQUIRE(worst1 < 0.03);


        dP = Adv(v, P);
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
        //domain.write("v2");
        BOOST_REQUIRE(worst2 < 0.03);

    }


    BOOST_AUTO_TEST_CASE(dcpse_slice) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;

        vector_dist<2, double, aggregate<double,VectorS<2, double>,double,double[3][3]>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;

            Particles.getLastProp<1>()[0] = sin(x+y);
            Particles.getLastProp<1>()[1] = cos(x+y);

            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto S = getV<2>(Particles);
        auto Sig = getV<3>(Particles);


        Derivative_x Dx(Particles, 2, rCut,2);

        P = Dx(V[0]);
        S = V[0]*V[0] + V[1]*V[1];

        Sig[0][1] = V[0]*V[0] + V[1]*V[1];
        Sig[1][0] = P;
        Sig[2][2] = 5.0;

        auto it2 = Particles.getDomainIterator();

        double err = 0.0;

        while (it2.isNext())
        {
            auto p = it2.get();

            if (fabs(Particles.getProp<0>(p) - Particles.getProp<1>(p)[1]) >= err )
            {
                err = fabs(Particles.getProp<0>(p) - Particles.getProp<1>(p)[1]);
            }

            if (fabs(Particles.getProp<2>(p) - 1.0) >= err )
            {
                err = fabs(Particles.getProp<2>(p) - 1.0);
            }

            if (fabs(Particles.getProp<3>(p)[0][1] - 1.0) >= err )
            {
                err = fabs(Particles.getProp<3>(p)[0][1] - 1.0);
            }

            if (fabs(Particles.getProp<3>(p)[1][0] - Particles.getProp<1>(p)[1]) >= err )
            {
                err = fabs(Particles.getProp<3>(p)[0][1] - Particles.getProp<1>(p)[1]);
            }

            if (fabs(Particles.getProp<3>(p)[2][2] - 5.0) >= err )
            {
                err = fabs(Particles.getProp<3>(p)[2][2] - 5.0);
            }

            ++it2;
        }

        Particles.write("test_out");
        BOOST_REQUIRE(err < 0.03);

    }

    BOOST_AUTO_TEST_CASE(dcpse_slice_3d) {
        const size_t sz[3] = {17,17,17};
        Box<3, double> box({0, 0,0}, {1,1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<3, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;

        vector_dist<3, double, aggregate<double,VectorS<3, double>,double,double[3][3]>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            double z = key.get(2) * it.getSpacing(2);
            Particles.getLastPos()[2] = z;

            Particles.getLastProp<1>()[0] = sin(x+y);
            Particles.getLastProp<1>()[1] = cos(x+y);
            Particles.getLastProp<1>()[2] = 1.0;

            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto S = getV<2>(Particles);
        auto Sig = getV<3>(Particles);


        Derivative_x Dx(Particles, 2, rCut,2);

        P = Dx(V[0]);
        S = V[0]*V[0] + V[1]*V[1]+V[2]*V[2];

        Sig[0][1] = V[0]*V[0] + V[1]*V[1]+V[2]*V[2];
        Sig[1][0] = P;
        Sig[2][2] = 5.0;

        auto it2 = Particles.getDomainIterator();

        double err1 = 0.0;
        double err2 = 0.0;
        double err3 = 0.0;
        double err4 = 0.0;
        double err5 = 0.0;

        while (it2.isNext())
        {
            auto p = it2.get();

            // Here we check that P = Dx(V[0]) works   P is Prop 0 = sin(x+y) V[0] = sin(x+y) and V[1] is cos(x+y)
            if (fabs(Particles.getProp<0>(p) - Particles.getProp<1>(p)[1]) >= err1 )
            {
                err1 = fabs(Particles.getProp<0>(p) - Particles.getProp<1>(p)[1]);
            }

            if (fabs(Particles.getProp<2>(p) - 2.0) >= err2 )
            {
                err2 = fabs(Particles.getProp<2>(p) - 2.0);
            }

            if (fabs(Particles.getProp<3>(p)[0][1] - 2.0) >= err3 )
            {
                err3 = fabs(Particles.getProp<3>(p)[0][1] - 2.0);
            }

            // V[0]*V[0] + V[1]*V[1]+V[2]*V[2] = 2 and
            if (fabs(Particles.getProp<3>(p)[0][1] - Particles.getProp<2>(p)) >= err4 )
            {
                err4 = fabs(Particles.getProp<3>(p)[0][1] - Particles.getProp<2>(p));
            }

            if (fabs(Particles.getProp<3>(p)[2][2] - 5.0) >= err5 )
            {
                err5 = fabs(Particles.getProp<3>(p)[2][2] - 5.0);
            }

            ++it2;
        }

        //std::cout << err1 << " " << err2 << " " << err3 << " " << err4 << " " << err5 << std::endl;
        BOOST_REQUIRE(err1 < 0.08);
        BOOST_REQUIRE(err2 < 0.03);
        BOOST_REQUIRE(err3 < 0.03);
        BOOST_REQUIRE(err4 < 0.03);
        BOOST_REQUIRE(err5 < 0.03);
    }


BOOST_AUTO_TEST_SUITE_END()






