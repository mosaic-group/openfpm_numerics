//
// Created by Abhinav Singh on 01.05.20.
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


BOOST_AUTO_TEST_SUITE(fd_op_suite_tests)
    BOOST_AUTO_TEST_CASE(dcpse_op_tests) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);

        BOOST_TEST_MESSAGE("Init vector_dist...");
        vector_dist<2, double, aggregate<double, double, double>> domain(0, box, bc,ghost);

        BOOST_TEST_MESSAGE("Init domain...");
        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        while (it.isNext()) {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;
            // Here fill the function value P
            domain.template getLastProp<0>() = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            // Here fill the validation value for Df/Dx in property 3
            domain.template getLastProp<3>() = cos(domain.getLastPos()[0]);
            ++counter;
            ++it;
        }


        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, Equations2d1, CENTRAL);

        auto v = getV<1>(domain);
        auto P = getV<0>(domain);

        v = Dx(P);
        auto it2 = domain.getDomainIterator();

        double worst = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();

            if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) > worst) {
                worst = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
            }

            ++it2;
        }

        std::cout << "Maximum Error: " << worst << std::endl;

        domain.deleteGhost();
        domain.write("v");
    }



















BOOST_AUTO_TEST_SUITE_END()

