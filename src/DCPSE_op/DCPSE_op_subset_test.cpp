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
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 2.0 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, double, double, VectorS<2, double>, VectorS<2, double>>> domain(0, box,
                                                                                                                 bc,
                                                                                                                 ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");
        std::mt19937 rng{6666666};

        std::normal_distribution<> gaussian{0, sigma2};

        openfpm::vector<aggregate<int>> no_border_subset;

        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext())
        {
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
            domain.template getLastProp<4>()[0] = -sin(domain.getLastPos()[0]);
            domain.template getLastProp<4>()[1] = -sin(domain.getLastPos()[1]);

            if (k0 != 0 && k1 != 0 && k0 != sz[0] -1 && k1 != sz[1] - 1)
            {
            	no_border_subset.add();
            	no_border_subset.template get<0>(no_border_subset.size()-1) = domain.size_local() - 1;
            }

            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        vector_dist_subset<2, double, aggregate<double, double, double, VectorS<2, double>, VectorS<2, double>>> domain_bulk(domain,no_border_subset);

        // move particles

        auto P = getV<0>(domain);
        auto P_bulk = getV<2>(domain_bulk);
        auto Out_bulk = getV<1>(domain_bulk);
	auto Out_V_bulk = getV<3>(domain_bulk);

        P_bulk = 5;

//        domain.write("Test_output_subset");

        // Create the subset

        Derivative_x Dx(domain, 2, rCut);
        Derivative_y Dy(domain, 2, rCut);
        Derivative_x Dx_bulk(domain_bulk, 2, rCut);

        Out_bulk = Dx_bulk(P);
	Out_V_bulk[0] = P + Dx_bulk(P);

//        P_bulk = Dx_bulk(P_bulk);  <------------ Incorrect produce error message
//        P = Dx_bulk(P);   <------- Incorrect produce overflow

        domain.write("Out");
    }


BOOST_AUTO_TEST_SUITE_END()
