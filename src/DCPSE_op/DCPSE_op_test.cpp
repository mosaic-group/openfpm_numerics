/*
 * DCPSE_op_test.cpp
 *
 *  Created on: Jan 7, 2020
 *      Author: Abhinav Singh, Pietro Incardona
 *
 */

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "DCPSE_op.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE( dcpse_op_suite_tests )

BOOST_AUTO_TEST_CASE( dcpse_op_tests )
{
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size_t edgeSemiSize = 160;
    const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
    Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
    size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
    double spacing[2];
    spacing[0] = 2*M_PI / (sz[0] - 1);
    spacing[1] = 2*M_PI / (sz[1] - 1);
    Ghost<2, double> ghost(spacing[0]*3);
    double rCut = 2.0*spacing[0];
    BOOST_TEST_MESSAGE("Init vector_dist...");
    double sigma2 = spacing[0] * spacing[1] / (2 * 4);

    vector_dist<2, double, aggregate<double, double, double>> domain(0, box, bc, ghost);

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
    while (it.isNext())
    {
        domain.add();
        auto key = it.get();
        mem_id k0 = key.get(0);
        double x = k0 * spacing[0];
        domain.getLastPos()[0] = x + gaussian(rng);
        mem_id k1 = key.get(1);
        double y = k1 * spacing[1];
        domain.getLastPos()[1] = y + gaussian(rng);
        // Here fill the function value
        domain.template getLastProp<0>() = sin(domain.getLastPos()[0]);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
        // Here fill the validation value for Df/Dx
        domain.template getLastProp<2>() = cos(domain.getLastPos()[0]);
//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

        ++counter;
        ++it;
    }
    BOOST_TEST_MESSAGE("Sync domain across processors...");

    domain.map();
    domain.ghost_get<0>();

    Derivative_x Dx(domain,2,rCut);
	//Dy Dy;
	auto v = getV<1>(domain);
	auto P = getV<0>(domain);


	v = Dx(P);//+ Dy(F);

	auto it2 = domain.getDomainIterator();

	double worst = 0.0;

	while (it2.isNext())
    {
	    auto p = it2.get();

	    if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) > worst)
        {
	        worst = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
        }

	    ++it2;
    }

	std::cout << "WORST " << worst << std::endl;

	domain.deleteGhost();
	domain.write("Dx");

	std::cout << demangle(typeid(decltype(v)).name()) << "\n";

	//Debug<decltype(expr)> a;

	//typedef decltype(expr)::blabla blabla;

	//auto err = Dx + Dx;
}

BOOST_AUTO_TEST_SUITE_END()


