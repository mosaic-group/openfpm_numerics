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


BOOST_AUTO_TEST_CASE( dcpse_Navier)
    {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 0.1 / (sz[0] - 1);
        spacing[1] = 0.1 / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0]*3);
        double rCut = 2.0*spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        constexpr int Velocity          =     0;
        constexpr int Pressure          =     1;
        constexpr int bflag             =     2;
        double nu=1/3200.0;
        //                                    v              P    bflag
        vector_dist<2, double, aggregate<VectorS<2,double>,double,int>> domain(0, box, bc, ghost);
        //Init_DCPSE(domain)
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
            // Initial Velocity field and Pressure
            if (x == 0 || y==0 || x==0.1)
            {
                domain.template getLastProp<bflag>() = 1;
                domain.template getLastProp<Velocity>()[0] = 0.0;
                domain.template getLastProp<Velocity>()[1] = 0.0;
                domain.template getLastProp<Pressure>()    = 0.01;
                if(x==0 && y==0)
                {
                    domain.template getLastProp<bflag>() = 3;
                    domain.template getLastProp<Velocity>()[0] = 0.0;
                    domain.template getLastProp<Velocity>()[1] = 0.0;
                    domain.template getLastProp<Pressure>()    = 0.0;
                }
            }
            else{
                domain.template getLastProp<bflag>() = 0;
                domain.template getLastProp<Velocity>()[0] = 0.0;
                domain.template getLastProp<Velocity>()[1] = 0.0;
                domain.template getLastProp<Pressure>()    = 0.01;
            }
            if(y==0.1)
            {
                domain.template getLastProp<bflag>() = 2;
                domain.template getLastProp<Velocity>()[0] = 1.0;
                domain.template getLastProp<Velocity>()[1] = 0.0;
                domain.template getLastProp<Pressure>()    = 0.01;
            }
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain,2,rCut);
        Derivative_y Dy(domain,2,rCut);
        Gradient Grad(domain,2,rCut);
        Laplacian Lap(domain,2,rCut);
        auto v = getV<0>(domain);
        auto P = getV<1>(domain);
        std::cout<<"RHS Initiated"<<std::endl;

        v = -Grad(P);//-Grad(v)+Lap(v);

        std::cout<<"RHS Computed"<<std::endl;
        auto it2 = domain.getDomainIterator();
        while (it2.isNext())
        {   auto p = it2.get();
            if(domain.getProp<bflag>(p)==2)
            {
                domain.getProp<Velocity>(p)[0] = 1.0;
                domain.getProp<Velocity>(p)[1] = 0.0;
            }
            ++it2;
        }
        domain.deleteGhost();
        domain.write("Nav");

        //std::cout << demangle(typeid(decltype(v)).name()) << "\n";

        //Debug<decltype(expr)> a;

        //typedef decltype(expr)::blabla blabla;

        //auto err = Dx + Dx;
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_op_vec)
    {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 40;
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

        vector_dist<2, double, aggregate<double, VectorS<2,double>, VectorS<2,double>, VectorS<2,double>,VectorS<2,double>>> domain(0, box, bc, ghost);

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
            domain.getLastPos()[0] = x ;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y ;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<1>()[0] = sin(domain.getLastPos()[0])+sin(domain.getLastPos()[1]);
            domain.template getLastProp<1>()[1] = cos(domain.getLastPos()[0])+cos(domain.getLastPos()[1]);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
            // Here fill the validation value for Df/Dx
            domain.template getLastProp<2>()[0] = 0;//cos(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<2>()[1] = 0;//-sin(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<4>()[0] = cos(domain.getLastPos()[0])*(sin(domain.getLastPos()[0])+sin(domain.getLastPos()[1]))+cos(domain.getLastPos()[1])*(cos(domain.getLastPos()[0])+cos(domain.getLastPos()[1]));
            domain.template getLastProp<4>()[1] = -sin(domain.getLastPos()[0])*(sin(domain.getLastPos()[0])+sin(domain.getLastPos()[1]))-sin(domain.getLastPos()[1])*(cos(domain.getLastPos()[0])+cos(domain.getLastPos()[1]));


//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain,2,rCut);
        Derivative_y Dy(domain,2,rCut);
        Gradient Grad(domain,2,rCut);
        Laplacian Lap(domain,2,rCut,3);
        Advection Adv(domain,2,rCut,3);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);
        auto dv = getV<3>(domain);

//        typedef boost::mpl::int_<std::is_fundamental<point_expression_op<Point<2U, double>, point_expression<double>, Point<2U, double>, 3>>::value>::blabla blabla;

//        std::is_fundamental<decltype(o1.value(key))>

        //vv=Lap(P);
        //dv=Lap(v);//+Dy(P);
        dv=Adv(v,v);//+Dy(P);
        auto it2 = domain.getDomainIterator();

        double worst1 = 0.0;

        while (it2.isNext())
        {
            auto p = it2.get();

            std::cout << "VALS: " << domain.getProp<3>(p)[0] << " " << domain.getProp<4>(p)[0] << std::endl;
            std::cout << "VALS: " << domain.getProp<3>(p)[1] << " " << domain.getProp<4>(p)[1] << std::endl;

            if (fabs(domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0]) > worst1)
            {
                worst1 = fabs(domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0]);
            }

            ++it2;
        }

        std::cout << "Maximum Error: " << worst1 << std::endl;

        domain.deleteGhost();
        domain.write("v");

        //std::cout << demangle(typeid(decltype(v)).name()) << "\n";

        //Debug<decltype(expr)> a;

        //typedef decltype(expr)::blabla blabla;

        //auto err = Dx + Dx;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE( dcpse_op_tests )
{
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size_t edgeSemiSize = 40;
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

    vector_dist<2, double, aggregate<double, double, double, VectorS<2,double>,VectorS<2,double>>> domain(0, box, bc, ghost);

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
        domain.getLastPos()[0] = x ;//+ gaussian(rng);
        mem_id k1 = key.get(1);
        double y = k1 * spacing[1];
        domain.getLastPos()[1] = y ;//+gaussian(rng);
        // Here fill the function value
        domain.template getLastProp<0>() = sin(domain.getLastPos()[0])+sin(domain.getLastPos()[1]);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
        // Here fill the validation value for Df/Dx
        domain.template getLastProp<2>() = cos(domain.getLastPos()[0])+cos(domain.getLastPos()[1]);
        domain.template getLastProp<4>()[0] = -sin(domain.getLastPos()[0]);
        domain.template getLastProp<4>()[1] = -sin(domain.getLastPos()[1]);


//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

        ++counter;
        ++it;
    }
    BOOST_TEST_MESSAGE("Sync domain across processors...");

    domain.map();
    domain.ghost_get<0>();

    Derivative_x Dx(domain,2,rCut);
    Derivative_y Dy(domain,2,rCut);
    Gradient Grad(domain,2,rCut);
    Laplacian Lap(domain,2,rCut);
	auto v = getV<1>(domain);
	auto P = getV<0>(domain);
	auto vv = getV<3>(domain);

	vv=Lap(P);
	v=Dx(P)+Dy(P);
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

	std::cout << "Maximum Error: " << worst << std::endl;

	domain.deleteGhost();
	domain.write("v");

	//std::cout << demangle(typeid(decltype(v)).name()) << "\n";

	//Debug<decltype(expr)> a;

	//typedef decltype(expr)::blabla blabla;

	//auto err = Dx + Dx;
}

BOOST_AUTO_TEST_SUITE_END()


