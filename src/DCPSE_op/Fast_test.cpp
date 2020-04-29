//
// Created by Abhinav Singh on 25.04.20.
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
BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)
    BOOST_AUTO_TEST_CASE(dcpse_Kernels) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 30;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {10.0, 10.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = box.getHigh(0) / (sz[0] - 1);
        spacing[1] = box.getHigh(1) / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3.0);
        double rCut = 2.0* spacing[0];
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
            domain.template getLastProp<1>()[1] = cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]);
            domain.template getLastProp<0>()= cos(domain.getLastPos()[0]) - sin(domain.getLastPos()[1]);

//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
            // Here fill the validation value for Df/Dx
            domain.template getLastProp<2>()[0] = 0;//cos(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<2>()[1] = 0;//-sin(domain.getLastPos()[0]);//+cos(domain.getLastPos()[1]);
            domain.template getLastProp<4>()[0] =
                    cos(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) +
                    cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
            domain.template getLastProp<4>()[1] =
                    -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                    sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));


//            domain.template getLastProp<2>() = 2 * x;
//            domain.template getLastProp<2>() = 1;

            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        //Derivative_x Dx(domain, 2, rCut);
        //Derivative_y Dy(domain, 2, rCut);
        //Gradient Grad(domain, 2, rCut);
        //Laplacian Lap(domain, 2, rCut, 3);
        //Advection Adv(domain, 3, rCut, 2);
        Divergence Div(domain, 2, rCut, 1.9,support_options::RADIUS);
        auto v = getV<1>(domain);
        auto anasol = getV<0>(domain);
        auto div = getV<5>(domain);

//        typedef boost::mpl::int_<std::is_fundamental<point_expression_op<Point<2U, double>, point_expression<double>, Point<2U, double>, 3>>::value>::blabla blabla;

//        std::is_fundamental<decltype(o1.value(key))>

        //vv=Lap(P);
        //dv=Lap(v);//+Dy(P);
        div = Div(v);//+Dy(P);
        auto it2 = domain.getDomainIterator();

        double worst1 = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();

            //std::cout << "VALS: " << domain.getProp<3>(p)[0] << " " << domain.getProp<4>(p)[0] << std::endl;
            //std::cout << "VALS: " << domain.getProp<3>(p)[1] << " " << domain.getProp<4>(p)[1] << std::endl;
            domain.getProp<6>(p)=fabs(domain.getProp<0>(p) - domain.getProp<5>(p));
            //domain.getProp<0>(p)=std::sqrt((domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0])*(domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0])+(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1])*(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1]));

            if (fabs(domain.getProp<0>(p) - domain.getProp<5>(p)) > worst1) {
                worst1 = fabs(domain.getProp<0>(p) - domain.getProp<5>(p));

            }
            ++it2;
        }
        std::cout << "Maximum Error: " << worst1 << std::endl;

        //Adv.checkMomenta(domain);
        //Adv.DrawKernel<2>(domain,0);

        domain.deleteGhost();
        domain.write("v");

        //std::cout << demangle(typeid(decltype(v)).name()) << "\n";

        //Debug<decltype(expr)> a;

        //typedef decltype(expr)::blabla blabla;

        //auto err = Dx + Dx;
    }



BOOST_AUTO_TEST_SUITE_END()
