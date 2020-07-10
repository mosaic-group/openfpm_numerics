//
// Created by Abhinav Singh on 01.05.20.
//


#include "config.h"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "Grid/grid_dist_id.hpp"
#include "FD_op.hpp"
#include "Grid/staggered_dist_grid.hpp"
//#include "FDScheme.hpp"
//#include "Solvers/petsc_solver.hpp"
//#include "DCPSE_op/EqnsStruct.hpp"

//! Specify the general characteristic of system to solve
/*struct equations2d1 {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef  grid_dist_id<2, double, aggregate<double, double, double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};*/

struct no_equation
{};

BOOST_AUTO_TEST_SUITE(fd_op_suite_tests)

    BOOST_AUTO_TEST_CASE(fd_op_tests) {
        size_t edgeSemiSize = 80;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        periodicity<2> bc({NON_PERIODIC, NON_PERIODIC});
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, long int> ghost(1);

        std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << std::endl;

        grid_dist_id<2, double, aggregate<double, double, double>> domain(sz, box,ghost,bc);

        BOOST_TEST_MESSAGE("Init domain...");
        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key_l = it.get();
            auto key = it.getGKey(key_l);
            mem_id i = key.get(0);
            double x = i * spacing[0];
            mem_id j = key.get(1);
            double y = j * spacing[1];
            // Here fill the function value P
            domain.template getProp<0>(key_l) = sin(x) + sin(y);
            domain.template getProp<1>(key_l) = 0;
            // Here fill the validation value for Df/Dx in property 3

            domain.template getProp<2>(key_l) = cos(x) + cos(y) + sin(x) + sin(y) + 5.0;

            ++it;
        }

        domain.ghost_get<0>();

        FD::Derivative_x Dx;
        FD::Derivative_y Dy;

        auto v = FD::getV<1>(domain);
        auto P = FD::getV<0>(domain);

        v = Dx(P) + Dy(P) + P + 5;

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

        domain.write("g");

        BOOST_REQUIRE(worst < 0.003);
    }


    BOOST_AUTO_TEST_CASE(lalpacian_test) {
        size_t edgeSemiSize = 80;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        periodicity<2> bc({NON_PERIODIC, NON_PERIODIC});
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, long int> ghost(1);

        std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << std::endl;

        grid_dist_id<2, double, aggregate<double, double, double>> domain(sz, box,ghost,bc);

        BOOST_TEST_MESSAGE("Init domain...");
        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key_l = it.get();
            auto key = it.getGKey(key_l);
            mem_id i = key.get(0);
            double x = M_PI / 2 + i * spacing[0];
            mem_id j = M_PI / 2 + key.get(1);
            double y = j * spacing[1];
            // Here fill the function value P
            domain.template getProp<0>(key_l) = sin(x) + sin(y);
            domain.template getProp<1>(key_l) = 0;
            // Here fill the validation value for Df/Dx in property 3

            domain.template getProp<2>(key_l) = -sin(x) -sin(y) + cos(x) + cos(y) + 5.0;

            ++it;
        }

        domain.ghost_get<0>();

        FD::Derivative_x Dx;
        FD::Derivative_y Dy;
        FD::Lap L;

        auto v = FD::getV<1>(domain);
        auto P = FD::getV<0>(domain);

        v = L(P) + Dx(P) + Dy(P) + 5;

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

        domain.write("g");

        BOOST_REQUIRE(worst < 0.003);
    }

    BOOST_AUTO_TEST_CASE(fd_op_tests_staggered) {

    	constexpr int phi_ = 0;
    	constexpr int ux_ = 1;

        size_t edgeSemiSize = 80;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        periodicity<2> bc({NON_PERIODIC, NON_PERIODIC});
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, long int> ghost(1);

        std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << std::endl;

        staggered_grid_dist<2, double, aggregate<double, double>> domain(sz, box,ghost,bc);

        comb<2> c({-1,-1});
        openfpm::vector<comb<2>> lc;
        lc.add(c);

        domain.setStagPosition<phi_>(lc);

        comb<2> c2({0,-1});
        lc.clear();
        lc.add(c2);

        domain.setStagPosition<ux_>(lc);

        BOOST_TEST_MESSAGE("Init domain...");
        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key_l = it.get();
            auto key = it.getGKey(key_l);
            mem_id i = key.get(0);
            double x = i * spacing[0];
            mem_id j = key.get(1);
            double y = j * spacing[1];
            // Here fill the function value P
            domain.template getProp<phi_>(key_l) = j % 2;

            ++it;
        }

        domain.ghost_get<0>();

        FD::Derivative_x Dx;
        FD::Derivative_y Dy;

        auto ux = FD::getV_stag<ux_>(domain);
        auto phi = FD::getV_stag<phi_>(domain);

        grid_key_dx<2> start({0ul,0u});
        grid_key_dx<2> stop({(long int)(sz[0]-2),(long int)(sz[1]-2)});

        auto it2 = domain.getSubDomainIterator(start,stop);
        while (it2.isNext())
        {
            auto key_l = it2.get();

            double test = phi.value(key_l,c2);

            BOOST_REQUIRE_EQUAL(test,0.5);

            ++it2;
        }
    }

    BOOST_AUTO_TEST_CASE(fd_op_tests_derivative_staggered) {

    	constexpr int x = 0;
    	constexpr int y = 1;

    	constexpr int phi_ = 0;
    	constexpr int ux_ = 1;

        size_t edgeSemiSize = 80;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        periodicity<2> bc({NON_PERIODIC, NON_PERIODIC});
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, long int> ghost(1);

        std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << std::endl;

        staggered_grid_dist<2, double, aggregate<double, double>> domain(sz, box,ghost,bc);

        comb<2> c({-1,-1});
        openfpm::vector<comb<2>> lc;
        lc.add(c);

        domain.setStagPosition<phi_>(lc);

        comb<2> c2({0,-1});
        lc.clear();
        lc.add(c2);

        domain.setStagPosition<ux_>(lc);

        BOOST_TEST_MESSAGE("Init domain...");
        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key_l = it.get();
            auto key = it.getGKey(key_l);
            mem_id i = key.get(0);
            double x = i * spacing[0];
            mem_id j = key.get(1);
            double y = j * spacing[1];
            // Here fill the function value P
            domain.template getProp<phi_>(key_l) = x*x + 3.0*y*y;

            ++it;
        }

        domain.ghost_get<0>();

        FD::Derivative_x Dx;
        FD::Derivative_y Dy;

        auto ux = FD::getV_stag<ux_>(domain);
        auto phi = FD::getV_stag<phi_>(domain);

        ux = Dx(phi) + Dy(phi) + phi + 5;

        grid_key_dx<2> start({1ul,1u});
        grid_key_dx<2> stop({(long int)(sz[0]-2),(long int)(sz[1]-2)});

        double worst_error = 0.0;

        auto it2 = domain.getSubDomainIterator(start,stop);
        while (it2.isNext())
        {
            auto key = it2.get();

            double test = domain.template getProp<ux_>(key);

            double dx = 0.5* (domain.template getProp<phi_>(key.move(x,1).move(y,1)) + domain.template getProp<phi_>(key.move(x,1)) - (domain.template getProp<phi_>(key.move(x,-1).move(y,1)) + domain.template getProp<phi_>(key.move(x,-1))) );
            double dy = 0.5* (domain.template getProp<phi_>(key.move(y,1).move(y,1)) + domain.template getProp<phi_>(key.move(y,1)) - (domain.template getProp<phi_>(key.move(y,-1).move(y,1)) + domain.template getProp<phi_>(key.move(y,-1))) );

            dx *= 0.5/spacing[0];
            dy *= 0.5/spacing[1];

            double phi_inte = 0.5*(domain.template getProp<phi_>(key.move(y,1)) + domain.template getProp<phi_>(key));

            if (fabs(test - (dx + dy + phi_inte + 5.0)) > worst_error )
            {worst_error = fabs(test - (dx + dy + phi_inte + 5.0));}

            ++it2;
        }

        BOOST_REQUIRE(worst_error < 5e-13);
    }


BOOST_AUTO_TEST_SUITE_END()
