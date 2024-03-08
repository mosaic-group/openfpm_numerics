//
// Created by Abhinav Singh on 15.11.21.
//

#include "config.h"
#ifdef HAVE_EIGEN
#ifdef HAVE_PETSC


#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

#define BOOST_TEST_DYN_LINK
#define HACKATHON
// #define SE_CLASS1

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "DCPSE/DCPSE_op/DCPSE_surface_op.hpp"
#include "DCPSE/DCPSE_op/DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include <iostream>

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)

    BOOST_AUTO_TEST_CASE(dcpse_surface_subset_circle) {
        // Surface prameters
        const double radius = 1.0;
        std::array<double,2> center{0.0,0.0};
        Point<2,double> coord;
        const double pi = 3.14159265358979323846;
        
        double boxP1 = -1.0, boxP2 = 1.0;
        double boxSize = boxP2 - boxP1;
        size_t n = 100;
        auto &v_cl = create_vcluster();
        size_t sz[2] = {n,n};
        double grid_spacing = 2.0*pi*radius/double(n);
        double rCut = 3.1 * grid_spacing;

        Box<2,double> domain({boxP1,boxP1},{boxP2,boxP2});
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost(rCut + grid_spacing/8.0); 
        vector_dist_ws<2, double, aggregate<double[2], double, double>> Particles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // Place particles
        // double theta = 0.0;
        // double dtheta = 2*pi/double(n);
        if (v_cl.rank() == 0) {
            // 1. Particles in bulk
            // auto it = Particles.getGridIterator(sz);
            // while (it.isNext()) {
            //     auto key = it.get();
            //     double x_pos = key.get(0) * grid_spacing + domain.getLow(0);
            //     double y_pos = key.get(1) * grid_spacing + domain.getLow(0);

            //     if ( (x_pos*x_pos + y_pos*y_pos) >= (radius - rCut*1.2)*(radius - rCut*1.2)){
            //         ++it;
            //         continue;
            //     }

            //     Particles.add();
            //     Particles.getLastPos()[0] = x_pos;
            //     Particles.getLastPos()[1] = y_pos;
            //     Particles.getLastProp<0>()[0] = 0;
            //     Particles.getLastProp<0>()[1] = 0;
            //     Particles.getLastProp<1>() = 0;
            //     Particles.getLastSubset(0);
            //     ++it;
            // }

            // 3. Particles on surface (outer circle)
            double theta = 0.0;
            double dtheta = 2*pi/double(n);
            for (int i = 0; i < n; ++i) {
                coord[0] = center[0] + radius * std::cos(theta);
                coord[1] = center[1] + radius * std::sin(theta);
                Particles.add();
                Particles.getLastPos()[0] = coord[0];
                Particles.getLastPos()[1] = coord[1];
                Particles.getLastProp<0>()[0] = std::cos(theta);
                Particles.getLastProp<0>()[1] = std::sin(theta);
                Particles.getLastProp<1>() = 0.0;
                Particles.getLastSubset(1);
                theta += dtheta;
            }

            theta = 0.0;
            dtheta = 2*pi/double(n);
            // 2. Particles on surface (inner circle)
            for (int i = 0; i < n; ++i) {
                coord[0] = center[0] + radius/2 * std::cos(theta);
                coord[1] = center[1] + radius/2 * std::sin(theta);
                Particles.add();
                Particles.getLastPos()[0] = coord[0];
                Particles.getLastPos()[1] = coord[1];
                Particles.getLastProp<0>()[0] = std::cos(theta);
                Particles.getLastProp<0>()[1] = std::sin(theta);
                Particles.getLastProp<1>() = 0.0;
                Particles.getLastSubset(0);
                theta += dtheta;
            }

            // // 2. Particles on surface (inner circle)
            // for (int i = 0; i < n; ++i) {
            //     coord[0] = center[0] + radius/2 * std::cos(theta);
            //     coord[1] = center[1] + radius/2 * std::sin(theta);
            //     Particles.add();
            //     Particles.getLastPos()[0] = coord[0];
            //     Particles.getLastPos()[1] = coord[1];
            //     Particles.getLastProp<0>()[0] = std::cos(theta);
            //     Particles.getLastProp<0>()[1] = std::sin(theta);
            //     Particles.getLastProp<1>() = 0.0;
            //     Particles.getLastSubset(0);
                
            //     coord[0] = center[0] + radius * std::cos(theta);
            //     coord[1] = center[1] + radius * std::sin(theta);
            //     Particles.add();
            //     Particles.getLastPos()[0] = coord[0];
            //     Particles.getLastPos()[1] = coord[1];
            //     Particles.getLastProp<0>()[0] = std::cos(theta);
            //     Particles.getLastProp<0>()[1] = std::sin(theta);
            //     Particles.getLastProp<1>() = 0.0;
            //     Particles.getLastSubset(1);
            //     theta += dtheta;
            // }
        }
        Particles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Particles.ghost_get<0>();

        // Declare subsets 
        vector_dist_subset<2, double, aggregate<double[2], double, double>> Particles_bulk(Particles,0);
        // vector_dist_subset<2, double, aggregate<double[2], double, double>> Particles_surface(Particles,1);

        std::cout << "after ghost_get: " << Particles_bulk.size_local_with_ghost() << " " << Particles.size_local_with_ghost() << " " << grid_spacing << " " << rCut << std::endl;
        
        Particles.write("Particles");
        //Here template parameters are Normal property no.
        // SurfaceDerivative_x<0> SDx(Particles, 2, rCut, grid_spacing);
        SurfaceDerivative_x<0> SDx(Particles, 2, rCut, 2.0*pi*radius/double(n));

        SurfaceDerivative_x<0> SDxBulk(Particles_bulk, 2, 3.1 * (2.0*pi*radius/(2*double(n))), 2.0*pi*radius/(2*double(n)));
        // SurfaceDerivative_x<0> SDxSurface(Particles_surface, 2, rCut, grid_spacing);



        std::cout << "1!!!!!!\n";

        // Derivative_x SDx(Particles_surface, 2, rCut, grid_spacing, support_options::RADIUS);
        // Derivative_y SDy(Particles_surface, 2, rCut, grid_spacing, support_options::RADIUS);

        auto normal = getV<0>(Particles);
        std::cout << "2!!!!!!\n";
        auto SDxEval = getV<1>(Particles);
        std::cout << "3!!!!!!\n";
        auto SDxSubsetEval = getV<2>(Particles);

        // curvature = (SDx(normal[0]) + SDy(normal[1]));

        // // check if curvature of bulk particles actually is zero, i.e. they were not used for surface DCPSE 
        std::cout << "4!!!!!!\n";
        SDxEval = SDx(normal[0]);
        std::cout << "5!!!!!!\n";
        SDxSubsetEval = SDxBulk(normal[0]);

        std::cout << "6!!!!!!\n";
        Particles.write("Particles2");

        bool check = true;
        auto it2 = Particles_bulk.getDomainIterator();
        std::cout << "before " << Particles_bulk.size_local() << " " << Particles.size_local() << " " << rCut << std::endl;
        while (it2.isNext()) {
            auto p = it2.get();

            std::cout << it2.get().getKey() << " " << it2.getOrig().getKey() << " " << Particles.getProp<1>(it2.getOrig()) << " " << Particles.getProp<2>(it2.getOrig()) << std::endl;
            if (Particles.getProp<1>(it2.getOrig()) != Particles.getProp<2>(it2.getOrig()))
                check = false;

            ++it2;
        }

        std::cout << " check " << check << std::endl;

        // // check if curvature of surface particles actually unequal zero
        // auto it3 = Particles_surface.getDomainIterator();
        // int count2 = 0;
        // while (it3.isNext()) {
        //     auto p = it3.get();
        //     if (Particles.getProp<1>(p) == 0.0) {
        //         ++count2;
        //     }
        //     ++it3;
        // }
        // Particles.deleteGhost();
        // //std::cout<<worst;
        // BOOST_REQUIRE_EQUAL(count, 0);
        BOOST_REQUIRE_EQUAL(check, 1);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
