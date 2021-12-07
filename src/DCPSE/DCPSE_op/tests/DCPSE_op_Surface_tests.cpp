//
// Created by Abhinav Singh on 15.11.21.
//

#include "config.h"
#ifdef HAVE_EIGEN
#ifdef HAVE_PETSC


#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../DCPSE_surface_op.hpp"
#include "../DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include <iostream>


BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)
    BOOST_AUTO_TEST_CASE(dcpse_surface_simple) {

        double boxP1{-1.5}, boxP2{1.5};
        double boxSize{boxP2 - boxP1};
        size_t n=81;
        size_t sz[2] = {n,n};
        double grid_spacing{boxSize/(sz[0]-1)};
        double rCut{3.9 * grid_spacing};

        Box<2,double> domain{{boxP1,boxP1},{boxP2,boxP2}};
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost{rCut + grid_spacing/8.0};
        auto &v_cl=create_vcluster();

        vector_dist_ws<2, double, aggregate<double,double,double[2],double,double[2]>> Sparticles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // 1. Particles on a line
        if (v_cl.rank() == 0) {
            for (int i = 0; i < n; ++i) {
                double xp = -1.5+i*grid_spacing;
                Sparticles.add();
                Sparticles.getLastPos()[0] = xp;
                Sparticles.getLastPos()[1] = 0;
                Sparticles.getLastProp<3>() = std::sin(xp);
                Sparticles.getLastProp<2>()[0] = 0;
                Sparticles.getLastProp<2>()[1] = 1.0;
                Sparticles.getLastProp<1>() = -std::sin(xp);//sin(theta)*exp(-finalT/(radius*radius));
                Sparticles.getLastSubset(0);
            }
        }
        Sparticles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Sparticles.ghost_get<0>();
        Sparticles.write("Sparticles");
        //Here template parameters are Normal property no.
        SurfaceDerivative_xx<2> SDxx(Sparticles, 2, rCut,grid_spacing);
        SurfaceDerivative_yy<2> SDyy(Sparticles, 2, rCut,grid_spacing);
        //SurfaceDerivative_x<2> SDx(Sparticles, 4, rCut,grid_spacing);
        //SurfaceDerivative_y<2> SDy(Sparticles, 4, rCut,grid_spacing);
        auto INICONC = getV<3>(Sparticles);
        auto CONC = getV<0>(Sparticles);
        auto TEMP = getV<4>(Sparticles);
        auto normal = getV<2>(Sparticles);
        //auto ANASOL = getV<1>(domain);
        CONC=SDxx(INICONC)+SDyy(INICONC);
        //TEMP[0]=(-normal[0]*normal[0]+1.0) * SDx(INICONC) - normal[0]*normal[1] * SDy(INICONC);
        //TEMP[1]=(-normal[1]*normal[1]+1.0) * SDy(INICONC) - normal[0]*normal[1] * SDx(INICONC);
        //Sparticles.ghost_get<4>();
        //CONC=SDxx(TEMP[0]) + SDyy(TEMP[1]);
        auto it2 = Sparticles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();

            if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            }

            ++it2;
        }
        Sparticles.deleteGhost();
        Sparticles.write("Sparticles");
        BOOST_REQUIRE(worst < 0.03);

}
    BOOST_AUTO_TEST_CASE(dcpse_surface_circle) {

        double boxP1{-1.5}, boxP2{1.5};
        double boxSize{boxP2 - boxP1};
        size_t n=512;
        size_t sz[2] = {n,n};
        double grid_spacing{boxSize/(sz[0]-1)};
        double rCut{3.9 * grid_spacing};

        Box<2,double> domain{{boxP1,boxP1},{boxP2,boxP2}};
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost{rCut + grid_spacing/8.0};
        auto &v_cl=create_vcluster();

        vector_dist_ws<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // Surface prameters
        const double radius{1.0};
        std::array<double,2> center{0.0,0.0};
        Point<2,double> coord;
        const double pi{3.14159265358979323846};

        // 1. Particles on surface
        double theta{0.0};
        double dtheta{2*pi/double(n)};
        if (v_cl.rank() == 0) {
            for (int i = 0; i < n; ++i) {
                coord[0] = center[0] + radius * std::cos(theta);
                coord[1] = center[1] + radius * std::sin(theta);
                Sparticles.add();
                Sparticles.getLastPos()[0] = coord[0];
                Sparticles.getLastPos()[1] = coord[1];
                Sparticles.getLastProp<3>() = std::sin(theta);
                Sparticles.getLastProp<2>()[0] = std::cos(theta);
                Sparticles.getLastProp<2>()[1] = std::sin(theta);
                Sparticles.getLastProp<1>() = -std::sin(theta);;//sin(theta)*exp(-finalT/(radius*radius));
                Sparticles.getLastSubset(0);
                theta += dtheta;
            }
        }
        Sparticles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Sparticles.ghost_get<0>();
        Sparticles.write("Sparticles");
        //Here template parameters are Normal property no.
        SurfaceDerivative_xx<2> SDxx(Sparticles, 2, rCut,grid_spacing);
        SurfaceDerivative_yy<2> SDyy(Sparticles, 2, rCut,grid_spacing);
        //SurfaceDerivative_xy<2> SDxy(Sparticles, 3, rCut,grid_spacing);
        //SurfaceDerivative_x<2> SDx(Sparticles, 3, rCut,grid_spacing);
        //SurfaceDerivative_y<2> SDy(Sparticles, 3, rCut,grid_spacing);
        auto INICONC = getV<3>(Sparticles);
        auto CONC = getV<0>(Sparticles);
        auto TEMP = getV<4>(Sparticles);
        auto normal = getV<2>(Sparticles);
        //auto ANASOL = getV<1>(domain);
        CONC=SDxx(INICONC)+SDyy(INICONC);
        //TEMP[0]=(-normal[0]*normal[0]+1.0) * SDx(INICONC) - normal[0]*normal[1] * SDy(INICONC);
        //TEMP[1]=(-normal[1]*normal[1]+1.0) * SDy(INICONC) - normal[0]*normal[1] * SDx(INICONC);
        //TEMP[0]=(-normal[0]*normal[0]+1.0);
        //TEMP[1]=normal[0]*normal[1];
        //Sparticles.ghost_get<2,4>();
        //CONC=SDx(TEMP[0]) + SDy(TEMP[1]);
        //CONC= (SDx(TEMP[0])*SDx(INICONC)+TEMP[0]*SDxx(INICONC)-(SDx(TEMP[1])*SDy(INICONC)+TEMP[1]*SDxy(INICONC)))+
        //        (SDy((-normal[1]*normal[1]+1.0))*SDy(INICONC)+(-normal[1]*normal[1]+1.0)*SDyy(INICONC)-(SDy(TEMP[1])*SDx(INICONC)+TEMP[1]*SDxy(INICONC)));
        auto it2 = Sparticles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            }
            ++it2;
        }
        Sparticles.deleteGhost();
        Sparticles.write("Sparticles");
        std::cout<<worst;
        BOOST_REQUIRE(worst < 0.03);
}
    BOOST_AUTO_TEST_CASE(dcpse_surface_solver_circle) {
        double boxP1{-1.5}, boxP2{1.5};
        double boxSize{boxP2 - boxP1};
        size_t n=512;
        size_t sz[2] = {n,n};
        double grid_spacing{boxSize/(sz[0]-1)};
        double rCut{3.9 * grid_spacing};

        Box<2,double> domain{{boxP1,boxP1},{boxP2,boxP2}};
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost{rCut + grid_spacing/8.0};
        auto &v_cl=create_vcluster();

        vector_dist_ws<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // Surface prameters
        const double radius{1.0};
        std::array<double,2> center{0.0,0.0};
        Point<2,double> coord;
        const double pi{3.14159265358979323846};

        // 1. Particles on surface
        double theta{0.0};
        double dtheta{2*pi/double(n)};
        if (v_cl.rank() == 0) {
            for (int i = 0; i < n; ++i) {
                coord[0] = center[0] + radius * std::cos(theta);
                coord[1] = center[1] + radius * std::sin(theta);
                Sparticles.add();
                Sparticles.getLastPos()[0] = coord[0];
                Sparticles.getLastPos()[1] = coord[1];
                Sparticles.getLastProp<3>() = -4*std::sin(2*theta);
                Sparticles.getLastProp<2>()[0] = std::cos(theta);
                Sparticles.getLastProp<2>()[1] = std::sin(theta);
                Sparticles.getLastProp<1>() = std::sin(2*theta);;//sin(theta)*exp(-finalT/(radius*radius));
                Sparticles.getLastSubset(0);
                if(coord[0]==1. && coord[1]==0.)
                {Sparticles.getLastSubset(1);}
                theta += dtheta;
            }
        }
        Sparticles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Sparticles.ghost_get<0>();
        Sparticles.write("Sparticles");
        vector_dist_subset<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles_bulk(Sparticles,0);
        vector_dist_subset<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles_boundary(Sparticles,1);
        auto & bulk=Sparticles_bulk.getIds();
        auto & boundary=Sparticles_boundary.getIds();
        //Here template parameters are Normal property no.
        SurfaceDerivative_xx<2> SDxx(Sparticles, 2, rCut,grid_spacing);
        SurfaceDerivative_yy<2> SDyy(Sparticles, 2, rCut,grid_spacing);
        auto INICONC = getV<3>(Sparticles);
        auto CONC = getV<0>(Sparticles);
        auto TEMP = getV<4>(Sparticles);
        auto normal = getV<2>(Sparticles);
        auto ANASOL = getV<1>(Sparticles);
        DCPSE_scheme<equations2d1,decltype(Sparticles)> Solver(Sparticles);
        Solver.impose(SDxx(CONC)+SDyy(CONC), bulk, INICONC);
        Solver.impose(CONC, boundary, ANASOL);
        Solver.solve(CONC);


        auto it2 = Sparticles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            }
            ++it2;
        }
        Sparticles.deleteGhost();
        Sparticles.write("Sparticles");
        std::cout<<worst;
        BOOST_REQUIRE(worst < 0.03);
}


BOOST_AUTO_TEST_SUITE_END()

#endif
#endif