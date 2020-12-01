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
#include "util/SphericalHarmonics.hpp"


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


    BOOST_AUTO_TEST_CASE(Sph_harm) {
        BOOST_REQUIRE(openfpm::math::Y(2,1,0.5,0)+0.459674<0.00001);
        //These would be a requirement once Boost releases their fix
        //
        //BOOST_REQUIRE(boost::math::legendre_p(0,-1,1)=?);
        double nu=1.0;
        size_t grd_sz=20;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({-1.0, -1.0,-1.0}, {1.0,1.0,1.0});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = 2.0 / (sz[0] - 1);
        double rCut = 3.9*spacing;
        double R=1.0;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_B           RHS            V_t         P_anal              RHS2            Polar cord
        vector_dist_ws<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles(0, box, bc, ghost);


        auto &v_cl = create_vcluster();

//        openfpm::vector<aggregate<int>> bulk;
//        openfpm::vector<aggregate<int>> Surface;

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            auto key = it.get();
            double x = -1.0+key.get(0) * it.getSpacing(0);
            double y = -1.0+key.get(1) * it.getSpacing(1);
            double z = -1.0+key.get(2) * it.getSpacing(2);
            double r=sqrt(x*x+y*y+z*z);
            if (r<R-spacing/2.0) {
                Particles.add();
                Particles.getLastPos()[0] = x;
                Particles.getLastPos()[1] = y;
                Particles.getLastPos()[2] = z;
                Particles.getLastProp<8>()[0] = r;
                if (r==0){
                    Particles.getLastProp<8>()[1] = 0.0;
                }
                else{
                    Particles.getLastProp<8>()[1] = std::atan2(sqrt(x*x+y*y),z);
                }
                Particles.getLastProp<8>()[2] = std::atan2(y,x);
            }
            ++it;
        }

        int n_sp=int(grd_sz)*int(grd_sz)*3;

        double Golden_angle=M_PI * (3.0 - sqrt(5.0));

        for(int i=1;i<=n_sp;i++)
        {
            double y = 1.0 - (i /double(n_sp - 1.0)) * 2.0;
            double radius = sqrt(1 - y * y);
            double Golden_theta = Golden_angle * i;
            double x = cos(Golden_theta) * radius;
            double z = sin(Golden_theta) * radius;

            if (acos(z)==0 || acos(z)==M_PI){
                std::cout<<"Theta 0/Pi "<<std::endl;
                continue;
            }

            Particles.add();
            Particles.getLastPos()[0] = x;
            Particles.getLastPos()[1] = y;
            Particles.getLastPos()[2] = z;
            Particles.getLastProp<8>()[0] = 1.0 ;
            Particles.getLastProp<8>()[1] = std::atan2(sqrt(x*x+y*y),z);
            Particles.getLastProp<8>()[2] = std::atan2(y,x);
        }
        Particles.map();
        Particles.ghost_get<0>();


        std::unordered_map<const lm,double,key_hash,key_equal> Vr;
        std::unordered_map<const lm,double,key_hash,key_equal> V1;
        std::unordered_map<const lm,double,key_hash,key_equal> V2;
        //Setting max mode l_max
        constexpr int K = 2;
        //Setting amplitudes to 0
        for(int l=0;l<=K;l++){
            for(int m=-l;m<=l;m++){
                Vr[std::make_tuple(l,m)]=0.0;
                V1[std::make_tuple(l,m)]=0.0;
                V2[std::make_tuple(l,m)]=0.0;
            }


        }
        //Setting some amplitude for boundary velocity
        V1[std::make_tuple(1,0)]=1.0;

        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Point<3, double> xP = Particles.getProp<8>(p);
            Particles.getProp<0>(p) =0;
            if (xP[0]==1.0) {
//                Surface.add();
//                Surface.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) =  0;
                std::vector<double> SVel;
                SVel=openfpm::math::sumY<K>(xP[0],xP[1],xP[2],Vr,V1,V2);
                double SP=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Vr);
                Particles.getProp<2>(p)[0] = SVel[0];
                Particles.getProp<2>(p)[1] = SVel[1];
                Particles.getProp<2>(p)[2] = SVel[2];
                Particles.getProp<9>(p)[0] = SVel[0];
                Particles.getProp<9>(p)[1] = SVel[1];
                Particles.getProp<9>(p)[2] = SVel[2];
                Particles.getProp<5>(p) = SP;
                Particles.setSubset(p,1);

            }
            else {
//                bulk.add();
//                bulk.last().get<0>() = p.getKey();
                Particles.setSubset(p,0);
                Particles.getProp<0>(p) =  0;
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            ++it2;
        }

        vector_dist_subset<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles_bulk(Particles,0);
        vector_dist_subset<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles_surface(Particles,1);

        auto & bulk = Particles_bulk.getIds();
        auto & Surface = Particles_surface.getIds();

        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            Point<3, double> xp = Particles.getPos(p);
            Point<3, double> xP = Particles.getProp<8>(p);

            std::unordered_map<const lm,double,key_hash,key_equal> Ur;
            std::unordered_map<const lm,double,key_hash,key_equal> U2;
            std::unordered_map<const lm,double,key_hash,key_equal> U1;
            std::unordered_map<const lm,double,key_hash,key_equal> Plm;

            for (int l = 0; l <= K; l++) {
                for (int m = -l; m <= l; m++) {
                    auto Er= Vr.find(std::make_tuple(l,m));
                    auto E1= V1.find(std::make_tuple(l,m));
                    auto E2= V2.find(std::make_tuple(l,m));
                    std::vector<double> Sol=openfpm::math::sph_anasol_u(nu,l,m,Er->second,E1->second,E2->second,xP[0]);
                    Ur[std::make_tuple(l,m)]=Sol[0];
                    U1[std::make_tuple(l,m)]=Sol[1];
                    U2[std::make_tuple(l,m)]=Sol[2];
                    Plm[std::make_tuple(l,m)]=Sol[3];
                }

            }

            if(fabs(xP[0])>=1e-5 && xP[1]>1e-5 && (M_PI-xP[1])>=1e-5)
            {
                std::vector<double> SVel = openfpm::math::sumY<K>(xP[0], xP[1], xP[2], Ur, U1, U2);
                Particles.getProp<9>(p)[0] = SVel[0];
                Particles.getProp<9>(p)[1] = SVel[1];
                Particles.getProp<9>(p)[2] = SVel[2];
                Particles.getProp<5>(p) = openfpm::math::sumY_Scalar<K>(xP[0], xP[1], xP[2], Plm);
            }
        }

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_B = getV<2>(Particles);
        V.setVarId(0);
        auto DIV = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto P_anal = getV<5>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS = getV<7>(Particles);
        auto P_bulk = getV<0>(Particles_bulk);
        auto RHS_bulk = getV<7>(Particles_bulk);
        auto V_anal = getV<9>(Particles);

        V_t=V;
        P=0;
        P_bulk=0;
        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        double sampling=3.1;
        double sampling2=1.9;
        double rCut2=3.9*spacing;

        Derivative_x Dx(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dx(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dy(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dz(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_xx Dxx(Particles, 2, rCut2,sampling2,support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut2,sampling2,support_options::RADIUS);
        Derivative_zz Dzz(Particles, 2, rCut2,sampling2,support_options::RADIUS);

        //std::cout << "DCPSE KERNELS DONE" << std::endl;
        petsc_solver<double> solverPetsc;
        solverPetsc.setPreconditioner(PCNONE);
        timer tt;
        double sum=0,sum1=0;
        V_t=V;
        double V_err_eps = 1e-5;

        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 30;
        int ctr = 0, errctr, Vreset = 0;
        V_err = 1;
        n = 0;
        tt.start();
        while (V_err >= V_err_eps && n <= nmax) {
            //Particles.write_frame("StokesSphere",n);
            Particles.ghost_get<0>(SKIP_LABELLING);
            RHS_bulk[0] = B_Dx(P);
            RHS_bulk[1] = B_Dy(P);
            RHS_bulk[2] = B_Dz(P);
            DCPSE_scheme<equations3d3, decltype(Particles)> Solver(Particles);
            auto Stokes1 = nu * (Dxx(V[0])+Dyy(V[0])+Dzz(V[0]));
            auto Stokes2 = nu * (Dxx(V[1])+Dyy(V[1])+Dzz(V[1]));
            auto Stokes3 = nu * (Dxx(V[2])+Dyy(V[2])+Dzz(V[2]));
            Solver.impose(Stokes1, bulk, RHS[0], vx);
            Solver.impose(Stokes2, bulk, RHS[1], vy);
            Solver.impose(Stokes3, bulk, RHS[2], vz);
            Solver.impose(V[0], Surface, V_B[0], vx);
            Solver.impose(V[1], Surface, V_B[1], vy);
            Solver.impose(V[2], Surface, V_B[2], vz);
            Solver.solve_with_solver(solverPetsc, V[0], V[1], V[2]);
            //Solver.solve(V[0],V[1],V[2]);
            //std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<1>();
            DIV = -(Dx(V[0])+Dy(V[1])+Dz(V[2]));
            P_bulk = P + DIV;
            sum = 0;
            sum1 = 0;
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<4>(p)[0] - Particles.getProp<1>(p)[0]) *
                       (Particles.getProp<4>(p)[0] - Particles.getProp<1>(p)[0]) +
                       (Particles.getProp<4>(p)[1] - Particles.getProp<1>(p)[1]) *
                       (Particles.getProp<4>(p)[1] - Particles.getProp<1>(p)[1]) +
                       (Particles.getProp<4>(p)[2] - Particles.getProp<1>(p)[2]) *
                       (Particles.getProp<4>(p)[2] - Particles.getProp<1>(p)[2]);
                sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                        Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1] +
                        Particles.getProp<1>(p)[2] * Particles.getProp<1>(p)[2];
            }
            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            V_t = V;
            Particles.ghost_get<1>(SKIP_LABELLING);
            V_err_old = V_err;
            V_err = sum / sum1;
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-14) {
                errctr++;
            } else {
                errctr = 0;
            }
            if (n > 3) {
                if (errctr > 1) {
                    std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN ERROR" << std::endl;
                    Vreset = 1;
                    break;
                } else {
                    Vreset = 0;
                }
            }
            n++;

        }

        V_t=0;

        double worst=0;
        double L2=0;
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            Point<3,double> xP=Particles.getProp<8>(p);
            if(xP[0]>=1e-5 && xP[1]>1e-5 && (M_PI-xP[1])>=1e-5)
            {
                double dx=Particles.getProp<9>(p)[0] - Particles.getProp<1>(p)[0];
                double dy=Particles.getProp<9>(p)[1] - Particles.getProp<1>(p)[1];
                double dz=Particles.getProp<9>(p)[2] - Particles.getProp<1>(p)[2];
                Particles.getProp<4>(p)[0]=fabs(dx);
                Particles.getProp<4>(p)[1]=fabs(dy);
                Particles.getProp<4>(p)[2]=fabs(dz);
                L2 += dx*dx+dy*dy+dz*dz;
                if (std::max({fabs(dx),fabs(dy),fabs(dz)}) > worst) {
                    worst = std::max({fabs(dx),fabs(dy),fabs(dz)});
                }

            }
        }

        v_cl.sum(worst);
        v_cl.sum(L2);
        v_cl.execute();
/*        if (v_cl.rank() == 0) {
            std::cout<<"Gd,Surf,Bulk Size: "<<grd_sz<<","<<Surface.size()<<","<<bulk.size()<<std::endl;
            std::cout << "L2_Final: " <<sqrt(L2)<<","<<sqrt(L2/(bulk.size()+Surface.size()))
                      << std::endl;
            std::cout << "L_inf_Final: " << worst
                      << std::endl;
        }*/

        //Particles.write("StokesSphere");
        BOOST_REQUIRE(worst<1e-4);

    }

BOOST_AUTO_TEST_SUITE_END()


