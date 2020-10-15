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

//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests3)
    BOOST_AUTO_TEST_CASE(dcpse_op_vec3d) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 31;
        const size_t sz[3] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1,2 * edgeSemiSize+1};
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
            domain.getLastPos()[1] = z;//+gaussian(rng);
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

            domain.template getLastProp<4>()[0] = cos(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) +
                                                  cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
            domain.template getLastProp<4>()[1] = -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                                                  sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));
            domain.template getLastProp<4>()[2] = -sin(domain.getLastPos()[0]) * (sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1])) -
                                                  sin(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));

            domain.template getLastProp<5>()    = -(sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]) + sin(domain.getLastPos()[2])) ;


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

        std::cout << "Maximum Error in component 2: " << worst1 << std::endl;

        //Adv.checkMomenta(domain);
        //Adv.DrawKernel<2>(domain,0);

        //domain.deleteGhost();
        domain.write("v1");

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

        //std::cout << "Maximum Error: " << worst2 << std::endl;

        domain.deleteGhost();
        domain.write("v2");
        BOOST_REQUIRE(worst2 < 0.03);


    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_poisson_Robin_anal3d) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {161,161};
        Box<2, double> box({0, 0}, {0.5, 0.5});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3.1);
        double rCut = 3.1 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,double>> domain(0, box, bc, ghost);


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

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut,1.9,support_options::RADIUS);
        Derivative_y Dy(domain, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(domain, 2, rCut,1.9,support_options::RADIUS);



        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            //domain.getProp<3>(p)=1+xp[0]*xp[0]+2*xp[1]*xp[1];
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            }
            ++it2;
        }

        v=abs(DCPSE_sol-RHS);
        double worst1 = 0.0;

        it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(domain.getProp<1>(p) - domain.getProp<5>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<1>(p) - domain.getProp<5>(p));
            }
            ++it2;
        }
        BOOST_REQUIRE(worst1 < 0.03);




    }

    BOOST_AUTO_TEST_CASE(stokes_3d_petscP) {
        size_t grd_sz=21;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({0, 0,0}, {1,1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t   Helmholtz
        vector_dist<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,    double,VectorS<3, double>>> Particles(0, box, bc, ghost);
        vector_dist<3, double, aggregate<double,double,VectorS<3, double>>> Particles_subset(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            double z = key.get(2) * it.getSpacing(1);
            Particles.getLastPos()[2] = z;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> f_p;
        openfpm::vector<aggregate<int>> b_p;


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        V.setVarId(0);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        H.setVarId(0);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);
        auto P_bulk = getV<0>(Particles_subset);
        auto H_bulk = getV<1>(Particles_subset); //Pressure only on inside
        auto Grad_bulk = getV<2>(Particles_subset);

        // Here fill up the boxes for particle detection.

        Box<3, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1)- spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> down({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        Box<3, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)- spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)+ spacing / 2.0});

        Box<3, double> front({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                             {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getLow(2) + spacing / 2.0});

        Box<3, double> back({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getHigh(2) - spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(front);
        boxes.add(back);
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes_3d.vtk");
        auto it2 = Particles.getDomainIterator();
        Particles.ghost_get<0,1,2,3,4,5>();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;

            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (front.isInside(xp) == true) {
                f_p.add();
                f_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (back.isInside(xp) == true) {
                b_p.add();
                b_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            ++it2;
        }

        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
            Particles_subset.getLastPos()[2] = Particles.getPos(bulk.template get<0>(i))[2];
        }
        V_t=V;


        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS),B_Dx(Particles_subset, 2, rCut,1.9,support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS),B_Dy(Particles_subset, 2, rCut,1.9,support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,1.9,support_options::RADIUS),B_Dz(Particles_subset, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS),B_Lap(Particles_subset, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);


        petsc_solver<double> solverPetsc;
        solverPetsc.setRestart(250);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setRestart(250);


        double nu=1e-2;
        Particles.ghost_get<0,1,2,3,4,5,6,7>();
        double sum=0,sum2=0;
        int n=50;
        Particles.write_frame("Stokes3d",0);
        V_t=V;
        for(int i=1; i<=n ;i++)
        {
            Grad_bulk[0]=-B_Dx(P_bulk);
            Grad_bulk[1]=-B_Dy(P_bulk);
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles.template getProp<7>(bulk.template get<0>(i))[0] = Particles_subset.getProp<2>(i)[0];
                Particles.template getProp<7>(bulk.template get<0>(i))[1] = Particles_subset.getProp<2>(i)[1];
                Particles.template getProp<7>(bulk.template get<0>(i))[2] = Particles_subset.getProp<2>(i)[2];
            }
            DCPSE_scheme<equations3d3,decltype(Particles)> Solver(Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            auto Stokes3 = Adv(V[2],V_star[2])-nu*Lap(V_star[2]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(Stokes3,bulk,RHS2[2],vz);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[2], up_p,0,vz);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[2], r_p, 0,vz);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[2], dw_p,0,vz);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.impose(V_star[2], l_p, 0,vz);
            Solver.impose(V_star[0], f_p, 0,vx);
            Solver.impose(V_star[1], f_p, 0,vy);
            Solver.impose(V_star[2], f_p, 0,vz);
            Solver.impose(V_star[0], b_p, 0,vx);
            Solver.impose(V_star[1], b_p, 0,vy);
            Solver.impose(V_star[2], b_p, 0,vz);

            Solver.solve_with_solver(solverPetsc,V_star[0],V_star[1],V_star[2]);
            std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<2>();
            RHS=-Div(V_star);
            DCPSE_scheme<equations3d1,decltype(Particles)> SolverH(Particles);
            auto Helmholtz = Lap(H);
            auto D_x=Dx(H);
            auto D_y=Dy(H);
            auto D_z=Dz(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(H, up_p,0);
            SolverH.impose(H, r_p, 0);
            SolverH.impose(H, dw_p,0);
            SolverH.impose(H, l_p,0);
            SolverH.impose(H, f_p,0);
            SolverH.impose(H, b_p,0);

            SolverH.solve_with_solver(solverPetsc2,H);
            Particles.ghost_get<5>();
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_subset.getProp<1>(i) = Particles.template getProp<5>(bulk.template get<0>(i));
            }
            std::cout << "Helmholtz Solved" << std::endl;
            V=V_star+Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<f_p.size();j++)
            {   auto p=f_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<b_p.size();j++)
            {   auto p=b_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            P_bulk=P_bulk+B_Lap(H_bulk);//-0.5*B_Adv(V_t,H);
            for (int i = 0 ; i < bulk.size() ; i++) {

                Particles.template getProp<0>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            }
            Particles.ghost_get<0,1>();
            std::cout << "V,P Corrected" << std::endl;
            sum=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])+(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2])*(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1]+Particles.getProp<1>(p)[2]*Particles.getProp<1>(p)[2];
            }
            sum=sqrt(sum);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "Relative l2 convergence error = " <<sum/sum2<< std::endl;
            Particles.write_frame("Stokes3d",i);
        }
    }


    BOOST_AUTO_TEST_CASE(stokes_3d_petsc) {
        size_t grd_sz=31;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({0, 0,0}, {1,1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t   Helmholtz
        vector_dist<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,    double,VectorS<3, double>>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            double z = key.get(2) * it.getSpacing(1);
            Particles.getLastPos()[2] = z;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> f_p;
        openfpm::vector<aggregate<int>> b_p;


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        V.setVarId(0);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        H.setVarId(0);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);

        // Here fill up the boxes for particle detection.

        Box<3, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1)- spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> down({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        Box<3, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)- spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)+ spacing / 2.0});

        Box<3, double> front({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                             {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getLow(2) + spacing / 2.0});

        Box<3, double> back({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getHigh(2) - spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(front);
        boxes.add(back);
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes_3d.vtk");
        auto it2 = Particles.getDomainIterator();
        Particles.ghost_get<0,1,2,3,4,5>();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;

            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (front.isInside(xp) == true) {
                f_p.add();
                f_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (back.isInside(xp) == true) {
                b_p.add();
                b_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            ++it2;
        }
        V_t=V;


        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);


        petsc_solver<double> solverPetsc;
        solverPetsc.setRestart(250);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setRestart(250);


        double nu=1e-2;
        Particles.ghost_get<0,1,2,3,4,5,6,7>();
        double sum=0,sum2=0;
        int n=50;
        //Particles.write_frame("Stokes3d",0);
        V_t=V;
        for(int i=1; i<=n ;i++)
        {   RHS2=-Grad(P);
            DCPSE_scheme<equations3d3,decltype(Particles)> Solver(Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            auto Stokes3 = Adv(V[2],V_star[2])-nu*Lap(V_star[2]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(Stokes3,bulk,RHS2[2],vz);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[2], up_p,0,vz);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[2], r_p, 0,vz);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[2], dw_p,0,vz);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.impose(V_star[2], l_p, 0,vz);
            Solver.impose(V_star[0], f_p, 0,vx);
            Solver.impose(V_star[1], f_p, 0,vy);
            Solver.impose(V_star[2], f_p, 0,vz);
            Solver.impose(V_star[0], b_p, 0,vx);
            Solver.impose(V_star[1], b_p, 0,vy);
            Solver.impose(V_star[2], b_p, 0,vz);

            Solver.solve_with_solver(solverPetsc,V_star[0],V_star[1],V_star[2]);
            //std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<2>();
            RHS=-Div(V_star);
            DCPSE_scheme<equations3d1,decltype(Particles)> SolverH(Particles);
            auto Helmholtz = Lap(H);
            auto D_x=Dx(H);
            auto D_y=Dy(H);
            auto D_z=Dz(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(H, up_p,0);
            SolverH.impose(H, r_p, 0);
            SolverH.impose(H, dw_p,0);
            SolverH.impose(H, l_p,0);
            SolverH.impose(H, f_p,0);
            SolverH.impose(H, b_p,0);

            SolverH.solve_with_solver(solverPetsc2,H);
            Particles.ghost_get<5>();
            //std::cout << "Helmholtz Solved" << std::endl;
            V=V_star+Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<f_p.size();j++)
            {   auto p=f_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<b_p.size();j++)
            {   auto p=b_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            P=P+Lap(H)-0.5*Adv(V_t,H);
            Particles.ghost_get<0,1>();
            //std::cout << "V,P Corrected" << std::endl;
            sum=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])+(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2])*(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1]+Particles.getProp<1>(p)[2]*Particles.getProp<1>(p)[2];
            }
            sum=sqrt(sum);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "Relative l2 convergence error = " <<sum/sum2<< std::endl;
            return;
            Particles.write_frame("Stokes3d",i);

        }
    }

    BOOST_AUTO_TEST_CASE(stokes_3d) {
        size_t grd_sz=51;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({0, 0,0}, {1,1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t   Helmholtz
        vector_dist<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,    double,VectorS<3, double>>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            double z = key.get(2) * it.getSpacing(1);
            Particles.getLastPos()[2] = z;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> f_p;
        openfpm::vector<aggregate<int>> b_p;


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);

        // Here fill up the boxes for particle detection.

        Box<3, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1)- spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> down({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        Box<3, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)- spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)+ spacing / 2.0});

        Box<3, double> front({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getLow(2) + spacing / 2.0});

        Box<3, double> back({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getHigh(2) - spacing / 2.0},
                             {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(front);
        boxes.add(back);
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes_3d.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;

            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (front.isInside(xp) == true) {
                f_p.add();
                f_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (back.isInside(xp) == true) {
                b_p.add();
                b_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            ++it2;
        }
        V_t=V;


        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);

        double nu=1e-2;

        double sum=0,sum2=0;
        int n=50;
        Particles.write_frame("Stokes3d",0);
        V_t=V;
        for(int i=1; i<=n ;i++)
        {   RHS2=-Grad(P);
            DCPSE_scheme<equations3d3,decltype(Particles)> Solver(Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            auto Stokes3 = Adv(V[2],V_star[2])-nu*Lap(V_star[2]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(Stokes3,bulk,RHS2[2],vz);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[2], up_p,0,vz);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[2], r_p, 0,vz);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[2], dw_p,0,vz);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.impose(V_star[2], l_p, 0,vz);
            Solver.impose(V_star[0], f_p, 0,vx);
            Solver.impose(V_star[1], f_p, 0,vy);
            Solver.impose(V_star[2], f_p, 0,vz);
            Solver.impose(V_star[0], b_p, 0,vx);
            Solver.impose(V_star[1], b_p, 0,vy);
            Solver.impose(V_star[2], b_p, 0,vz);
            Solver.solve(V_star[0],V_star[1],V_star[2]);
            std::cout << "Stokes Solved" << std::endl;
            RHS=-Div(V_star);
            DCPSE_scheme<equations3d1,decltype(Particles)> SolverH(Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Lap(H);
            auto D_x=Dx(H);
            auto D_y=Dy(H);
            auto D_z=Dz(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(H, up_p,0);
            SolverH.impose(H, r_p, 0);
            SolverH.impose(H, dw_p,0);
            SolverH.impose(H, l_p,0);
            SolverH.impose(H, f_p,0);
            SolverH.impose(H, b_p,0);
            SolverH.solve(H);
            std::cout << "Helmholtz Solved" << std::endl;
            Particles.ghost_get<5>();
            V=V_star+Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<f_p.size();j++)
            {   auto p=f_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<b_p.size();j++)
            {   auto p=b_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            P=P+Lap(H)-0.5*Adv(V_t,H);
            std::cout << "V,P Corrected" << std::endl;
            sum=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])+(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2])*(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1]+Particles.getProp<1>(p)[2]*Particles.getProp<1>(p)[2];
            }
            sum=sqrt(sum);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "Relative l2 convergence error = " <<sum/sum2<< std::endl;
            Particles.write_frame("Stokes3d",i);

        }
    }

    BOOST_AUTO_TEST_CASE(stokes_3d_anal_petsc) {
        /*
      In 3D we use exact solution:

       u = x^2 + y^2
       v = y^2 + z^2
       w = x^2 + y^2 - 2(x+y)z
       p = x + y + z - 3/2
       f_x = f_y = f_z = 3

      -\Delta u + \nabla p + f = <-4, -4, -4> + <1, 1, 1> + <3, 3, 3> = 0
    \nabla \cdot u           = 2x + 2y - 2(x + y)                   = 0


    */
        size_t grd_sz=31;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({0, 0,0}, {1,1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t   Helmholtz
        vector_dist<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,    double,VectorS<3, double>>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            double z = key.get(2) * it.getSpacing(1);
            Particles.getLastPos()[2] = z;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> f_p;
        openfpm::vector<aggregate<int>> b_p;


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);

        // Here fill up the boxes for particle detection.

        Box<3, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1)- spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> down({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)+ spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,box.getHigh(2)- spacing / 2.0});

        Box<3, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        Box<3, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2)- spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2)+ spacing / 2.0});

        Box<3, double> front({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getLow(2) - spacing / 2.0},
                             {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getLow(2) + spacing / 2.0});

        Box<3, double> back({box.getLow(0) + spacing / 2.0, box.getLow(1) - spacing / 2.0,box.getHigh(2) - spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0,box.getHigh(2) + spacing / 2.0});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(front);
        boxes.add(back);
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes_3d.vtk");
        auto it2 = Particles.getDomainIterator();
        Particles.ghost_get<0,1,2,3,4,5>();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;

            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (front.isInside(xp) == true) {
                f_p.add();
                f_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else if (back.isInside(xp) == true) {
                b_p.add();
                b_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            ++it2;
        }
        V_t=V;


        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);

        double nu=1e-2;
        Particles.ghost_get<0,1,2,3,4,5,6,7>();
        double sum=0,sum2=0;
        int n=50;
        Particles.write_frame("Stokes3d",0);
        V_t=V;
        for(int i=1; i<=n ;i++)
        {   RHS2=-Grad(P);
            DCPSE_scheme<equations3d3,decltype(Particles)> Solver(Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            auto Stokes3 = Adv(V[2],V_star[2])-nu*Lap(V_star[2]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(Stokes3,bulk,RHS2[2],vz);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[2], up_p,0,vz);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[2], r_p, 0,vz);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[2], dw_p,0,vz);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.impose(V_star[2], l_p, 0,vz);
            Solver.impose(V_star[0], f_p, 0,vx);
            Solver.impose(V_star[1], f_p, 0,vy);
            Solver.impose(V_star[2], f_p, 0,vz);
            Solver.impose(V_star[0], b_p, 0,vx);
            Solver.impose(V_star[1], b_p, 0,vy);
            Solver.impose(V_star[2], b_p, 0,vz);
            petsc_solver<double> solverPetsc;
            solverPetsc.setRestart(500);
            Solver.solve_with_solver(solverPetsc,V_star[0],V_star[1],V_star[2]);
            std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<2>();
            RHS=-Div(V_star);
            DCPSE_scheme<equations3d1,decltype(Particles)> SolverH(Particles);
            auto Helmholtz = Lap(H);
            auto D_x=Dx(H);
            auto D_y=Dy(H);
            auto D_z=Dz(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(H, up_p,0);
            SolverH.impose(H, r_p, 0);
            SolverH.impose(H, dw_p,0);
            SolverH.impose(H, l_p,0);
            SolverH.impose(H, f_p,0);
            SolverH.impose(H, b_p,0);
            petsc_solver<double> solverPetsc2;
            solverPetsc2.setRestart(500);
            SolverH.solve_with_solver(solverPetsc2,H);
            Particles.ghost_get<5>();
            std::cout << "Helmholtz Solved" << std::endl;
            V=V_star+Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<f_p.size();j++)
            {   auto p=f_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<b_p.size();j++)
            {   auto p=b_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            P=P+Lap(H)-0.5*Adv(V_t,H);
            Particles.ghost_get<0,1>();
            std::cout << "V,P Corrected" << std::endl;
            sum=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])+(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2])*(Particles.getProp<4>(p)[2]- Particles.getProp<1>(p)[2]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1]+Particles.getProp<1>(p)[2]*Particles.getProp<1>(p)[2];
            }
            sum=sqrt(sum);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "Relative l2 convergence error = " <<sum/sum2<< std::endl;
            Particles.write_frame("Stokes3d",i);
        }
    }

BOOST_AUTO_TEST_SUITE_END()


