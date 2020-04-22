/*
 * DCPSE_op_test.cpp
 *
 *  Created on: Feb 25, 2020
 *      Author: Abhinav Singh
 *
 */

#include "config.h"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "DCPSE_op.hpp"
#include "DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "EqnsStructPetsc.hpp"


//int vector_dist_expression_op<void,void,VECT_COPY_N_TO_N>::i = 0;
//int vector_dist_expression_op<void,void,VECT_COPY_1_TO_N>::i = 0;


//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests2)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_Lid_StokesPetsc) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t   Helmholtz
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,    double,VectorS<2, double>>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);

        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<3>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }
       // V_t=V;


        eq_id vx,vy;

        vx.setId(0);
        vy.setId(1);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);

        double nu=1e-2;

        double sum=0,sum2=0;
        int n=100;
        Particles.write_frame("Stokes",0);
        V_t=V;
        petsc_solver<double> solverPetsc;
        solverPetsc.setRestart(300);
        solverPetsc.setSolver(KSPGMRES);
        solverPetsc.setPreconditioner(PCKSP);


        petsc_solver<double> solverPetsc2;
        solverPetsc2.setRestart(300);
        solverPetsc2.setSolver(KSPGMRES);
        solverPetsc2.setPreconditioner(PCKSP);

        //Particles.ghost_get<0,1,2,3,4,5,6,7>();
        for(int i=1; i<=n ;i++)
        {   Particles.ghost_get<0,1>();
            RHS2=-Grad(P);
            DCPSE_scheme<equations2d2,decltype(Particles)> Solver( Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            return;
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.solve_with_solver(solverPetsc,V_star[0],V_star[1]);
            Particles.ghost_get<2>();

            //std::cout << "Stokes Solved" << std::endl;
            RHS=-Div(V_star);
            DCPSE_scheme<equations2d1,decltype(Particles)> SolverH( Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Lap(H);
            auto D_y=Dy(H);
            auto D_x=Dx(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(Dy(H), up_p,0);
            SolverH.impose(Dx(H), r_p, 0);
            SolverH.impose(Dy(H), dw_p,0);
            SolverH.impose(Dx(H), l_p,0);
            SolverH.solve_with_solver(solverPetsc2,H);

            Particles.ghost_get<5>();

            //std::cout << "Helmholtz Solved" << std::endl;
            V=V_star+Grad(H);

            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            Particles.ghost_get<4>();
            P=P+Lap(H)-0.5*Adv(V_t,H);

            //std::cout << "V,P Corrected" << std::endl;
            sum=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1];
            }
            sum=sqrt(sum);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "Rel l2 cgs err at "<<i<<"= " <<sum/sum2<< std::endl;
            if(i%5==0)
                Particles.write_frame("Stokes",i);
        }

    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_Lid_periodic) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(rCut);

        //                                  P        V                 v_star           RHS            V_t   Helmholtz                 RHS2
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,    double,VectorS<2, double>>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);

        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<3>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }
        V_t=V;


        eq_id vx,vy;

        vx.setId(0);
        vy.setId(1);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut, 1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut, 1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut, 1.9,support_options::RADIUS);


/*      starting the simulation at a nice *continuous* place
        V_t=1e-3*(1e-2*Lap(V)-Adv(V,V));
        RHS=Div(V_t);
        DCPSE_scheme<equations1d,decltype(Particles)> Solver( Particles,options_solver::LAGRANGE_MULTIPLIER);
        auto Pressure_Poisson = Lap(P);
        auto D_y=Dy(P);
        auto D_x=Dx(P);
        Solver.impose(Pressure_Poisson,bulk,prop_id<3>());
        Solver.impose(D_y, up_p,0);
        Solver.impose(D_x, r_p, 0);
        Solver.impose(-D_y, dw_p,0);
        Solver.impose(-D_x, l_p,0);
        Solver.solve(P);
        std::cout << "Poisson Solved" << std::endl;
        V_star = V + (V_t - 1e-3*Grad(P));
        V = V_star;*/

        double sum1=0,sum2=0;
        int n=10;
        double nu=1e-2;
        Particles.write_frame("Stokes",0);
        for(int i=1; i<=n ;i++)
        {   RHS2=-Grad(P);
            DCPSE_scheme<equations2d2,decltype(Particles)> Solver( Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.solve(V_star[0],V_star[1]);
            //std::cout << "Stokes Solved" << std::endl;
            RHS=Div(V_star);
            DCPSE_scheme<equations2d1,decltype(Particles)> SolverH( Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Lap(H);
            auto D_y=Dy(H);
            auto D_x=Dx(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(D_y, up_p,0);
            SolverH.impose(D_x, r_p, 0);
            SolverH.impose(-D_y, dw_p,0);
            SolverH.impose(-D_x, l_p,0);
            SolverH.solve(H);
            //std::cout << "Helmholtz Solved" << std::endl;
            V=V_star-Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            P=P+Lap(H);
            //std::cout << "V,P Corrected" << std::endl;
            sum1=0;
            sum2=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum1+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1];
            }
            sum1=sqrt(sum1);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "eps RMS=" <<sum1/sum2<< std::endl;
            Particles.write_frame("Stokes",i);

        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_Lid_Stokes) {
        const size_t sz[2] = {161,161};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t   Helmholtz
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,    double,VectorS<2, double>>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_star = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto H = getV<5>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS2 = getV<7>(Particles);

        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<3>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }
        V_t=V;


        eq_id vx,vy;

        vx.setId(0);
        vy.setId(1);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);

        double nu=1e-2;

        double sum=0,sum2=0;
        int n=100;
        Particles.write_frame("Stokes",0);
        V_t=V;
        for(int i=1; i<=n ;i++)
        {   RHS2=-Grad(P);
            DCPSE_scheme<equations2d2,decltype(Particles)> Solver( Particles);
            auto Stokes1 = Adv(V[0],V_star[0])-nu*Lap(V_star[0]);
            auto Stokes2 = Adv(V[1],V_star[1])-nu*Lap(V_star[1]);
            Solver.impose(Stokes1,bulk,RHS2[0],vx);
            Solver.impose(Stokes2,bulk,RHS2[1],vy);
            Solver.impose(V_star[0], up_p,1.0,vx);
            Solver.impose(V_star[1], up_p,0,vy);
            Solver.impose(V_star[0], r_p, 0,vx);
            Solver.impose(V_star[1], r_p, 0,vy);
            Solver.impose(V_star[0], dw_p,0,vx);
            Solver.impose(V_star[1], dw_p,0,vy);
            Solver.impose(V_star[0], l_p, 0,vx);
            Solver.impose(V_star[1], l_p, 0,vy);
            Solver.solve(V_star[0],V_star[1]);

            //std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<2>();
            RHS=-Div(V_star);

            DCPSE_scheme<equations2d1,decltype(Particles)> SolverH( Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Lap(H);
            auto D_y=Dy(H);
            auto D_x=Dx(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(Dy(H), up_p,0);
            SolverH.impose(Dx(H), r_p, 0);
            SolverH.impose(Dy(H), dw_p,0);
            SolverH.impose(Dx(H), l_p,0);
            SolverH.solve(H);

            Particles.write("Debug_out");
            return;

            //std::cout << "Helmholtz Solved" << std::endl;
            Particles.ghost_get<4,5>();
            V=V_star+Grad(H);

            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            P=P+Lap(H)-0.5*Adv(V_t,H);
            //std::cout << "V,P Corrected" << std::endl;
            sum=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<4>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<4>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<4>(p)[1]- Particles.getProp<1>(p)[1]);
                sum2+= Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1];
            }
            sum=sqrt(sum);
            sum2=sqrt(sum2);
            V_t=V;
            std::cout << "Rel l2 cgs err at "<<i<<"= " <<sum/sum2<< std::endl;
            if(i%5==0)
            Particles.write_frame("Stokes",i);
        }

    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_Lid_normal) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.5 * spacing;
        //                                  P        V                 Dv              RHS    Vtemp                   Proj_lap
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,double>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto dV = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Vtemp = getV<4>(Particles);
        auto divV = getV<5>(Particles);
        auto Proj_lap=getV<6>(Particles);


        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<3>(p) =0;
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }

        Derivative_x Dx(Particles, 2, 3.1*spacing,1.9,support_options::RADIUS);
        Derivative_y Dy(Particles, 2, 3.1*spacing,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, 3.1*spacing,1.9,support_options::RADIUS);
        Laplacian Lap(Particles, 2, 3.1*spacing,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, 3.1*spacing,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, 3.1*spacing,1.9,support_options::RADIUS);
        double dt=1e-3;
        int n=50;
        double nu=1e-2;
        dV=dt*(nu*Lap(V)-Adv(V,V));
        DCPSE_scheme<equations2d1,decltype(Particles)> Solver( Particles,options_solver::LAGRANGE_MULTIPLIER);
        auto Pressure_Poisson = Lap(P);
        auto D_y=Dy(P);
        auto D_x=Dx(P);
        Solver.impose(Pressure_Poisson,bulk,prop_id<3>());
        Solver.impose(D_y, up_p,0);
        Solver.impose(D_x, r_p, 0);
        Solver.impose(-D_y, dw_p,0);
        Solver.impose(-D_x, l_p,0);
        Solver.solve(P);
        std::cout << "Poisson Solved" << std::endl;
        Vtemp = V + (dV - dt*Grad(P));
        V = Vtemp;
        divV=Div(V);
        Particles.write_frame("Re100-3e-3-Lid",0);
        for(int i=1; i<=n ;i++)
        {   dV=dt*(nu*Lap(V)-Adv(V,V));
            RHS=1/dt*Div(dV);
            DCPSE_scheme<equations2d1,decltype(Particles)> Solver( Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Pressure_Poisson = Lap(P);
            auto D_y=Dy(P);
            auto D_x=Dx(P);
            Solver.impose(Pressure_Poisson,bulk,prop_id<3>());
            Solver.impose(D_y, up_p,0);
            Solver.impose(D_x, r_p, 0);
            Solver.impose(-D_y, dw_p,0);
            Solver.impose(-D_x, l_p,0);
            Solver.solve(P);
            Vtemp = V + (dV - dt*Grad(P));
            V = Vtemp;
            divV = Div(V);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            //if (i%10==0)
            Particles.write_frame("Re100-3e-3-Lid",i);
            std::cout<<i<<std::endl;
            if (i==100)
                dt=1e-4;
        }
    }


/*    BOOST_AUTO_TEST_CASE(dcpse_Lid_mirror) {
        const size_t sz[2] = {81,81};
        Box<2, double> box_mirror({-0.3, -0.3}, {1.3,1.3});
        Box<2, double> box_domain({0.0, 0.0}, {1.0,1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box_domain.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        //                                 P        V                 Dv              RHS    Vtemp
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles(0, box_mirror, bc, ghost);
//        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles_V(0, box, bc, ghost);
//        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles_P(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * spacing;
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * spacing;
            Particles.getLastPos()[1] = y;
            ++it;
        }
        Particles.map();
        Particles.ghost_get<0>();

        // Here fill up the boxes for particle detection.

        Box<2, double> up({box_domain.getLow(0) - spacing / 2.0, box_domain.getHigh(1) - spacing / 2.0},
                          {box_domain.getHigh(0) +  spacing / 2.0, box_domain.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box_domain.getLow(0)  - spacing / 2.0, box_domain.getLow(1) - spacing / 2.0},
                            {box_domain.getHigh(0) +  spacing / 2.0, box_domain.getLow(1) + spacing / 2.0});

        Box<2, double> left({box_domain.getLow(0) - spacing / 2.0, box_domain.getLow(1) + spacing / 2.0},
                            {box_domain.getLow(0) + spacing / 2.0, box_domain.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box_domain.getHigh(0) - spacing / 2.0, box_domain.getLow(1) + spacing / 2.0},
                             {box_domain.getHigh(0) + spacing / 2.0, box_domain.getHigh(1) - spacing / 2.0});

        *//*Box<2, double> up_m({box_domain.getLow(0) + 7*spacing / 2.0, box_domain.getHigh(1) - 5*spacing / 2.0},
                            {box_domain.getHigh(0) - 7* spacing / 2.0, box_domain.getHigh(1) + spacing / 2.0});

        Box<2, double> down_m({box_domain.getLow(0)  + 7*spacing / 2.0, box_domain.getLow(1) - spacing / 2.0},
                              {box_domain.getHigh(0) - 7* spacing / 2.0, box_domain.getLow(1) + 5*spacing / 2.0});

        Box<2, double> left_m({box_domain.getLow(0) - spacing / 2.0, box_domain.getLow(1) + 7*spacing / 2.0},
                              {box_domain.getLow(0) + 5*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> right_m({box_domain.getHigh(0) - 5*spacing / 2.0, box_domain.getLow(1) + 7*spacing / 2.0},
                               {box_domain.getHigh(0) + spacing / 2.0, box_domain.getHigh(1) - 7*spacing / 2.0});
*//*
        Box<2, double> up_mI({box_domain.getLow(0) + spacing / 2.0, box_domain.getHigh(1) - 7*spacing / 2.0},
                             {box_domain.getHigh(0) - spacing / 2.0, box_domain.getHigh(1) - spacing / 2.0});

        Box<2, double> down_mI({box_domain.getLow(0)  + spacing / 2.0, box_domain.getLow(1) + spacing / 2.0},
                               {box_domain.getHigh(0) - spacing / 2.0, box_domain.getLow(1) + 7*spacing / 2.0});

        Box<2, double> left_mI({box_domain.getLow(0) + spacing / 2.0, box_domain.getLow(1) + spacing / 2.0},
                               {box_domain.getLow(0) + 7*spacing / 2.0, box_domain.getHigh(1) - spacing / 2.0});

        Box<2, double> right_mI({box_domain.getHigh(0) - 7*spacing / 2.0, box_domain.getLow(1) + spacing / 2.0},
                                {box_domain.getHigh(0) - spacing / 2.0, box_domain.getHigh(1) - spacing / 2.0});


        Box<2, double> Bulk_box({box_domain.getLow(0) + spacing / 2.0, box_domain.getLow(1) + spacing / 2.0},
                                {box_domain.getHigh(0) - spacing / 2.0, box_domain.getHigh(1)  -spacing / 2.0});
        Box<2, double> bulk_box({box_domain.getLow(0) + 7*spacing / 2.0, box_domain.getLow(1) + 7*spacing / 2.0},
                                {box_domain.getHigh(0) - 7*spacing / 2.0, box_domain.getHigh(1) - 7*spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        *//*boxes.add(up_m);
        boxes.add(down_m);
        boxes.add(left_m);
        boxes.add(right_m);*//*
        boxes.add(up_mI);
        boxes.add(down_mI);
        boxes.add(left_mI);
        boxes.add(right_mI);
        boxes.add(Bulk_box) ;
        boxes.add(bulk_box);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("box.vtk");


        openfpm::vector<aggregate<int>> Bulk;
        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_pM;
        openfpm::vector<aggregate<int>> dw_pM;
        openfpm::vector<aggregate<int>> l_pM;
        openfpm::vector<aggregate<int>> r_pM;
        openfpm::vector<aggregate<int>> up_pMI;
        openfpm::vector<aggregate<int>> dw_pMI;
        openfpm::vector<aggregate<int>> l_pMI;
        openfpm::vector<aggregate<int>> r_pMI;
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> V_p;
        openfpm::vector<aggregate<int>> P_p;

        openfpm::vector<aggregate<int>> corner_ul;
        openfpm::vector<openfpm::vector<aggregate<int>>> corner_ref_ul;

        openfpm::vector<aggregate<int>> corner_ur;
        openfpm::vector<openfpm::vector<aggregate<int>>> corner_ref_ur;

        openfpm::vector<aggregate<int>> corner_dr;
        openfpm::vector<openfpm::vector<aggregate<int>>> corner_ref_dr;

        openfpm::vector<aggregate<int>> corner_dl;
        openfpm::vector<openfpm::vector<aggregate<int>>> corner_ref_dl;

        bool bulk_ref = false;
        bool Bulk_ref = false;

        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (up.isInside(xp) == true){
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                if (xp[0]==0 && xp[1]==0) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                }
                else {
                    dw_p.add();
                    dw_p.last().get<0>() = p.getKey();
                }
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            }
*//*            else if (up_m.isInside(xp) == true){
                up_pM.add();
                up_pM.last().get<0>() = p.getKey();
            }
            else if (down_m.isInside(xp) == true) {
                dw_pM.add();
                dw_pM.last().get<0>() = p.getKey();
            }
            else if (left_m.isInside(xp) == true) {
                l_pM.add();
                l_pM.last().get<0>() = p.getKey();
            }
            else if (right_m.isInside(xp) == true) {
                r_pM.add();
                r_pM.last().get<0>() = p.getKey();
            }*//*


            else if(Bulk_box.isInside(xp) == true) {
                Bulk.add();
                Bulk.last().get<0>() = p.getKey();
                if (up_mI.isInside(xp) == true){
                    if (left_mI.isInside(xp) == true)
                    {
                        corner_ref_ul.add();

                        Particles.add();
                        Particles.getLastPos()[0] = Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = 2 * box_domain.getHigh(1) - Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_ul.last().add();
                        corner_ref_ul.last().last().get<0>() = Particles.size_local() - 1;

                        Particles.add();
                        Particles.getLastPos()[0] = box_domain.getLow(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_ul.last().add();
                        corner_ref_ul.last().last().get<0>() = Particles.size_local() - 1;

                        corner_ul.add();
                        corner_ul.last().get<0>() = p.getKey();
                    }
                    else if (right_mI.isInside(xp) == true)
                    {
                        corner_ref_ur.add();

                        Particles.add();
                        Particles.getLastPos()[0] = Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = 2 * box_domain.getHigh(1) - Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_ur.last().add();
                        corner_ref_ur.last().last().get<0>() = Particles.size_local() - 1;

                        Particles.add();
                        Particles.getLastPos()[0] = 2 * box_domain.getHigh(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_ur.last().add();
                        corner_ref_ur.last().last().get<0>() = Particles.size_local() - 1;

                        corner_ur.add();
                        corner_ur.last().get<0>() = p.getKey();
                    }
                    else
                    {
                        Particles.add();
                        Particles.getLastPos()[0] = Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = 2 * box_domain.getHigh(1) - Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        up_pM.add();
                        up_pM.last().get<0>() = Particles.size_local()-1;

                        up_pMI.add();
                        up_pMI.last().get<0>() = p.getKey();
                    }
                }
                if (down_mI.isInside(xp) == true) {
                    if (left_mI.isInside(xp) == true)
                    {
                        corner_ref_dl.add();

                        Particles.add();
                        Particles.getLastPos()[0] = Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = box_domain.getLow(1) - Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_dl.last().add();
                        corner_ref_dl.last().last().get<0>() = Particles.size_local() - 1;

                        Particles.add();
                        Particles.getLastPos()[0] = box_domain.getLow(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_dl.last().add();
                        corner_ref_dl.last().last().get<0>() = Particles.size_local() - 1;

                        corner_dl.add();
                        corner_dl.last().get<0>() = p.getKey();
                    }
                    else if (right_mI.isInside(xp) == true)
                    {
                        corner_ref_dr.add();

                        Particles.add();
                        Particles.getLastPos()[0] = Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = box_domain.getLow(1) - Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_dr.last().add();
                        corner_ref_dr.last().last().get<0>() = Particles.size_local() - 1;

                        Particles.add();
                        Particles.getLastPos()[0] = 2 * box_domain.getHigh(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        corner_ref_dr.last().add();
                        corner_ref_dr.last().last().get<0>() = Particles.size_local() - 1;

                        corner_dr.add();
                        corner_dr.last().get<0>() = p.getKey();
                    }
                    else {

                        Particles.add();
                        Particles.getLastPos()[0] = Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = box_domain.getLow(1) - Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        dw_pM.add();
                        dw_pM.last().get<0>() = Particles.size_local()-1;

                        dw_pMI.add();
                        dw_pMI.last().get<0>() = p.getKey();
                    }
                }
                if (left_mI.isInside(xp) == true) {
                    if (down_mI.isInside(xp) == true)
                    {
*//*
                        Particles.add();
                        Particles.getLastPos()[0] = box_domain.getLow(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;
*//*
                    }
                    else if (up_mI.isInside(xp) == true)
                    {

                    }
                    else {
                        Particles.add();
                        Particles.getLastPos()[0] = box_domain.getLow(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        l_pM.add();
                        l_pM.last().get<0>() = Particles.size_local()-1;

                        l_pMI.add();
                        l_pMI.last().get<0>() = p.getKey();
                    }
                }
                if (right_mI.isInside(xp) == true) {
                    if (down_mI.isInside(xp) == true)
                    {

                    }
                    else if (up_mI.isInside(xp) == true)
                    {

                    }
                    else {

                        Particles.add();
                        Particles.getLastPos()[0] = 2 * box_domain.getHigh(0) - Particles.getPos(p)[0];
                        Particles.getLastPos()[1] = Particles.getPos(p)[1];

                        Particles.getLastProp<1>()[0] = 0.0;
                        Particles.getLastProp<1>()[1] = 0.0;

                        r_pM.add();
                        r_pM.last().get<0>() = Particles.size_local()-1;

                        r_pMI.add();
                        r_pMI.last().get<0>() = p.getKey();
                    }
                }
                if(bulk_box.isInside(xp)==true){
*//*                    if (bulk_ref == false)
                    {
                        bulk_ref = true;
                        ref_p.add();
                        ref_p.last().get<0>() = p.getKey();
                    }*//*
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                }
            }

            ++it2;
        }

        Particles.write("Particles");
        Derivative_x Dx(Particles, 2, rCut, 2);
        Derivative_y Dy(Particles, 2, rCut, 2);

        Gradient Grad(Particles, 2, rCut, 2);
        Laplacian Lap(Particles, 2, rCut, 2);
        Advection Adv(Particles, 2, rCut, 2);
        Divergence Div(Particles, 2, rCut, 2);

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto dV = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Vtemp = getV<4>(Particles);
        auto divV = getV<5>(Particles);

        double dt=0.005;
        int n=3;
        double nu=1e-2;

        Particles.write_frame("Re1000-1e-4-Lid",0);
        std::cout<<"Init Done"<<std::endl;

        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_up(up_pM, up_pMI);
        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_l(l_pM, l_pMI);
        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_r(r_pM, r_pMI);
        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_dw(dw_pM, dw_pMI);

        for(int i =1; i<=n ;i++)
        {   dV=V+dt*(nu*Lap(V)-Adv(V,V));
            std::cout<<"dV Done"<<std::endl;
            RHS=1.0/dt*Div(dV);

            std::cout<<"RHS Done"<<std::endl;
            DCPSE_scheme<equations2,decltype(Particles)> Solver(2*rCut, Particles);
            auto Pressure_Poisson = Lap(P);
            auto D_x=Dx(P);
            auto D_y=Dy(P);
            Solver.impose(Pressure_Poisson, Bulk, prop_id<3>());
            Solver.impose(D_y, up_p, prop_id<3>());
            Solver.impose(-D_y, dw_p, prop_id<3>());
            Solver.impose(-D_x, l_p, prop_id<3>());
            Solver.impose(D_x, r_p, prop_id<3>());
            Solver.impose(vcp_up, up_pM, 0);
            Solver.impose(vcp_l, l_pM,0);
            Solver.impose(vcp_r, r_pM, 0);
            Solver.impose(vcp_dw, dw_pM, 0);

            for (int i = 0 ; i < corner_ur.size() ; i++)
            {
                vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> cor(corner_ref_ur.template get<0>(i),corner_ur.template get<0>(i));

                Solver.impose(cor,corner_ref_ur.template get<0>(i), 0);
            }

            for (int i = 0 ; i < corner_ul.size() ; i++)
            {
                vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> cor(corner_ref_ul.template get<0>(i),corner_ul.template get<0>(i));

                Solver.impose(cor,corner_ref_ul.template get<0>(i), 0);
            }

            for (int i = 0 ; i < corner_dl.size() ; i++)
            {
                vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> cor(corner_ref_dl.template get<0>(i),corner_dl.template get<0>(i));

                Solver.impose(cor,corner_ref_dl.template get<0>(i), 0);
            }

            for (int i = 0 ; i < corner_dr.size() ; i++)
            {
                vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> cor(corner_ref_dr.template get<0>(i),corner_dr.template get<0>(i));

                Solver.impose(cor,corner_ref_dr.template get<0>(i), 0);
            }
            Solver.impose(P,ref_p,0);
            Solver.solve(P);
            std::cout<<"Poisson Solved"<<std::endl;
            divV=P-Lap(P);
            Particles.write_frame("Re1000-1e-4-Lid",i);
            return;
            Vtemp=dV-dt*Grad(P);
            V=Vtemp;
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            Particles.getProp<1>(0)[0] =  0;
            Particles.getProp<1>(0)[1] =  0;
            std::cout<<"Velocity updated"<<std::endl;
            Particles.write_frame("Re1000-1e-4-Lid",i);
            std::cout<<i<<std::endl;
        }

    }*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////














   /* BOOST_AUTO_TEST_CASE(dcpse_Lid_mirr) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        //                                 P        V                 Dv              RHS    Vtemp
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles(0, box, bc, ghost);
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles_V(0, box, bc, ghost);
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles_P(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }
        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> Bulk;
        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_pM;
        openfpm::vector<aggregate<int>> dw_pM;
        openfpm::vector<aggregate<int>> l_pM;
        openfpm::vector<aggregate<int>> r_pM;
        openfpm::vector<aggregate<int>> up_pMI;
        openfpm::vector<aggregate<int>> dw_pMI;
        openfpm::vector<aggregate<int>> l_pMI;
        openfpm::vector<aggregate<int>> r_pMI;
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> V_p;
        openfpm::vector<aggregate<int>> P_p;

        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + 5*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0},
                          {box.getHigh(0) - 5* spacing / 2.0, box.getHigh(1) - 5*spacing / 2.0});

        Box<2, double> down({box.getLow(0)  + 5*spacing / 2.0, box.getLow(1) + 5*spacing / 2.0},
                            {box.getHigh(0) - 5* spacing / 2.0, box.getLow(1) + 7*spacing / 2.0});

        Box<2, double> left({box.getLow(0) + 5*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                            {box.getLow(0) + 7*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - 7*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                             {box.getHigh(0) - 5*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> up_m({box.getLow(0) + 7*spacing / 2.0, box.getHigh(1) - 5*spacing / 2.0},
                            {box.getHigh(0) - 7* spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down_m({box.getLow(0)  + 7*spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) - 7* spacing / 2.0, box.getLow(1) + 5*spacing / 2.0});

        Box<2, double> left_m({box.getLow(0) - spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                            {box.getLow(0) + 5*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> right_m({box.getHigh(0) - 5*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> up_mI({box.getLow(0) + 7*spacing / 2.0, box.getHigh(1) - 13*spacing / 2.0},
                            {box.getHigh(0) - 7* spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> down_mI({box.getLow(0)  + 7*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                              {box.getHigh(0) - 7* spacing / 2.0, box.getLow(1) + 13*spacing / 2.0});

        Box<2, double> left_mI({box.getLow(0) + 7*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                              {box.getLow(0) + 13*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});

        Box<2, double> right_mI({box.getHigh(0) - 13*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                               {box.getHigh(0) - 7* spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});


        Box<2, double> Bulk_box({box.getLow(0) + 7*spacing / 2.0, box.getLow(1) + 7*spacing / 2.0},
                             {box.getHigh(0) - 7*spacing / 2.0, box.getHigh(1) - 7*spacing / 2.0});
        Box<2, double> bulk_box({box.getLow(0) + 13*spacing / 2.0, box.getLow(1) + 13*spacing / 2.0},
                                {box.getHigh(0) - 13*spacing / 2.0, box.getHigh(1) - 13*spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(up_m);
        boxes.add(down_m);
        boxes.add(left_m);
        boxes.add(right_m);
        boxes.add(up_mI);
        boxes.add(down_mI);
        boxes.add(left_mI);
        boxes.add(right_mI);
        boxes.add(Bulk_box) ;
        boxes.add(bulk_box);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("box.vtk");
        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (up.isInside(xp) == true){
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            }
            else if (up_m.isInside(xp) == true){
                up_pM.add();
                up_pM.last().get<0>() = p.getKey();
            }
            else if (down_m.isInside(xp) == true) {
                dw_pM.add();
                dw_pM.last().get<0>() = p.getKey();
            }
            else if (left_m.isInside(xp) == true) {
                l_pM.add();
                l_pM.last().get<0>() = p.getKey();
            }
            else if (right_m.isInside(xp) == true) {
                r_pM.add();
                r_pM.last().get<0>() = p.getKey();
            }
            else if(Bulk_box.isInside(xp) == true){
                Bulk.add();
                Bulk.last().get<0>() = p.getKey();
                if (up_mI.isInside(xp) == true){
                    up_pMI.add();
                    up_pMI.last().get<0>() = p.getKey();
                }
                if (down_mI.isInside(xp) == true) {
                    dw_pMI.add();
                    dw_pMI.last().get<0>() = p.getKey();
                }
                if (left_mI.isInside(xp) == true) {
                    l_pMI.add();
                    l_pMI.last().get<0>() = p.getKey();
                }
                if (right_m.isInside(xp) == true) {
                    r_pMI.add();
                    r_pMI.last().get<0>() = p.getKey();
                }
                if(bulk_box.isInside(xp)==true){
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                }
            }
            ++it2;
        }
        for (int i = 0 ; i < up_p.size() ; i++) {
            Particles_V.add();
            Particles_V.getLastPos()[0] = Particles.getPos(up_p.template get<0>(i))[0];
            Particles_V.getLastPos()[1] = Particles.getPos(up_p.template get<0>(i))[1];
            Particles_V.getLastProp<1>() = Particles.template getProp<1>(up_p.template get<0>(i));
            V_p.add();
            V_p.last().get<0>() = up_p.template get<0>(i);
        }

        for (int i = 0 ; i < dw_p.size() ; i++) {
            Particles_V.add();
            Particles_V.getLastPos()[0] = Particles.getPos(dw_p.template get<0>(i))[0];
            Particles_V.getLastPos()[1] = Particles.getPos(dw_p.template get<0>(i))[1];
            Particles_V.getLastProp<1>() = Particles.template getProp<1>(dw_p.template get<0>(i));
            V_p.add();
            V_p.last().get<0>() = dw_p.template get<0>(i);
        }
        for (int i = 0 ; i < l_p.size() ; i++) {
            Particles_V.add();
            Particles_V.getLastPos()[0] = Particles.getPos(l_p.template get<0>(i))[0];
            Particles_V.getLastPos()[1] = Particles.getPos(l_p.template get<0>(i))[1];
            Particles_V.getLastProp<1>() = Particles.template getProp<1>(l_p.template get<0>(i));
            V_p.add();
            V_p.last().get<0>() = l_p.template get<0>(i);
        }
        for (int i = 0 ; i < r_p.size() ; i++) {
            Particles_V.add();
            Particles_V.getLastPos()[0] = Particles.getPos(r_p.template get<0>(i))[0];
            Particles_V.getLastPos()[1] = Particles.getPos(r_p.template get<0>(i))[1];
            Particles_V.getLastProp<1>() = Particles.template getProp<1>(r_p.template get<0>(i));
            V_p.add();
            V_p.last().get<0>() = r_p.template get<0>(i);
        }

        for (int i = 0 ; i < Bulk.size() ; i++) {
            Particles_V.add();
            Particles_V.getLastPos()[0] = Particles.getPos(Bulk.template get<0>(i))[0];
            Particles_V.getLastPos()[1] = Particles.getPos(Bulk.template get<0>(i))[1];
            Particles_V.getLastProp<1>() = Particles.template getProp<1>(Bulk.template get<0>(i));
            V_p.add();
            V_p.last().get<0>() = Bulk.template get<0>(i);

            Particles_P.add();
            Particles_P.getLastPos()[0] = Particles.getPos(Bulk.template get<0>(i))[0];
            Particles_P.getLastPos()[1] = Particles.getPos(Bulk.template get<0>(i))[1];
            Particles_P.getLastProp<1>() = Particles.template getProp<1>(Bulk.template get<0>(i));
            P_p.add();
            P_p.last().get<0>() = Bulk.template get<0>(i);
        }

        for (int i = 0 ; i < up_pM.size() ; i++) {
            Particles_P.add();
            Particles_P.getLastPos()[0] = Particles.getPos(up_pM.template get<0>(i))[0];
            Particles_P.getLastPos()[1] = Particles.getPos(up_pM.template get<0>(i))[1];
            Particles_P.getLastProp<1>() = Particles.template getProp<1>(up_pM.template get<0>(i));
            P_p.add();
            P_p.last().get<0>() = Bulk.template get<0>(i);
        }
        for (int i = 0 ; i < dw_pM.size() ; i++) {
            Particles_P.add();
            Particles_P.getLastPos()[0] = Particles.getPos(dw_pM.template get<0>(i))[0];
            Particles_P.getLastPos()[1] = Particles.getPos(dw_pM.template get<0>(i))[1];
            Particles_P.getLastProp<1>() = Particles.template getProp<1>(dw_pM.template get<0>(i));
            P_p.add();
            P_p.last().get<0>() = Bulk.template get<0>(i);
        }
        for (int i = 0 ; i < l_pM.size() ; i++) {
            Particles_P.add();
            Particles_P.getLastPos()[0] = Particles.getPos(l_pM.template get<0>(i))[0];
            Particles_P.getLastPos()[1] = Particles.getPos(l_pM.template get<0>(i))[1];
            Particles_P.getLastProp<1>() = Particles.template getProp<1>(l_pM.template get<0>(i));
            P_p.add();
            P_p.last().get<0>() = Bulk.template get<0>(i);
        }
        for (int i = 0 ; i < r_pM.size() ; i++) {
            Particles_P.add();
            Particles_P.getLastPos()[0] = Particles.getPos(r_pM.template get<0>(i))[0];
            Particles_P.getLastPos()[1] = Particles.getPos(r_pM.template get<0>(i))[1];
            Particles_P.getLastProp<1>() = Particles.template getProp<1>(r_pM.template get<0>(i));
            P_p.add();
            P_p.last().get<0>() = Bulk.template get<0>(i);
        }
        Particles_V.map();
        Particles_V.ghost_get<0>();
        Particles_P.map();
        Particles_P.ghost_get<0>();
        Particles_V.write_frame("V",0);
        Particles_P.write_frame("P",0);

        Gradient Grad_P(Particles_P, 2, rCut, 2);
        Laplacian Lap_P(Particles_P, 2, rCut, 2),Lap_V(Particles_V, 2, rCut, 2);
        Advection Adv(Particles_V, 2, rCut, 2);
        Divergence Div_V(Particles_V, 2, rCut, 2);

        auto P = getV<0>(Particles_P);
        auto V = getV<1>(Particles_V);
        auto dV = getV<2>(Particles_V);
        auto RHS = getV<3>(Particles_V);
        auto Vtemp = getV<4>(Particles_V);
        auto divV = getV<5>(Particles_V);
        
        double dt=0.005;
        int n=3;
        double nu=1e-2;

        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_up(up_pM, up_pMI);
        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_l(l_pM, l_pMI);
        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_r(r_pM, r_pMI);
        vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_dw(dw_pM, dw_pMI);

        for (int i = 0 ; i < V_p.size() ; i++) {
            Particles_V.getProp<1>(i) = Particles.template getProp<1>(V_p.template get<0>(i));
        }
        Particles_V.write_frame("V",1);
        Particles.write_frame("Re1000-1e-4-Lid",0);
        std::cout<<"Init Done"<<std::endl;
        for(int i =1; i<=n ;i++)
        {   dV=V+dt*(nu*Lap_V(V)-Adv(V,V));
            std::cout<<"dV Done"<<std::endl;
            RHS=1.0/dt*Div_V(dV);
            for (int i = 0 ; i < V_p.size() ; i++) {
                Particles.getProp<3>(V_p.template get<0>(i)) = Particles_V.template getProp<3>(i);
            }
            for (int i = 0 ; i < P_p.size() ; i++) {
                Particles_P.getProp<3>(i) = Particles.template getProp<3>(P_p.template get<0>(i));
            }
            std::cout<<"RHS Done"<<std::endl;
            DCPSE_scheme<equations2,decltype(Particles_P)> Solver(2*rCut, Particles_P);
            auto Pressure_Poisson = Lap_P(P);
            Solver.impose(Pressure_Poisson, bulk, prop_id<3>());
            Solver.impose(vcp_up, up_pM, 0);
            Solver.impose(vcp_l, l_pM,0);
            Solver.impose(vcp_r, r_pM, 0);
            Solver.impose(vcp_dw, dw_pM, 0);
            Solver.solve(P);
            std::cout<<"Poisson Solved"<<std::endl;
            for (int i = 0 ; i < P_p.size() ; i++) {
                Particles.getProp<3>(i) = Particles_P.template getProp<3>(P_p.template get<0>(i));
            }
            Particles.write_frame("Re1000-1e-4-Lid",i);
            return;
*//*          k1=dt*(dV-Grad(P));
            Vtemp=V+k1*0.5;
            k2=dt*(nu*Lap(Vtemp)-Adv(Vtemp,Vtemp)-Grad(P));
            Vtemp=V+k2*0.5;
            k3=dt*(nu*Lap(Vtemp)-Adv(Vtemp,Vtemp)-Grad(P));
            Vtemp=V+k3;
            k4=dt*(nu*Lap(Vtemp)-Adv(Vtemp,Vtemp)-Grad(P));
            Vtemp=V+0.16666666666*(k1+2*k2+2*k3+k4);*//*
            Vtemp=dV-dt*Grad_P(P);
            V=Vtemp;
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            Particles.getProp<1>(0)[0] =  0;
            Particles.getProp<1>(0)[1] =  0;
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_P.getProp<1>(i) = Particles.template getProp<1>(bulk.template get<0>(i));
            }
//            divV_bulk=Div_bulk(V_bulk);
            for (int i = 0 ; i < bulk.size() ; i++)
            {
                Particles.template getProp<5>(bulk.template get<0>(i)) = Particles_P.getProp<0>(i);
            }
            //divV=Div(V);
            std::cout<<"Velocity updated"<<std::endl;
            Particles.write_frame("Re1000-1e-4-Lid",i);
            std::cout<<i<<std::endl;
        }
    }*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////











/*    BOOST_AUTO_TEST_CASE(dcpse_Lid2) {
        const size_t sz[2] = {31, 31};
        Box<2, double> box({0, 0}, {1, 1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        //                                  P        V                 Dv              RHS    Vtemp
        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>, double, VectorS<2, double>, double>> Particles(
                0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> ref_p2;
        openfpm::vector<aggregate<int>> FullBulk;
        openfpm::vector<aggregate<int>> UL_p;
        openfpm::vector<aggregate<int>> UL_p_ref;
        openfpm::vector<aggregate<int>> UR_p;
        openfpm::vector<aggregate<int>> UR_p_ref;
        openfpm::vector<aggregate<int>> LL_p;
        openfpm::vector<aggregate<int>> LL_p_ref;
        openfpm::vector<aggregate<int>> LR_p;
        openfpm::vector<aggregate<int>> LR_p_ref;

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto dV = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Vtemp = getV<4>(Particles);
        auto divV = getV<5>(Particles);


        // Here fill up the boxes for particle detection.

*//*
        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) + spacing / 2.0, box.getHigh(1) - 3*spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) - spacing / 2.0, box.getLow(1) + 3*spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) + spacing / 2.0, box.getLow(1) + 3 *spacing / 2.0},
                              {box.getLow(0) + 3*spacing / 2.0, box.getHigh(1) - 3*spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - 3*spacing / 2.0, box.getLow(1) + 3*spacing / 2.0},
                               {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - 3*spacing / 2.0});
*//*
        Box<2, double> up({box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) - 3 * spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                            {box.getHigh(0) - 3 * spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> down({box.getLow(0) + 3 * spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) + 3 * spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) + spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                              {box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                               {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> UL({box.getLow(0) - spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                          {box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> UL_ref({box.getLow(0) + spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                              {box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - spacing / 2.0});


        Box<2, double> UR({box.getHigh(0) - 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> UR_ref({box.getHigh(0) - 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                              {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> LL({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                          {box.getLow(0) + 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> LL_ref({box.getLow(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getLow(0) + 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> LR({box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> LR_ref({box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) - spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(up_d);
        boxes.add(down);
        boxes.add(down_u);
        boxes.add(left);
        boxes.add(left_r);
        boxes.add(right);
        boxes.add(right_l);
        boxes.add(UL);
        boxes.add(UL_ref);
        boxes.add(UR_ref);
        boxes.add(UR);
        boxes.add(LL_ref);
        boxes.add(LL);
        boxes.add(LR_ref);
        boxes.add(LR);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (xp[0] != 0 || xp[1] != 0) {
                if (xp[0] != 0 + spacing || xp[1] != 0) {
                    FullBulk.add();
                    FullBulk.last().get<0>() = p.getKey();
                }
            }
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] = 1;
                Particles.getProp<1>(p)[1] = 0;
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            } else if (up_d.isInside(xp) == true) {
                up_p1.add();
                up_p1.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            } else if (down_u.isInside(xp) == true) {
                dw_p1.add();
                dw_p1.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();

            } else if (left_r.isInside(xp) == true) {
                l_p1.add();
                l_p1.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            } else if (right_l.isInside(xp) == true) {
                r_p1.add();
                r_p1.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            } else if (UL.isInside(xp) == true) {
                if (UL_ref.isInside(xp) == true) {
                    UL_p_ref.add();
                    UL_p_ref.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                } else {
                    UL_p.add();
                    UL_p.last().get<0>() = p.getKey();
                }
            } else if (UR.isInside(xp) == true) {
                if (UR_ref.isInside(xp) == true) {
                    UR_p_ref.add();
                    UR_p_ref.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                } else {
                    UR_p.add();
                    UR_p.last().get<0>() = p.getKey();
                }
            } else if (LR.isInside(xp) == true) {
                if (LR_ref.isInside(xp) == true) {
                    LR_p_ref.add();
                    LR_p_ref.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                } else {
                    LR_p.add();
                    LR_p.last().get<0>() = p.getKey();
                }
            } else if (LL.isInside(xp) == true) {
                if (LL_ref.isInside(xp) == true) {
                    LL_p_ref.add();
                    LL_p_ref.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                } else if (xp[0] == 0 && xp[1] == 0) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                } else {
                    LL_p.add();
                    LL_p.last().get<0>() = p.getKey();
                }
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) = 0.0;
            }

            ++it2;
        }
        Particles.getProp<1>(959)[1] = 0;
        Particles.getProp<1>(931)[0] = 1;

        *//* for(int j=0;j<up_p1.size();j++)
         {
             auto p1=up_p1.get<0>(j);
             Point<2, double> xp = Particles.getPos(p1);
             Particles.getProp<3>(p1) =  xp[0];
         }
         for(int j=0;j<l_p1.size();j++)
         {
             auto p1=l_p1.get<0>(j);
             Point<2, double> xp = Particles.getPos(p1);
             Particles.getProp<3>(p1) =  xp[1];
         }
         for(int j=0;j<r_p1.size();j++)
         {
             auto p1=r_p1.get<0>(j);
             Point<2, double> xp = Particles.getPos(p1);
             Particles.getProp<3>(p1) =  xp[1];
         }
         for(int j=0;j<dw_p1.size();j++)
         {
             auto p1=dw_p1.get<0>(j);
             Point<2, double> xp = Particles.getPos(p1);
             Particles.getProp<3>(p1)=  xp[0];
         }*//*
        Particles.write_frame("Re1000-1e-4-Lid", 0);
*//*
        for(int j=0;j<up_p1.size();j++)
        { if(j==0)
            {   auto p=up_p.get<0>(j);
                auto p1=up_p1.get<0>(j);
                auto p2=l_p.get<0>(l_p.size()-1);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
            auto p=up_p.get<0>(j+1);
            auto p1=up_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
            if(j==up_p1.size()-1)
            {   auto p=up_p.get<0>(j+2);
                auto p1=up_p1.get<0>(j);
                auto p2=r_p.get<0>(r_p.size()-1);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
        }
        for(int j=0;j<l_p1.size();j++)
        {
            auto p=l_p.get<0>(j+1);
            auto p1=l_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
        }
        for(int j=0;j<r_p1.size();j++)
        {
            auto p=r_p.get<0>(j+1);
            auto p1=r_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
        }
        for(int j=0;j<dw_p1.size();j++)
        {   auto p=dw_p.get<0>(j);
            auto p1=dw_p1.get<0>(j);
            Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
            if(j==0)
            {
                auto p = ref_p.get<0>(j);
                auto p2=l_p.get<0>(0);
                Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
            }
            if(j==dw_p1.size()-1)
            {   auto p=dw_p.get<0>(j+1);
                auto p2=r_p.get<0>(0);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
            }
        }*//*
        Particles.write_frame("Re1000-1e-4-Lid", 1);
        Derivative_x Dx(Particles, 2, rCut, 2);
        Derivative_y Dy(Particles, 2, rCut, 2);
        Gradient Grad(Particles, 2, rCut, 2);
        Laplacian Lap(Particles, 2, rCut, 2);
        Advection Adv(Particles, 2, rCut, 2);
        Divergence Div(Particles, 2, rCut, 2);
        double dt = 5e-4;
        int n = 3;
        double nu = 1e-2;
*//*        divV=Div(V);
        DCPSE_scheme<equations2,decltype(Particles)> Solver(2*rCut, Particles);
        auto Pressure_Poisson = Lap(P);
        auto D_y=Dy(P);
        auto D_x=Dx(P);
        Solver.impose(Pressure_Poisson,bulk,prop_id<5>());
        Solver.impose(D_y, up_p,0);
        Solver.impose(D_y, dw_p,0);
        Solver.impose(D_x, l_p,0);
        Solver.impose(D_x, r_p, 0);
        Solver.impose(P, ref_p, 0);
        Solver.solve(P);

        std::cout<<"Poisson Solved"<<std::endl;
        V=V-Grad(P);
        Particles.write_frame("Re1000-1e-4-Lid",0);
        std::cout<<"Init Done"<<std::endl;
        return;*//*
        for (int i = 2; i <= n; i++) {
            dV = dt * (nu * Lap(V) - Adv(V, V));
            RHS = 1 / dt * Div(dV);
            std::cout << "RHS Done" << std::endl;

            //Neumann Copy
*//*
            for(int j=0;j<up_p1.size();j++)
            { if(j==0)
                {   auto p=up_p.get<0>(j);
                    auto p1=up_p1.get<0>(j);
                    auto p2=l_p.get<0>(l_p.size()-1);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                    Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
                }
                auto p=up_p.get<0>(j+1);
                auto p1=up_p1.get<0>(j);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                if(j==up_p1.size()-1)
                {   auto p=up_p.get<0>(j+2);
                    auto p1=up_p1.get<0>(j);
                    auto p2=r_p.get<0>(r_p.size()-1);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                    Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
                }
            }
            for(int j=0;j<l_p1.size();j++)
            {
                auto p=l_p.get<0>(j+1);
                auto p1=l_p1.get<0>(j);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
            }
            for(int j=0;j<r_p1.size();j++)
            {
                auto p=r_p.get<0>(j+1);
                auto p1=r_p1.get<0>(j);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
            }
            for(int j=0;j<dw_p1.size();j++)
            {   auto p=dw_p.get<0>(j);
                auto p1=dw_p1.get<0>(j);
                Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                if(j==0)
                {
                    auto p = ref_p.get<0>(j);
                    auto p2=l_p.get<0>(0);
                    Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                    Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
                }
                if(j==dw_p1.size()-1)
                {   auto p=dw_p.get<0>(j+1);
                    auto p2=r_p.get<0>(0);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                    Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
                }
            }*//*

            vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_up(up_p, up_p1);
            vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_l(l_p, l_p1);
            vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_r(r_p, r_p1);
            vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> vcp_dw(dw_p, dw_p1);

            vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> vcp_ul(UL_p, UL_p_ref.get<0>(0));
            vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> vcp_ur(UR_p, UR_p_ref.get<0>(0));
            vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> vcp_dl(LL_p, LL_p_ref.get<0>(0));
            vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> vcp_dr(LR_p, LR_p_ref.get<0>(0));

            DCPSE_scheme<equations2, decltype(Particles)> Solver( Particles);
            auto Pressure_Poisson = Lap(P);
            Solver.impose(Pressure_Poisson, bulk, prop_id<3>());
            Solver.impose(vcp_up, up_p, 0);
            Solver.impose(vcp_l, l_p, 0);
            Solver.impose(vcp_r, r_p, 0);
            Solver.impose(vcp_dw, dw_p, 0);
            Solver.impose(vcp_ul, UL_p, 0);
            Solver.impose(vcp_ur, UR_p, 0);
            Solver.impose(vcp_dl, LL_p, 0);
            Solver.impose(vcp_dr, LR_p, 0);
*//*            Solver.impose(P, up_p, up_p1);
            Solver.impose(P, dw_p,prop_id<3>());
            Solver.impose(P, l_p,prop_id<3>());
            Solver.impose(P, r_p, prop_id<3>());*//*
            Solver.impose(P, ref_p, 0);
            Solver.solve(P);
            std::cout << "Poisson Solved" << std::endl;
            Vtemp = V + dV - dt * Grad(P);
            V = Vtemp;
            divV = Div(V);
            for (int j = 0; j < up_p.size(); j++) {
                auto p = up_p.get<0>(j);
                Particles.getProp<1>(p)[0] = 1;
                Particles.getProp<1>(p)[1] = 0;
            }
            for (int j = 0; j < l_p.size(); j++) {
                auto p = l_p.get<0>(j);
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            }
            for (int j = 0; j < r_p.size(); j++) {
                auto p = r_p.get<0>(j);
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            }
            for (int j = 0; j < dw_p.size(); j++) {
                auto p = dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] = 0;
                Particles.getProp<1>(p)[1] = 0;
            }
            Particles.getProp<1>(0)[0] = 0;
            Particles.getProp<1>(0)[1] = 0;
            Particles.getProp<1>(959)[1] = 0;
            Particles.getProp<1>(931)[0] = 1;
            std::cout << "Velocity updated" << std::endl;
            Particles.write_frame("Re1000-1e-4-Lid", i);
            std::cout << i << std::endl;
        }
    }*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





/*
        BOOST_AUTO_TEST_CASE(dcpse_Lid_copy_RHS) {
            const size_t sz[2] = {31,31};
            Box<2, double> box({0, 0}, {1,1});
            size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
            double spacing = box.getHigh(0) / (sz[0] - 1);
            Ghost<2, double> ghost(spacing * 3);
            double rCut = 2.0 * spacing;
            //                                  P        V                 Dv              RHS    Vtemp
            vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double>> Particles(0, box, bc, ghost);
            auto it = Particles.getGridIterator(sz);
            while (it.isNext()) {
                Particles.add();
                auto key = it.get();
                double x = key.get(0) * it.getSpacing(0);
                Particles.getLastPos()[0] = x;
                double y = key.get(1) * it.getSpacing(1);
                Particles.getLastPos()[1] = y;
                ++it;
            }

            Particles.map();
            Particles.ghost_get<0>();

            openfpm::vector<aggregate<int>> bulk;
            openfpm::vector<aggregate<int>> up_p;
            openfpm::vector<aggregate<int>> dw_p;
            openfpm::vector<aggregate<int>> l_p;
            openfpm::vector<aggregate<int>> r_p;
            openfpm::vector<aggregate<int>> up_p1;
            openfpm::vector<aggregate<int>> dw_p1;
            openfpm::vector<aggregate<int>> l_p1;
            openfpm::vector<aggregate<int>> r_p1;
            openfpm::vector<aggregate<int>> ref_p;
            openfpm::vector<aggregate<int>> ref_p2;
            openfpm::vector<aggregate<int>> FullBulk;
            *//*openfpm::vector<aggregate<int>> UL_p;
            openfpm::vector<aggregate<int>> UL_p_ref;
            openfpm::vector<aggregate<int>> UR_p;
            openfpm::vector<aggregate<int>> UR_p_ref;
            openfpm::vector<aggregate<int>> LL_p;
            openfpm::vector<aggregate<int>> LL_p_ref;
            openfpm::vector<aggregate<int>> LR_p;
            openfpm::vector<aggregate<int>> LR_p_ref;*//*

            auto P = getV<0>(Particles);
            auto V = getV<1>(Particles);
            auto dV = getV<2>(Particles);
            auto RHS = getV<3>(Particles);
            auto Vtemp = getV<4>(Particles);
            auto divV = getV<5>(Particles);


            // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) + spacing / 2.0, box.getHigh(1) - 3*spacing / 2.0},
                            {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) - spacing / 2.0, box.getLow(1) + 3*spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) + spacing / 2.0, box.getLow(1) + 3 *spacing / 2.0},
                              {box.getLow(0) + 3*spacing / 2.0, box.getHigh(1) - 3*spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - 3*spacing / 2.0, box.getLow(1) + 3*spacing / 2.0},
                               {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - 3*spacing / 2.0});


            openfpm::vector<Box<2, double>> boxes;
            boxes.add(up);
            boxes.add(up_d);
            boxes.add(down);
            boxes.add(down_u);
            boxes.add(left);
            boxes.add(left_r);
            boxes.add(right);
            boxes.add(right_l);
            VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
            vtk_box.add(boxes);
            vtk_box.write("boxes.vtk");
            auto it2 = Particles.getDomainIterator();

            while (it2.isNext()) {
                auto p = it2.get();
                Point<2, double> xp = Particles.getPos(p);
                if (xp[0]!=0 || xp[1]!=0){
                        FullBulk.add();
                        FullBulk.last().get<0>() = p.getKey();
               }
                if (up.isInside(xp) == true) {
                    up_p.add();
                    up_p.last().get<0>() = p.getKey();
                    Particles.getProp<1>(p)[0] =  1;
                    Particles.getProp<1>(p)[1] =  0;
                }
                else if (down.isInside(xp) == true) {
                    if (xp[0]==0 && xp[1]==0) {
                        ref_p.add();
                        ref_p.last().get<0>() = p.getKey();
                    }
                    else {
                        dw_p.add();
                        dw_p.last().get<0>() = p.getKey();
                    }
                }
                else if (left.isInside(xp) == true) {
                    l_p.add();
                    l_p.last().get<0>() = p.getKey();
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                else if (right.isInside(xp) == true) {
                    r_p.add();
                    r_p.last().get<0>() = p.getKey();
                }
                else if (up_d.isInside(xp) == true) {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                }
                else if (down_u.isInside(xp) == true) {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();

                }
                else if (left_r.isInside(xp) == true) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                }
                else if (right_l.isInside(xp) == true) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                }
                else {
                    bulk.add();
                    bulk.last().get<0>() = p.getKey();
                    Particles.getProp<0>(p) = 0.0;
                }

                ++it2;
            }

             for(int j=0;j<up_p1.size();j++)
             {
                 auto p1=up_p1.get<0>(j);
                 Point<2, double> xp = Particles.getPos(p1);
                 Particles.getProp<3>(p1) =  xp[0];
             }
             for(int j=0;j<l_p1.size();j++)
             {
                 auto p1=l_p1.get<0>(j);
                 Point<2, double> xp = Particles.getPos(p1);
                 Particles.getProp<3>(p1) =  xp[1];
             }
             for(int j=0;j<r_p1.size();j++)
             {
                 auto p1=r_p1.get<0>(j);
                 Point<2, double> xp = Particles.getPos(p1);
                 Particles.getProp<3>(p1) =  xp[1];
             }
             for(int j=0;j<dw_p1.size();j++)
             {
                 auto p1=dw_p1.get<0>(j);
                 Point<2, double> xp = Particles.getPos(p1);
                 Particles.getProp<3>(p1)=  xp[0];
             }
            Particles.write_frame("Re1000-1e-4-Lid",0);

        for(int j=0;j<up_p1.size();j++)
        { if(j==0)
            {   auto p=up_p.get<0>(j);
                auto p1=up_p1.get<0>(j);
                auto p2=l_p.get<0>(l_p.size()-1);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
            auto p=up_p.get<0>(j+1);
            auto p1=up_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
            if(j==up_p1.size()-1)
            {   auto p=up_p.get<0>(j+2);
                auto p1=up_p1.get<0>(j);
                auto p2=r_p.get<0>(r_p.size()-1);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
        }
        for(int j=0;j<l_p1.size();j++)
        {
            auto p=l_p.get<0>(j+1);
            auto p1=l_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
        }
        for(int j=0;j<r_p1.size();j++)
        {
            auto p=r_p.get<0>(j+1);
            auto p1=r_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
        }
        for(int j=0;j<dw_p1.size();j++)
        {   auto p=dw_p.get<0>(j);
            auto p1=dw_p1.get<0>(j);
            Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
            if(j==0)
            {
                auto p = ref_p.get<0>(j);
                auto p2=l_p.get<0>(0);
                Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
            }
            if(j==dw_p1.size()-1)
            {   auto p=dw_p.get<0>(j+1);
                auto p2=r_p.get<0>(0);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
            }
        }
            Particles.write_frame("Re1000-1e-4-Lid",1);
            Derivative_x Dx(Particles, 2, rCut,2);
            Derivative_y Dy(Particles, 2, rCut,2);
            Gradient Grad(Particles, 2, rCut,2);
            Laplacian Lap(Particles, 3, rCut, 3);
            Advection Adv(Particles, 2, rCut, 2);
            Divergence Div(Particles, 2, rCut, 2);
            double dt=5e-4;
            int n=3;
            double nu=1e-2;
*//*        divV=Div(V);
        DCPSE_scheme<equations2,decltype(Particles)> Solver(2*rCut, Particles);
        auto Pressure_Poisson = Lap(P);
        auto D_y=Dy(P);
        auto D_x=Dx(P);
        Solver.impose(Pressure_Poisson,bulk,prop_id<5>());
        Solver.impose(D_y, up_p,0);
        Solver.impose(D_y, dw_p,0);
        Solver.impose(D_x, l_p,0);
        Solver.impose(D_x, r_p, 0);
        Solver.impose(P, ref_p, 0);
        Solver.solve(P);

        std::cout<<"Poisson Solved"<<std::endl;
        V=V-Grad(P);
        Particles.write_frame("Re1000-1e-4-Lid",0);
        std::cout<<"Init Done"<<std::endl;
        return;*//*
            for(int i=2; i<=n ;i++)
            {   dV=V+dt*(nu*Lap(V)-Adv(V,V));
                RHS=1/dt*Div(dV);
                std::cout<<"RHS Done"<<std::endl;

  *//*              //Neumann Copy

                for(int j=0;j<up_p1.size();j++)
                { if(j==0)
                    {   auto p=up_p.get<0>(j);
                        auto p1=up_p1.get<0>(j);
                        auto p2=l_p.get<0>(l_p.size()-1);
                        Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
                    }
                    auto p=up_p.get<0>(j+1);
                    auto p1=up_p1.get<0>(j);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                    if(j==up_p1.size()-1)
                    {   auto p=up_p.get<0>(j+2);
                        auto p1=up_p1.get<0>(j);
                        auto p2=r_p.get<0>(r_p.size()-1);
                        Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
                    }
                }
                for(int j=0;j<l_p1.size();j++)
                {
                    auto p=l_p.get<0>(j+1);
                    auto p1=l_p1.get<0>(j);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                }
                for(int j=0;j<r_p1.size();j++)
                {
                    auto p=r_p.get<0>(j+1);
                    auto p1=r_p1.get<0>(j);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                }
                for(int j=0;j<dw_p1.size();j++)
                {   auto p=dw_p.get<0>(j);
                    auto p1=dw_p1.get<0>(j);
                    Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                    if(j==0)
                    {
                        auto p = ref_p.get<0>(j);
                        auto p2=l_p.get<0>(0);
                        Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
                    }
                    if(j==dw_p1.size()-1)
                    {   auto p=dw_p.get<0>(j+1);
                        auto p2=r_p.get<0>(0);
                        Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
                    }
                }        for(int j=0;j<up_p1.size();j++)
                { if(j==0)
                    {   auto p=up_p.get<0>(j);
                        auto p1=up_p1.get<0>(j);
                        auto p2=l_p.get<0>(l_p.size()-1);
                        Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
                    }
                    auto p=up_p.get<0>(j+1);
                    auto p1=up_p1.get<0>(j);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                    if(j==up_p1.size()-1)
                    {   auto p=up_p.get<0>(j+2);
                        auto p1=up_p1.get<0>(j);
                        auto p2=r_p.get<0>(r_p.size()-1);
                        Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
                    }
                }
                for(int j=0;j<l_p1.size();j++)
                {
                    auto p=l_p.get<0>(j+1);
                    auto p1=l_p1.get<0>(j);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                }
                for(int j=0;j<r_p1.size();j++)
                {
                    auto p=r_p.get<0>(j+1);
                    auto p1=r_p1.get<0>(j);
                    Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                }
                for(int j=0;j<dw_p1.size();j++)
                {   auto p=dw_p.get<0>(j);
                    auto p1=dw_p1.get<0>(j);
                    Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                    if(j==0)
                    {
                        auto p = ref_p.get<0>(j);
                        auto p2=l_p.get<0>(0);
                        Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
                    }
                    if(j==dw_p1.size()-1)
                    {   auto p=dw_p.get<0>(j+1);
                        auto p2=r_p.get<0>(0);
                        Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                        Particles.getProp<3>(p2)=  Particles.getProp<3>(p1);
                    }
                }
*//*
                DCPSE_scheme<equations2,decltype(Particles)> Solver(2*rCut, Particles);
                auto Pressure_Poisson = Lap(P);
                auto D_y=Dy(P);
                auto D_x=Dx(P);
                Solver.impose(Pressure_Poisson,FullBulk,prop_id<3>());
*//*                Solver.impose(P, up_p,0);
                Solver.impose(P, dw_p,0);
                Solver.impose(P, l_p,0);
                Solver.impose(P, r_p, 0);*//*
                Solver.impose(P, ref_p, 0);
                Solver.solve(P);
                std::cout<<"Poisson Solved"<<std::endl;
                Vtemp=dV-dt*Grad(P);
                V=Vtemp;
                for(int j=0;j<up_p.size();j++)
                {   auto p=up_p.get<0>(j);
                    Particles.getProp<1>(p)[0] =  1;
                    Particles.getProp<1>(p)[1] =  0;
                }
                for(int j=0;j<l_p.size();j++)
                {   auto p=l_p.get<0>(j);
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                for(int j=0;j<r_p.size();j++)
                {   auto p=r_p.get<0>(j);
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                for(int j=0;j<dw_p.size();j++)
                {   auto p=dw_p.get<0>(j);
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                Particles.getProp<1>(0)[0] =  0;
                Particles.getProp<1>(0)[1] =  0;
                divV=Div(V);
                std::cout<<"Velocity updated"<<std::endl;
                Particles.write_frame("Re1000-1e-4-Lid",i);
                std::cout<<i<<std::endl;
            }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////












    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





*//*

    BOOST_AUTO_TEST_CASE(dcpse_Lid_RotProj) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        //                                  P        V                 V1                      V2           Dv              RHS    Vtemp                     Proj_lap
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,double>> Particles(0, box, bc, ghost);
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> ref_p2;
        openfpm::vector<aggregate<int>> FullBulk;


        // Here fill up the boxes for particle detection.

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
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (xp[0]!=0 || xp[1]!=0){
                FullBulk.add();
                FullBulk.last().get<0>() = p.getKey();
            }
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (down.isInside(xp) == true) {
                if (xp[0]==0 && xp[1]==0) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                }
                else {
                    dw_p.add();
                    dw_p.last().get<0>() = p.getKey();
                }
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) = 0.0;
            }

            ++it2;
        }
        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V1 = getV<2>(Particles);
        auto V2 = getV<3>(Particles);
        auto dV = getV<4>(Particles);
        auto RHS = getV<5>(Particles);
        auto Vtemp = getV<6>(Particles);
        auto divV = getV<7>(Particles);
        auto Proj_lap=getV<8>(Particles);

        Derivative_x Dx(Particles, 2, rCut,2);
        Derivative_y Dy(Particles, 2, rCut,2);
        Gradient Grad(Particles, 2, rCut,2);
        Laplacian Lap(Particles, 2, rCut, 2);
        Advection Adv(Particles, 2, rCut, 2);
        Divergence Div(Particles, 2, rCut, 2);
        double dt=1e-3;
        int n=3;
        double nu=1e-2;
        Particles.write_frame("Re1000-1e-4-Lid",0);
        for(int i=1; i<=n ;i++)
        {   dV=2*dt*(nu*Lap(V)-Adv(V,V)-Grad(P));

            RHS=1/dt*Div(dV);
            std::cout<<"RHS Done"<<std::endl;

            DCPSE_scheme<equations2,decltype(Particles)> Solver(2*rCut, Particles);
            auto Pressure_Poisson = Lap(P);
            auto D_y=Dy(P);
            auto D_x=Dx(P);
            Solver.impose(Pressure_Poisson,bulk,prop_id<3>());
            Solver.impose(D_y, up_p,prop_id<3>());
            Solver.impose(D_y, dw_p,prop_id<3>());
            Solver.impose(D_x, l_p,prop_id<3>());
            Solver.impose(D_x, r_p, prop_id<3>());
            Solver.impose(P, ref_p, 0);
            Solver.solve(P);
            std::cout<<"Poisson Solved"<<std::endl;
            Proj_lap=Dx(P)-RHS;
            Vtemp=dV-dt*Grad(P);
            V=Vtemp;
            //divV=Div(V);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  1;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<l_p.size();j++)
            {   auto p=l_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<r_p.size();j++)
            {   auto p=r_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            Particles.getProp<1>(0)[0] =  0;
            Particles.getProp<1>(0)[1] =  0;
            std::cout<<"Velocity updated"<<std::endl;
            Particles.write_frame("Re1000-1e-4-Lid",i);
            std::cout<<i<<std::endl;
        }
    }
*/





BOOST_AUTO_TEST_SUITE_END()


