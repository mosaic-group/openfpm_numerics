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
#include "DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "EqnsStructPetsc.hpp"
//#include "EqnsStructEigen.hpp"
const bool equations3d1::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations3d3::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations2d1::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations2d2::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations2d1p::boundary[] = {PERIODIC, NON_PERIODIC};
const bool equations2d2p::boundary[] = {PERIODIC, NON_PERIODIC};

//template<typename T>
//struct Debug;



BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)
    BOOST_AUTO_TEST_CASE(dcpse_Active2DPetsc) {
        const size_t sz[2] = {41, 41};
        Box<2, double> box({0, 0}, {10,10});
        double Lx=box.getHigh(0);
        double Ly=box.getHigh(1);
        size_t bc[2] = {PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(rCut);

/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div   H_t   */
        vector_dist<2, double, aggregate<VectorS<2, double>,VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>,VectorS<2, double>, VectorS<2, double> , VectorS<2, double>, double,double,double,double,double,double,double,VectorS<2, double>,double,double>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y*0.99999;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkP;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        constexpr int x      =     0;
        constexpr int y      =     1;

        constexpr int Polarization      =     0;
        constexpr int Velocity          =     1;
        constexpr int Vorticity         =     2;
        constexpr int ExtForce          =     3;
        constexpr int Pressure          =     4;
        constexpr int Strain_rate       =     5;
        constexpr int Stress            =     6;
        constexpr int MolField          =     7;

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dP = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);

        double eta       =     1.0;
        double nu        =     -0.5;
        double gama      =     0.1;
        double zeta      =     0.07;
        double Ks        =     1;
        double Kb        =     1;
        double lambda    =     0.1;
        double delmu     =     -1;
        g    =     0;

        // Here fill up the boxes for particle boundary detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        Box<2, double> mid({box.getHigh(0)/2.0 - spacing, box.getHigh(1)/2.0 - spacing / 2.0},
                           {box.getHigh(0)/2.0 , box.getHigh(1)/2.0 + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                //Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                //Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                //  Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                bulkP.add();
                bulkP.last().get<0>() = p.getKey();
                // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            } else if (right.isInside(xp) == true){
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                bulkP.add();
                bulkP.last().get<0>() = p.getKey();
                // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                //Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }
            else {
                if(mid.isInside(xp) == true) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                    // Particles.getProp<0>(p)[x]= cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                    // Particles.getProp<0>(p)[y] = sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                    Particles.getProp<4>(p) = 0;
                }
                else{
                    bulkP.add();
                    bulkP.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }


            ++it2;
        }



/*        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_xy Dxy(Particles, 2, rCut,1.9,support_options::RADIUS);
        auto Dyx=Dxy;
        Derivative_xx Dxx(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);*/

        int ord=2;
        double sampling_factor=1.9;
        Derivative_x Dx(Particles, ord, rCut,sampling_factor);
        Derivative_y Dy(Particles, ord, rCut,sampling_factor);
        Derivative_xy Dxy(Particles, ord, rCut,sampling_factor);
        auto Dyx=Dxy;
        Derivative_xx Dxx(Particles, ord, rCut,sampling_factor);
        Derivative_yy Dyy(Particles, ord, rCut,sampling_factor);
        Gradient Grad(Particles, ord, rCut,sampling_factor);
        Laplacian Lap(Particles, ord, rCut,sampling_factor);
        Advection Adv(Particles, ord, rCut,sampling_factor);
        Divergence Div(Particles, ord, rCut,sampling_factor);

        sigma[x][x] =    -Ks * Dx(Pol[x]) * Dx(Pol[x])- Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
        sigma[x][y] =    -Ks * Dy(Pol[y]) * Dx(Pol[y])- Kb * Dy(Pol[y]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[x]) * Dx(Pol[y]);
        sigma[y][x] =    -Ks * Dx(Pol[x]) * Dy(Pol[x])- Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
        sigma[y][y] =    -Ks * Dy(Pol[y]) * Dy(Pol[y])- Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);

        h[y] = Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) - Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y]));


        f1 =     gama*nu*Pol[x]*Pol[x]* (Pol[x]*Pol[x] - Pol[y]*Pol[y]) / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f2 = 2.0*gama*nu*Pol[x]*Pol[y]* (Pol[x]*Pol[x] - Pol[y]*Pol[y]) / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f3 =     gama*nu*Pol[y]*Pol[y]* (Pol[x]*Pol[x] - Pol[y]*Pol[y]) / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f4 = 2.0*gama*nu*Pol[x]*Pol[x]*Pol[x]*Pol[y] / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f5 = 4.0*gama*nu*Pol[x]*Pol[x]*Pol[y]*Pol[y] / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f6 = 2.0*gama*nu*Pol[x]*Pol[x]*Pol[x]*Pol[y] / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);

        dV[x] = -0.5*Dy(h[y]) + zeta*Dx(delmu*Pol[x]*Pol[x]) + zeta*Dy(delmu*Pol[x]*Pol[y]) - zeta*Dx(0.5*delmu*(Pol[x]*Pol[x]+Pol[y]*Pol[y])) - 0.5*nu*Dx(-2*h[y]*Pol[x]*Pol[y])
                - 0.5*nu*Dy(h[y]*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) - g[x]
                - 0.5*nu*Dx(-gama*lambda*delmu*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - 0.5*Dy(-2*gama*lambda*delmu*(Pol[x]*Pol[y]));


        dV[y] = -0.5*Dx(-h[y]) + zeta*Dy(delmu*Pol[y]*Pol[y]) + zeta*Dx(delmu*Pol[x]*Pol[y]) - zeta*Dy(0.5*delmu*(Pol[x]*Pol[x]+Pol[y]*Pol[y])) - 0.5*nu*Dy(-2*h[y]*Pol[x]*Pol[y])
                - 0.5*nu*Dx(h[y]*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) - g[y]
                - 0.5*nu*Dy(gama*lambda*delmu*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - 0.5*Dx(-2*gama*lambda*delmu*(Pol[x]*Pol[y]));

        Particles.write_frame("Polar",0);
        //Velocity Solution n iterations
        eq_id vx,vy;
        vx.setId(0);
        vy.setId(1);
        double sum=0,sum1=0;
        int n=30;
/*        auto Stokes1 = nu*Lap(V[x]) + 0.5*nu*(f1*Dxx(V[x])+Dx(f1)*Dx(V[x])) + (Dx(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxx(V[y])+Dyx(V[x]))) + (Dx(f3)*Dy(V[y])+ f3*Dyx(V[y])) + (Dy(f4)*Dx(V[x])+f4*Dxy(V[x])) + (Dy(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxy(V[y])+Dyy(V[x]))) + ( Dy(f6)*Dy(V[y])+f6*Dyy(V[y]) );
        auto Stokes2 = nu*Lap(V[y]) + 0.5*nu*(f1*Dxy(V[x])+Dy(f1)*Dx(V[x])) + (Dy(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxy(V[y])+Dyy(V[x]))) + (Dy(f3)*Dy(V[y])+ f3*Dyy(V[y])) + (Dx(f4)*Dx(V[x])+f4*Dxx(V[x])) + (Dx(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxx(V[y])+Dyx(V[x]))) + ( Dx(f6)*Dy(V[y])+f6*Dyx(V[y]) );
        auto Helmholtz = Lap(H);
        auto D_y=Dy(H);
        auto D_x=Dx(H);*/
        petsc_solver<double> solverPetsc;
        solverPetsc.setRestart(500);
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setPreconditioner(PCKSP);


        petsc_solver<double> solverPetsc2;
        solverPetsc2.setRestart(500);
        solverPetsc2.setSolver(KSPGMRES);
        //solverPetsc2.setPreconditioner(PCKSP);
        /*
              solverPetsc.setPreconditioner(PCHYPRE_BOOMERAMG);
              solverPetsc2.setPreconditioner(PCHYPRE_BOOMERAMG);

              solver.setPreconditionerAMG_nl(6);
              solver.setPreconditionerAMG_maxit(10);
              solver.setPreconditionerAMG_relax("SOR/Jacobi");
              solver.setPreconditionerAMG_cycleType("V",0,4);
              solver.setPreconditionerAMG_coarsenNodalType(0);
              solver.setPreconditionerAMG_coarsen("HMIS");*/

        for(int i=1; i<=n ;i++)
        {   RHS[x]=Dx(P)+dV[x];
            RHS[y]=Dy(P)+dV[y];
            Particles.ghost_get<10>();
            DCPSE_scheme<equations2d2p,decltype(Particles)> Solver( Particles);
            auto Stokes1 = nu*Lap(V[x]) + 0.5*nu*(f1*Dxx(V[x]) + Dx(f1)*Dx(V[x])) + (Dx(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxx(V[y])+Dyx(V[x]))) + (Dx(f3)*Dy(V[y])+ f3*Dyx(V[y])) + (Dy(f4)*Dx(V[x])+f4*Dxy(V[x])) + (Dy(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxy(V[y])+Dyy(V[x]))) + ( Dy(f6)*Dy(V[y])+f6*Dyy(V[y]) );
            auto Stokes2 = nu*Lap(V[y]) + 0.5*nu*(f1*Dxy(V[x]) + Dy(f1)*Dx(V[x])) + (Dy(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxy(V[y])+Dyy(V[x]))) + (Dy(f3)*Dy(V[y])+ f3*Dyy(V[y])) + (Dx(f4)*Dx(V[x])+f4*Dxx(V[x])) + (Dx(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxx(V[y])+Dyx(V[x]))) + ( Dx(f6)*Dy(V[y])+f6*Dyx(V[y]) );
            Solver.impose(Stokes1,bulk,RHS[0],vx);
            Solver.impose(Stokes2,bulk,RHS[1],vy);
            Solver.impose(V[x], up_p,0,vx);
            Solver.impose(V[y], up_p,0,vy);
            Solver.impose(V[x], dw_p,0,vx);
            Solver.impose(V[y], dw_p,0,vy);

            //Solver.try_solve_with_solver(solverPetsc,V[x],V[y]);
            //Solver.solve_with_solver(solverPetsc,V[x],V[y]);
            return;
            std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<Velocity>();
            div=-Div(V);
            Particles.ghost_get<19>();
            DCPSE_scheme<equations2d1p,decltype(Particles)> SolverH(Particles,options_solver::LAGRANGE_MULTIPLIER);//,
            auto Helmholtz = Lap(H);
            auto D_y=Dy(H);
            SolverH.impose(Helmholtz,bulk,prop_id<19>());
            SolverH.impose(D_y, up_p,0);
            SolverH.impose(D_y, dw_p,0);

            //SolverH.try_solve_with_solver(solverPetsc2,H);
            //SolverH.solve_with_solver(solverPetsc2,H);

            Particles.ghost_get<17>();
            Particles.ghost_get<Velocity>();
            //std::cout << "Helmholtz Solved" << std::endl;
            V=V+Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            P=P+Lap(H);
            Particles.ghost_get<Velocity>();
            Particles.ghost_get<Pressure>();
            std::cout << "V,P Corrected" << std::endl;
            sum=0;
            sum1=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<18>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<18>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<18>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<18>(p)[1]- Particles.getProp<1>(p)[1]);
                sum1+=Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1];
            }
            sum=sqrt(sum);
            sum1=sqrt(sum1);
            V_t=V;
            std::cout << "Rel l2 cgs err in V at "<<i<<"= " <<sum/sum1<< std::endl;
            Particles.write_frame("Polar",i);

        }
        Particles.deleteGhost();
        Particles.write_frame("Polar",n+1);


        //u[x][x] =   Dx(V[x]);
        //u[x][y] =   0.5*(Dx(V[y])+Dy(V[x]));
        //u[y][x] =   0.5*(Dy(V[x])+Dx(V[y]));
        //u[y][y] =   Dy(V[y]);


/*
              double dt=5e-4;
              int n=5;
              double Re=1/3e-3;
              std::cout<<"Init Done"<<std::endl;
              for(int i =1; i<=n ;i++)
              {
                  dV=(1/Re*Lap(V)-Adv(V,V));
                  std::cout<<"dV Done"<<std::endl;
                  RHS = Div(dV);
                  std::cout<<"RHS Done"<<std::endl;
                  DCPSE_scheme<equations,decltype(Particles)> Solver( Particles);
                  auto Pressure_Poisson = Lap(H);
                  auto D_x = Dx(H);
                  auto D_y = Dy(H);
                  Solver.impose(Pressure_Poisson, bulk, prop_id<3>());
                  Solver.impose(D_y, up_p, 0);
                  Solver.impose(-D_y, dw_p,0);
                  Solver.impose(-D_x, l_p,0);
                  Solver.impose(D_x, r_p, 0);
                  Solver.impose(H, ref_p, 0);
                  Solver.solve(P);
                  std::cout<<"Poisson Solved"<<std::endl;
                  V=dt*(dV-Grad(P));
                  std::cout<<"Velocity updated"<<std::endl;

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
                  std::cout<<"Boundary done"<<std::endl;
                  //           if (i%10==0)
                  Particles.write_frame("Lid",i);
                  std::cout<<i<<std::endl;
              */
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_Active2D) {
        const size_t sz[2] = {31, 31};
        Box<2, double> box({0, 0}, {10,10});
        double Lx=box.getHigh(0);
        double Ly=box.getHigh(1);
        size_t bc[2] = {PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(rCut);

/*                                          pol                             V         vort                 Ext    Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div   H_t   */
        vector_dist<2, double, aggregate<VectorS<2, double>,VectorS<2, double>, double[2][2], VectorS<2, double>, double, double[2][2], double[2][2], VectorS<2, double>,VectorS<2, double>, VectorS<2, double> , VectorS<2, double>, double,double,double,double,double,double,double,VectorS<2, double>,double,double>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y*0.99999;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulkP;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        constexpr int x      =     0;
        constexpr int y      =     1;

        constexpr int Polarization      =     0;
        constexpr int Velocity          =     1;
        constexpr int Vorticity         =     2;
        constexpr int ExtForce          =     3;
        constexpr int Pressure          =     4;
        constexpr int Strain_rate       =     5;
        constexpr int Stress            =     6;
        constexpr int MolField          =     7;

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dP = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);

        double eta       =     1.0;
        double nu        =     -0.5;
        double gama      =     0.1;
        double zeta      =     0.07;
        double Ks        =     1;
        double Kb        =     1;
        double lambda    =     0.1;
        double delmu     =     -1;
        g    =     0;

        // Here fill up the boxes for particle boundary detection.

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});
        Box<2, double> mid({box.getHigh(0)/2.0 - spacing, box.getHigh(1)/2.0 - spacing / 2.0},
                             {box.getHigh(0)/2.0 , box.getHigh(1)/2.0 + spacing / 2.0});

        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(mid);

        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                //Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                //Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }
            else if (down.isInside(xp) == true) {
                    dw_p.add();
                    dw_p.last().get<0>() = p.getKey();
                  //  Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                   // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                bulkP.add();
                bulkP.last().get<0>() = p.getKey();
               // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
               // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            } else if (right.isInside(xp) == true){
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                bulkP.add();
                bulkP.last().get<0>() = p.getKey();
               // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                //Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }
            else {
                if(mid.isInside(xp) == true) {
                    ref_p.add();
                    ref_p.last().get<0>() = p.getKey();
                   // Particles.getProp<0>(p)[x]= cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                   // Particles.getProp<0>(p)[y] = sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
                    Particles.getProp<4>(p) = 0;
                }
                else{
                    bulkP.add();
                    bulkP.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
               // Particles.getProp<0>(p)[x]=  cos(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
               // Particles.getProp<0>(p)[y] =  sin(2 * M_PI * (cos((2 * xp[x]- sz[x]) / sz[x]) - sin((2 * xp[y] - sz[y]) / sz[y])));
            }


            ++it2;
        }


        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_xy Dxy(Particles, 2, rCut,1.9,support_options::RADIUS);
        auto Dyx=Dxy;
        Derivative_xx Dxx(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);

        sigma[x][x] =    -Ks * Dx(Pol[x]) * Dx(Pol[x])- Kb * Dx(Pol[y]) * Dx(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);
        sigma[x][y] =    -Ks * Dy(Pol[y]) * Dx(Pol[y])- Kb * Dy(Pol[y]) * Dx(Pol[x]) + (Kb - Ks) * Dx(Pol[x]) * Dx(Pol[y]);
        sigma[y][x] =    -Ks * Dx(Pol[x]) * Dy(Pol[x])- Kb * Dx(Pol[y]) * Dy(Pol[y]) + (Kb - Ks) * Dy(Pol[x]) * Dy(Pol[y]);
        sigma[y][y] =    -Ks * Dy(Pol[y]) * Dy(Pol[y])- Kb * Dy(Pol[x]) * Dy(Pol[x]) + (Kb - Ks) * Dy(Pol[x]) * Dx(Pol[y]);

        h[y] = Pol[x] * (Ks * Dyy(Pol[y]) + Kb * Dxx(Pol[y]) + (Ks - Kb) * Dxy(Pol[x])) - Pol[y] * (Ks * Dxx(Pol[x]) + Kb * Dyy(Pol[x]) + (Ks - Kb) * Dxy(Pol[y]));


        f1 =     gama*nu*Pol[x]*Pol[x]* (Pol[x]*Pol[x] - Pol[y]*Pol[y]) / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f2 = 2.0*gama*nu*Pol[x]*Pol[y]* (Pol[x]*Pol[x] - Pol[y]*Pol[y]) / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f3 =     gama*nu*Pol[y]*Pol[y]* (Pol[x]*Pol[x] - Pol[y]*Pol[y]) / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f4 = 2.0*gama*nu*Pol[x]*Pol[x]*Pol[x]*Pol[y] / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f5 = 4.0*gama*nu*Pol[x]*Pol[x]*Pol[y]*Pol[y] / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);
        f6 = 2.0*gama*nu*Pol[x]*Pol[x]*Pol[x]*Pol[y] / (Pol[x]*Pol[x] + Pol[y]*Pol[y]);

        dV[x] = -0.5*Dy(h[y]) + zeta*Dx(delmu*Pol[x]*Pol[x]) + zeta*Dy(delmu*Pol[x]*Pol[y]) - zeta*Dx(0.5*delmu*(Pol[x]*Pol[x]+Pol[y]*Pol[y])) - 0.5*nu*Dx(-2*h[y]*Pol[x]*Pol[y])
                  - 0.5*nu*Dy(h[y]*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - Dx(sigma[x][x]) - Dy(sigma[x][y]) - g[x]
                  - 0.5*nu*Dx(-gama*lambda*delmu*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - 0.5*Dy(-2*gama*lambda*delmu*(Pol[x]*Pol[y]));


        dV[y] = -0.5*Dx(-h[y]) + zeta*Dy(delmu*Pol[y]*Pol[y]) + zeta*Dx(delmu*Pol[x]*Pol[y]) - zeta*Dy(0.5*delmu*(Pol[x]*Pol[x]+Pol[y]*Pol[y])) - 0.5*nu*Dy(-2*h[y]*Pol[x]*Pol[y])
                  - 0.5*nu*Dx(h[y]*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - Dx(sigma[y][x]) - Dy(sigma[y][y]) - g[y]
                  - 0.5*nu*Dy(gama*lambda*delmu*(Pol[x]*Pol[x] - Pol[y]*Pol[y])) - 0.5*Dx(-2*gama*lambda*delmu*(Pol[x]*Pol[y]));

        Particles.write_frame("Polar",0);
        //Velocity Solution n iterations
        eq_id vx,vy;
        vx.setId(0);
        vy.setId(1);
        double sum=0,sum1=0;
        int n=30;
/*        auto Stokes1 = nu*Lap(V[x]) + 0.5*nu*(f1*Dxx(V[x])+Dx(f1)*Dx(V[x])) + (Dx(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxx(V[y])+Dyx(V[x]))) + (Dx(f3)*Dy(V[y])+ f3*Dyx(V[y])) + (Dy(f4)*Dx(V[x])+f4*Dxy(V[x])) + (Dy(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxy(V[y])+Dyy(V[x]))) + ( Dy(f6)*Dy(V[y])+f6*Dyy(V[y]) );
        auto Stokes2 = nu*Lap(V[y]) + 0.5*nu*(f1*Dxy(V[x])+Dy(f1)*Dx(V[x])) + (Dy(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxy(V[y])+Dyy(V[x]))) + (Dy(f3)*Dy(V[y])+ f3*Dyy(V[y])) + (Dx(f4)*Dx(V[x])+f4*Dxx(V[x])) + (Dx(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxx(V[y])+Dyx(V[x]))) + ( Dx(f6)*Dy(V[y])+f6*Dyx(V[y]) );
        auto Helmholtz = Lap(H);
        auto D_y=Dy(H);
        auto D_x=Dx(H);*/

        for(int i=1; i<=n ;i++)
        {   RHS[x]=Dx(P)+dV[x];
            RHS[y]=Dy(P)+dV[y];
            Particles.ghost_get<10>();
            DCPSE_scheme<equations2d2p,decltype(Particles)> Solver( Particles);
            auto Stokes1 = nu*Lap(V[x]) + 0.5*nu*(f1*Dxx(V[x])+Dx(f1)*Dx(V[x])) + (Dx(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxx(V[y])+Dyx(V[x]))) + (Dx(f3)*Dy(V[y])+ f3*Dyx(V[y])) + (Dy(f4)*Dx(V[x])+f4*Dxy(V[x])) + (Dy(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxy(V[y])+Dyy(V[x]))) + ( Dy(f6)*Dy(V[y])+f6*Dyy(V[y]) );
            auto Stokes2 = nu*Lap(V[y]) + 0.5*nu*(f1*Dxy(V[x])+Dy(f1)*Dx(V[x])) + (Dy(f2)*0.5*(Dx(V[y])+Dy(V[x])) + f2*0.5*(Dxy(V[y])+Dyy(V[x]))) + (Dy(f3)*Dy(V[y])+ f3*Dyy(V[y])) + (Dx(f4)*Dx(V[x])+f4*Dxx(V[x])) + (Dx(f5)*0.5*(Dx(V[y])+Dy(V[x]))+f5*0.5*(Dxx(V[y])+Dyx(V[x]))) + ( Dx(f6)*Dy(V[y])+f6*Dyx(V[y]) );
            Solver.impose(Stokes1,bulk,RHS[0],vx);
            Solver.impose(Stokes2,bulk,RHS[1],vy);
            Solver.impose(V[x], up_p,0,vx);
            Solver.impose(V[y], up_p,0,vy);
            Solver.impose(V[x], dw_p,0,vx);
            Solver.impose(V[y], dw_p,0,vy);
            Solver.solve(V[x],V[y]);
            std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<Velocity>();
            div=Div(V);
            Particles.ghost_get<19>();
            DCPSE_scheme<equations2d1p,decltype(Particles)> SolverH(Particles,options_solver::LAGRANGE_MULTIPLIER);//,
            auto Helmholtz = Lap(H);
            auto D_y=Dy(H);
            SolverH.impose(Helmholtz,bulk,prop_id<19>());
            SolverH.impose(D_y, up_p,0);
            SolverH.impose(D_y, dw_p,0);
            SolverH.solve(H);
            Particles.ghost_get<17>();
            Particles.ghost_get<Velocity>();
            //std::cout << "Helmholtz Solved" << std::endl;
            V=V-Grad(H);
            for(int j=0;j<up_p.size();j++)
            {   auto p=up_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            for(int j=0;j<dw_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            P=P-Lap(H);
            Particles.ghost_get<Velocity>();
            Particles.ghost_get<Pressure>();
            std::cout << "V,P Corrected" << std::endl;
            sum=0;
            sum1=0;
            for(int j=0;j<bulk.size();j++)
            {   auto p=bulk.get<0>(j);
                sum+=(Particles.getProp<18>(p)[0]-Particles.getProp<1>(p)[0])*(Particles.getProp<18>(p)[0]- Particles.getProp<1>(p)[0])+(Particles.getProp<18>(p)[1]- Particles.getProp<1>(p)[1])*(Particles.getProp<18>(p)[1]- Particles.getProp<1>(p)[1]);
                sum1+=Particles.getProp<1>(p)[0]*Particles.getProp<1>(p)[0]+Particles.getProp<1>(p)[1]*Particles.getProp<1>(p)[1];
            }
            sum=sqrt(sum);
            sum1=sqrt(sum1);
            V_t=V;
            std::cout << "Rel l2 cgs err in V at "<<i<<"= " <<sum/sum1<< std::endl;
            Particles.write_frame("Polar",i);

        }
        Particles.deleteGhost();
        Particles.write_frame("Polar",n+1);


        //u[x][x] =   Dx(V[x]);
        //u[x][y] =   0.5*(Dx(V[y])+Dy(V[x]));
        //u[y][x] =   0.5*(Dy(V[x])+Dx(V[y]));
        //u[y][y] =   Dy(V[y]);


/*
              double dt=5e-4;
              int n=5;
              double Re=1/3e-3;
              std::cout<<"Init Done"<<std::endl;
              for(int i =1; i<=n ;i++)
              {
                  dV=(1/Re*Lap(V)-Adv(V,V));
                  std::cout<<"dV Done"<<std::endl;
                  RHS = Div(dV);
                  std::cout<<"RHS Done"<<std::endl;
                  DCPSE_scheme<equations,decltype(Particles)> Solver( Particles);
                  auto Pressure_Poisson = Lap(H);
                  auto D_x = Dx(H);
                  auto D_y = Dy(H);
                  Solver.impose(Pressure_Poisson, bulk, prop_id<3>());
                  Solver.impose(D_y, up_p, 0);
                  Solver.impose(-D_y, dw_p,0);
                  Solver.impose(-D_x, l_p,0);
                  Solver.impose(D_x, r_p, 0);
                  Solver.impose(H, ref_p, 0);
                  Solver.solve(P);
                  std::cout<<"Poisson Solved"<<std::endl;
                  V=dt*(dV-Grad(P));
                  std::cout<<"Velocity updated"<<std::endl;

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
                  std::cout<<"Boundary done"<<std::endl;
                  //           if (i%10==0)
                  Particles.write_frame("Lid",i);
                  std::cout<<i<<std::endl;
              */
        }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_Lid) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        //                                  P        V                 Dv              RHS    Vtemp
        vector_dist<2, double, aggregate<double,VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,double>> Particles(0, box, bc, ghost);
        vector_dist<2, double, aggregate<double,VectorS<2, double>,double> > Particles_subset(0, box, bc, ghost);

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
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> ref_bulk;
        openfpm::vector<aggregate<int>> ref_bulk2;
        openfpm::vector<aggregate<int>> bulk_1;


        auto P = getV<0>(Particles);
        auto P_bulk = getV<2>(Particles_subset);
        auto V = getV<1>(Particles);
        auto V_bulk = getV<1>(Particles_subset);
        auto divV_bulk = getV<0>(Particles_subset);
        auto dV = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Vtemp = getV<4>(Particles);
        auto k1 = getV<5>(Particles);
        auto k2 = getV<6>(Particles);
        auto k3 = getV<7>(Particles);
        auto k4 = getV<8>(Particles);
        auto divV = getV<9>(Particles);


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
        vtk_box.write("Re1000-1e-3-vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
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
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                else{
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                }
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
            }
            else {
                 bulk.add();
                 bulk.last().get<0>() = p.getKey();
                 Particles.getProp<0>(p) = 0.0;
            }

            ++it2;
        }


        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
            Particles_subset.getLastProp<1>() = Particles.template getProp<1>(bulk.template get<0>(i));
        }



        Particles_subset.map();
        Particles_subset.ghost_get<0>();

        auto it3 = Particles_subset.getDomainIterator();

        while (it3.isNext()) {
            auto p = it3.get();
            Point<2, double> xp = Particles.getPos(p);

            if (xp[0]==0.5 && xp[1]==0.5) {
                ref_bulk.add();
                ref_bulk.last().get<0>() = p.getKey();
            }
            else if(xp[0]==0.5+spacing && xp[1]==0.5+spacing)
            {
                ref_bulk2.add();
                ref_bulk2.last().get<0>() = p.getKey();
            }
            else
            {
                bulk_1.add();
                bulk_1.last().get<0>() = p.getKey();
            }



            ++it3;
        }

        Derivative_x Dx(Particles, 2, rCut,2);
        Derivative_y Dy(Particles, 2, rCut,2);
        Gradient Grad_sub(Particles_subset, 2, rCut, 2),Grad(Particles, 2, rCut,2);
        Laplacian Lap_sub(Particles_subset, 2, rCut, 2),Lap(Particles, 2, rCut, 2);
        Advection Adv(Particles, 2, rCut, 2);
        Divergence Div(Particles, 2, rCut, 2),Div_bulk(Particles_subset, 2, rCut, 2);


        double dt=0.005;
        int n=50;
        double nu=1e-3;
        divV_bulk=Div_bulk(V_bulk);
        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<3>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            Particles.template getProp<9>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
        }
        //subset_create<0,1,2,4>(Particles,Particles_subset,bulk);
        DCPSE_scheme<equations2d1,decltype(Particles_subset)> Solver(Particles_subset);
        auto Pressure_Poisson = Lap_sub(P_bulk);
        Solver.impose(Pressure_Poisson, bulk_1,prop_id<0>());
        Solver.impose(P_bulk, ref_bulk, 1);
        Solver.impose(P_bulk, ref_bulk2, 1);
//        Solver.impose(P_bulk, ref_bulk2, 1);

//       Solver.impose(-D_y, dw_p,0);
//       Solver.impose(-D_x, l_p,0);
//       Solver.impose(D_x, r_p, 0);
//       Solver.impose(P, ref_p, 0);
        Solver.solve(P_bulk);

        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<0>(bulk.template get<0>(i)) = Particles_subset.getProp<2>(i);
        }

        std::cout<<"Poisson Solved"<<std::endl;
        V_bulk=V_bulk-Grad_sub(P_bulk);
        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<1>(bulk.template get<0>(i)) = Particles_subset.getProp<1>(i);
        }

        divV_bulk=Div_bulk(V_bulk);
        for (int i = 0 ; i < bulk.size() ; i++) {
            Particles.template getProp<9>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
        }

        Particles.write_frame("Re1000-1e-4-Lid",0);
        std::cout<<"Init Done"<<std::endl;
        return;
        for(int i =1; i<=n ;i++)
        {   dV=V+dt*(nu*Lap(V)-Adv(V,V));
            std::cout<<"dV Done"<<std::endl;
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_subset.getProp<1>(i) = Particles.template getProp<2>(bulk.template get<0>(i));
            }
            divV_bulk=1.0/dt*Div_bulk(V_bulk);
            for (int i = 0 ; i < bulk.size() ; i++)
            {
                Particles.template getProp<3>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            }
            //RHS=1.0/dt*Div(dV);
            std::cout<<"RHS Done"<<std::endl;
            DCPSE_scheme<equations2d1,decltype(Particles)> Solver(Particles);
            auto Pressure_Poisson = Lap(P);
            auto D_x = Dx(P);
            auto D_y = Dy(P);
            Solver.impose(Pressure_Poisson, bulk, prop_id<3>());
            Solver.impose(D_y, up_p, 0);
            Solver.impose(-D_y, dw_p,0);
            Solver.impose(-D_x, l_p,0);
            Solver.impose(D_x, r_p, 0);
            Solver.impose(P, ref_p, 0);
            Solver.solve(P);
            std::cout<<"Poisson Solved"<<std::endl;
/*          k1=dt*(dV-Grad(P));
            Vtemp=V+k1*0.5;
            k2=dt*(nu*Lap(Vtemp)-Adv(Vtemp,Vtemp)-Grad(P));
            Vtemp=V+k2*0.5;
            k3=dt*(nu*Lap(Vtemp)-Adv(Vtemp,Vtemp)-Grad(P));
            Vtemp=V+k3;
            k4=dt*(nu*Lap(Vtemp)-Adv(Vtemp,Vtemp)-Grad(P));
            Vtemp=V+0.16666666666*(k1+2*k2+2*k3+k4);*/
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
            for (int i = 0 ; i < bulk.size() ; i++) {
                Particles_subset.getProp<1>(i) = Particles.template getProp<1>(bulk.template get<0>(i));
            }
            divV_bulk=Div_bulk(V_bulk);
            for (int i = 0 ; i < bulk.size() ; i++)
            {
                Particles.template getProp<9>(bulk.template get<0>(i)) = Particles_subset.getProp<0>(i);
            }
            //divV=Div(V);
            std::cout<<"Velocity updated"<<std::endl;
                Particles.write_frame("Re1000-1e-4-Lid",i);
            std::cout<<i<<std::endl;
        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_Lid_sf) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3.1);
        double rCut = 3.1*spacing;
        std::cout<<spacing<<std::endl;
        //                                  sf    W   DW        RHS    Wnew  V
        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> Particles(0, box, bc, ghost);

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

        openfpm::vector<aggregate<int>> BULK;



        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) -  spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) +  spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                            {box.getHigh(0) -  spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> down({box.getLow(0) +  spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) -  spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) +   spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) -  spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) -  spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) +  spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) + spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                              {box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) -  spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) +  spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                               {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});
        Box<2, double> Bulk_box({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                                {box.getHigh(0) + spacing / 2.0, box.getHigh(1)  +spacing / 2.0});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(up_d);
        boxes.add(down);
        boxes.add(down_u);
        boxes.add(left);
        boxes.add(left_r);
        boxes.add(right);
        boxes.add(right_l);
        boxes.add(Bulk_box);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("box_sf.vtk");
       // Particles.write_frame("Re1000-1e-3-Lid_sf",0);

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;//-M_PI*M_PI*sin(M_PI*xp[0])*sin(M_PI*xp[1]);
            Particles.getProp<1>(p) =0;//sin(M_PI*xp[0])*sin(M_PI*xp[1]);
            Particles.getProp<2>(p) =0.0;
            Particles.getProp<3>(p) =0.0;
            Particles.getProp<4>(p) =0.0;
            Particles.getProp<5>(p)[0] =0.0;
            Particles.getProp<5>(p)[1] =0.0;

            BULK.add();
            BULK.last().get<0>() = p.getKey();

            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                Particles. template getProp<1>(p) = 12.0/spacing;// -2.0*Particles.getProp<0>(p)/(spacing*spacing) - 12.0/spacing;
            }
            else if (down.isInside(xp) == true) {
                    dw_p.add();
                    dw_p.last().get<0>() = p.getKey();
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            }
            else if (Bulk_box.isInside(xp) == true)
            {
                if (up_d.isInside(xp) == true) {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                }
                else if (down_u.isInside(xp) == true) {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                }else if (left_r.isInside(xp) == true) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                } else if (right_l.isInside(xp) == true) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }
/*        //test for copy
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
        //Copy Module
        for(int j=0;j<up_p1.size();j++)
        {   auto p=up_p.get<0>(j);
            auto p1=up_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
            if(j==0)
            {   auto p=l_p.get<0>(l_p.size()-1);
                auto p2=l_p.get<0>(l_p.size()-2);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
            if(j==up_p1.size()-1)
            {   auto p=r_p.get<0>(r_p.size()-1);
                auto p2=r_p.get<0>(r_p.size()-2);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
        }

        for(int j=0;j<l_p1.size();j++)
        {
            auto p=l_p.get<0>(j+2);
            auto p1=l_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
        }
        for(int j=0;j<r_p1.size();j++)
        {
            auto p=r_p.get<0>(j+2);
            auto p1=r_p1.get<0>(j);
            Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
        }
        for(int j=0;j<dw_p1.size();j++)
        {   auto p=dw_p.get<0>(j);
            auto p1=dw_p1.get<0>(j);
            Particles.getProp<3>(p)=  Particles.getProp<3>(p1);
            if(j==0)
            {   auto p=l_p.get<0>(0);
                auto p2=l_p.get<0>(1);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
            if(j==dw_p1.size()-1)
            {   auto p=r_p.get<0>(0);
                auto p2=r_p.get<0>(1);
                Particles.getProp<3>(p) =  Particles.getProp<3>(p1);
                Particles.getProp<3>(p2) =  Particles.getProp<3>(p1);
            }
        }*/

/*
        auto Sf = getV<0>(Particles);
        auto W = getV<1>(Particles);
        auto dW = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Wnew = getV<4>(Particles);
        auto V = getV<5>(Particles);
*/

        for(int j=0;j<up_p.size();j++) {
            auto p = up_p.get<0>(j);
            Particles.getProp<1>(p) =  -12.0/spacing;
        }

//        Particles.write_frame("Re1000-1e-3-Lid_sf",0);
//        Derivative_x Dx(Particles, 2, rCut,1);
//        Derivative_y Dy(Particles, 2, rCut,1);
//        Gradient Grad(Particles, 2, rCut,1);
//        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Laplacian Lap(Particles, 2, rCut,1.9);
//        Curl2D Curl(Particles, 2, rCut, 1);

        auto its = Particles.getDomainIterator();
        int ctr=0;
        while (its.isNext()) {
            auto p = its.get();
            Lap.DrawKernel<3>(Particles, p.getKey());
            Particles.write_frame("LapKer",ctr);
            for(int j=0;j<BULK.size();j++) {
                auto p1 = BULK.get<0>(j);
                Particles.getProp<3>(p1) =  0;
            }
            ++its;
            ctr++;
        }

/*        double dt=0.003;
        double nu=0.01;
        int n=5;
        std::cout<<"Init Done"<<std::endl;
        for(int i =0; i<=n ;i++)
        {   dW=Dx(Sf)*Dy(W)-Dy(Sf)*Dx(W)+nu*Lap(W);
            //Lap.DrawKernel<3>(Particles,837);
            Wnew=W+dt*dW;
            W=Wnew;
            //Copy Module
            for(int j=0;j<up_p1.size();j++)
            {   auto p=up_p.get<0>(j);
                auto p1=up_p1.get<0>(j);
                Particles.getProp<1>(p) =  -2.0*Particles.getProp<0>(p1)/(spacing*spacing)-12.0/spacing;
                if(j==0)
                {   auto p=l_p.get<0>(l_p.size()-1);
                    auto p2=l_p.get<0>(l_p.size()-2);
                    //Particles.getProp<1>(p) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                    Particles.getProp<1>(p2) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                }
                if(j==up_p1.size()-1)
                {   auto p=r_p.get<0>(r_p.size()-1);
                    auto p2=r_p.get<0>(r_p.size()-2);
                    //Particles.getProp<1>(p) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                    Particles.getProp<1>(p2) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                }
            }
            for(int j=0;j<l_p1.size();j++)
            {
                auto p=l_p.get<0>(j+2);
                auto p1=l_p1.get<0>(j);
                Particles.getProp<1>(p) =  -2.0*Particles.getProp<0>(p1)/(spacing*spacing);
            }
            for(int j=0;j<r_p1.size();j++)
            {
                auto p=r_p.get<0>(j+2);
                auto p1=r_p1.get<0>(j);
                Particles.getProp<1>(p) =  -2.0*Particles.getProp<0>(p1)/(spacing*spacing);
            }
            for(int j=0;j<dw_p1.size();j++)
            {   auto p=dw_p.get<0>(j);
                auto p1=dw_p1.get<0>(j);
                Particles.getProp<1>(p) =  -2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                if(j==0)
                {   auto p=l_p.get<0>(0);
                    auto p2=l_p.get<0>(1);
                    //Particles.getProp<1>(p) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                    Particles.getProp<1>(p2) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                }
                if(j==dw_p1.size()-1)
                {   auto p=r_p.get<0>(0);
                    auto p2=r_p.get<0>(1);
                   // Particles.getProp<1>(p) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                    Particles.getProp<1>(p2) =  0;//-2.0*Particles.getProp<0>(p1)/(spacing*spacing);
                }
            }
            std::cout<<"W Done"<<std::endl;

            DCPSE_scheme<equations,decltype(Particles)> Solver(2*rCut, Particles);
            auto Sf_Poisson = -Lap(Sf);
            Solver.impose(Sf_Poisson, bulk, prop_id<1>());
            Solver.impose(Sf, up_p, 0);
            Solver.impose(Sf, dw_p,0);
            Solver.impose(Sf, l_p,0);
            Solver.impose(Sf, r_p, 0);
            Solver.solve(Sf);
            V=Curl(Sf);
            std::cout<<"Poisson Solved"<<std::endl;
            Particles.write_frame("Re1000-1e-3-Lid_sf",i);
            //if (i%10==0)
            std::cout<<i<<std::endl;
        }*/
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_poisson_Robin_anal) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {50,50};
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

        // Add multi res patch 1

        {
        const size_t sz2[2] = {40,40};
        Box<2,double> bx({0.25 + it.getSpacing(0)/4.0,0.25 + it.getSpacing(0)/4.0},{sz2[0]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0, sz2[1]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<2,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[1] = y;

            ++it2;
        }
        }

        // Add multi res patch 2

        {
        const size_t sz2[2] = {40,40};
        Box<2,double> bx({0.25 + 21.0*spacing/8.0,0.25 + 21.0*spacing/8.0},{sz2[0]*spacing/4.0 + 0.25 + 21.0*spacing/8.0, sz2[1]*spacing/4.0 + 0.25 + 21*spacing/8.0});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<2,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[1] = y;

            ++it2;
        }
        }

        ///////////////////////

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        Derivative_x Dx(domain, 2, rCut / 3.0 ,1.9/*,support_options::RADIUS*/);
        Derivative_y Dy(domain, 2, rCut / 3.0 ,1.9/*,support_options::RADIUS*/);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut / 3.0 ,1.9/*,support_options::RADIUS*/);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);


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

        domain.ghost_get<1,3>();

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);
        auto Poisson = Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(D_y, up_p, 0);
        Solver.impose(D_x, r_p, 0);
        Solver.impose(v, dw_p, 0);
        Solver.impose(v, l_p, 0);

        petsc_solver<double> solver;

        solver.setRestart(500);
	solver.log_monitor();

//        auto A = Solver.getA();
//        A.write("Matrix_anal_sol");

//        solver.setSolver(KSPGMRES);
//        solver.setPreconditioner(PCKSP);
//        solver.setRestart(300);

/*        solver.setPreconditioner(PCHYPRE_BOOMERAMG);
        solver.setPreconditionerAMG_nl(6);
        solver.setPreconditionerAMG_maxit(10);
        solver.setPreconditionerAMG_relax("SOR/Jacobi");
        solver.setPreconditionerAMG_cycleType("V",0,4);
        solver.setPreconditionerAMG_coarsenNodalType(0);
        solver.setPreconditionerAMG_coarsen("HMIS");*/

//        solver.log_monitor();

        Solver.solve_with_solver(solver,sol);

        domain.ghost_get<2>();

        DCPSE_sol=Lap(sol);
        domain.ghost_get<5>();

        double worst1 = 0.0;

        v=abs(DCPSE_sol-RHS);

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<3>(p) - domain.getProp<2>(p));

        }
        std::cout << "Maximum Analytic Error: " << worst1 << std::endl;

        domain.ghost_get<4>();
        domain.write("Robin_anasol");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_poisson_Dirichlet_anal) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {81,81};
        Box<2, double> box({0, 0}, {1, 1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
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

        Derivative_x Dx(domain, 2, rCut,2);
        Derivative_y Dy(domain, 2, rCut,2);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut, 2);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);

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
        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);
        auto Poisson = Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(v, up_p, 0);
        Solver.impose(v, r_p, 0);
        Solver.impose(v, dw_p, 0);
        Solver.impose(v, l_p, 0);
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
        std::cout << "Maximum Analytic Error: " << worst1 << std::endl;

        domain.write("Dirichlet_anasol");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_poisson_Periodic) {
        //https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/periodic/python/documentation.html
        //  int rank;
        //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3.1);
        double rCut = 3.1 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> domain(0, box, bc, ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1)*0.99999;
            domain.getLastPos()[1] = y;

            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();


        Laplacian Lap(domain, 2, rCut, 1.9, support_options::RADIUS);

        DCPSE_scheme<equations2d1p,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto u = getV<5>(domain);

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
                domain.getProp<1>(p) =  0;
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  0;
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = xp.get(0)*sin(5*M_PI*xp.get(1))+exp(-((xp.get(0)-0.5)*(xp.get(0)-0.5)+(xp.get(1)-0.5)*(xp.get(1)-0.5))/0.02);
            }

            ++it2;
        }

        domain.ghost_get<1>();
        auto Poisson = -Lap(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(v, up_p, 0);
        Solver.impose(v, dw_p, 0);
        Solver.solve(v);

        domain.ghost_get<0>();
        anasol=-Lap(v);
        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<1>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<1>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<3>(p));

        }
        std::cout << "Maximum Auto Error: " << worst1 << std::endl;

        domain.write("Poisson_Periodic");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_poisson_Robin) {
        //http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation#Full_Neumann_boundary_conditions
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {150,150};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> domain(0, box, bc, ghost);


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

        //Derivative_x Dx(domain, 2, rCut,2);
        Derivative_y Dy(domain, 2, rCut,2);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut, 3);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);
        //DCPSE_scheme<equations,decltype(domain)> Solver( domain);
        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        //auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto u = getV<5>(domain);

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
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5.0*xp.get(0));
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -10.0*exp(-((xp.get(0)-0.5)*(xp.get(0)-0.5)+(xp.get(1)-0.5)*(xp.get(1)-0.5))/0.02);
            }

            ++it2;
        }


        auto Poisson = Lap(v);
        //auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(D_y, up_p, prop_id<1>());
        Solver.impose(-D_y, dw_p, prop_id<1>());
        Solver.impose(v, l_p, 0);
        Solver.impose(v, r_p, 0);
        Solver.solve(sol);
        anasol=Lap(sol);
        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<1>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<1>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<3>(p));

        }
        std::cout << "Maximum Auto Error: " << worst1 << std::endl;

        domain.write("Mixed");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_poisson_Neumann) {
    //https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/neumann-poisson/python/documentation.html
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {81,81};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 0.5 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double>> domain(0, box, bc, ghost);


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

        Derivative_x Dx(domain, 2, rCut,2);
        Derivative_y Dy(domain, 2, rCut,2);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 3, rCut, 2);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);
        //DCPSE_scheme<equations,decltype(domain)> Solver( domain);
        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain,options_solver::LAGRANGE_MULTIPLIER);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        //auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);

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
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  sin(5*xp.get(0));
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -10*exp(-((xp.get(0)-0.5)*(xp.get(0)-0.5)+(xp.get(1)-0.5)*(xp.get(1)-0.5))/0.02);
            }

            ++it2;
        }


        auto Poisson = -Lap(v);
        auto D_x = Dx(v);
        auto D_y = Dy(v);
        Solver.impose(Poisson, bulk, prop_id<1>());
        Solver.impose(D_y, up_p, prop_id<1>());
        Solver.impose(-D_y, dw_p, prop_id<1>());
        Solver.impose(-D_x, l_p, prop_id<1>());
        Solver.impose(D_x, r_p, prop_id<1>());
        Solver.solve(sol);
        anasol=-Lap(sol);
        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p) - domain.getProp<1>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p) - domain.getProp<1>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<3>(p));

        }
        std::cout << "Maximum Auto Error: " << worst1 << std::endl;

        domain.write("Neumann");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    BOOST_AUTO_TEST_CASE(dcpse_op_solver) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31, 31};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double>> domain(0, box, bc, ghost);


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

        //Derivative_x Dx(domain, 2, rCut);
        //Derivative_y Dy(domain, 2, rCut);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut, 2);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);
        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);

        // Here fill me

        Box<2, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) - spacing / 2.0, box.getHigh(1) - 8*spacing / 2.0},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - 6*spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) - spacing / 2.0, box.getLow(1) + 3*spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + 4*spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(up_d);
        boxes.add(down);
        boxes.add(down_u);
        boxes.add(left);
        boxes.add(left_r);
        boxes.add(right);
        boxes.add(right_l);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);
            domain.getProp<2>(p)=1+xp[0]*xp[0]+2*xp[1]*xp[1];
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  3 + xp.get(0)*xp.get(0);
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  1 + xp.get(0)*xp.get(0);
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  1 + 2*xp.get(1)*xp.get(1);
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  2 + 2*xp.get(1)*xp.get(1);
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }


            ++it2;
        }

//        bulk;
//        boundaries1;
//        boundaries2;

        auto eq1 = Lap(v);
        //auto flux = Dx(v) + v;



        Solver.impose(eq1, bulk, 6);
        Solver.impose(v, up_p, prop_id<1>());
        Solver.impose(v, dw_p, prop_id<1>());
        Solver.impose(v, l_p, prop_id<1>());
        Solver.impose(v, r_p, prop_id<1>());
        Solver.solve(v);
        anasol=Lap(v);

        double worst1 = 0.0;

        it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(domain.getProp<0>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<0>(p) - domain.getProp<2>(p));
            }

            domain.getProp<1>(p) = fabs(domain.getProp<0>(p) - domain.getProp<2>(p));

            ++it2;
        }

        std::cout << "Maximum Error: " << worst1 << std::endl;

        domain.write("particles");
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    BOOST_AUTO_TEST_CASE(dcpse_op_vec) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = box.getHigh(0)  / (sz[0] - 1);
        spacing[1] = box.getHigh(1) / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3.1);
        double rCut = 3.1 * spacing[0];
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
            domain.template getLastProp<0>()    = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            domain.template getLastProp<1>()[0] = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            domain.template getLastProp<1>()[1] = cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]);
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

            domain.template getLastProp<5>()    = cos(domain.getLastPos()[0])*(sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]))+cos(domain.getLastPos()[1]) * (cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]));



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
        Advection Adv(domain, 2, rCut, 1.9,support_options::RADIUS);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);
        auto dv = getV<3>(domain);
        auto dP = getV<6>(domain);


//        typedef boost::mpl::int_<std::is_fundamental<point_expression_op<Point<2U, double>, point_expression<double>, Point<2U, double>, 3>>::value>::blabla blabla;

//        std::is_fundamental<decltype(o1.value(key))>

        //vv=Lap(P);
        //dv=Lap(v);//+Dy(P);
        dv = Adv(v, v);//+Dy(P);
        auto it2 = domain.getDomainIterator();

        double worst1 = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();

            //std::cout << "VALS: " << domain.getProp<3>(p)[0] << " " << domain.getProp<4>(p)[0] << std::endl;
            //std::cout << "VALS: " << domain.getProp<3>(p)[1] << " " << domain.getProp<4>(p)[1] << std::endl;

            //domain.getProp<0>(p)=std::sqrt((domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0])*(domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0])+(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1])*(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1]));

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

            //std::cout << "VALS: " << domain.getProp<3>(p)[0] << " " << domain.getProp<4>(p)[0] << std::endl;
            //std::cout << "VALS: " << domain.getProp<3>(p)[1] << " " << domain.getProp<4>(p)[1] << std::endl;

            //domain.getProp<0>(p)=std::sqrt((domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0])*(domain.getProp<3>(p)[0] - domain.getProp<4>(p)[0])+(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1])*(domain.getProp<3>(p)[1] - domain.getProp<4>(p)[1]));

            if (fabs(domain.getProp<6>(p) - domain.getProp<5>(p)) > worst2) {
                worst2 = fabs(domain.getProp<6>(p) - domain.getProp<5>(p));

            }

            ++it3;
        }

        std::cout << "Maximum Error: " << worst2 << std::endl;

        //Adv.checkMomenta(domain);
        //Adv.DrawKernel<2>(domain,0);

        domain.deleteGhost();
        domain.write("v2");



        //std::cout << demangle(typeid(decltype(v)).name()) << "\n";

        //Debug<decltype(expr)> a;

        //typedef decltype(expr)::blabla blabla;

        //auto err = Dx + Dx;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    BOOST_AUTO_TEST_CASE(dcpse_op_div) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 80;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 2.0 * spacing[0];
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
        Divergence Div(domain, 2, rCut, 3);
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_op_tests) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {2 * edgeSemiSize, 2 * edgeSemiSize};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 2.0 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, double, double, VectorS<2, double>, VectorS<2, double>>> domain(0, box,
                                                                                                                 bc,
                                                                                                                 ghost);

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
            domain.template getLastProp<0>() = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
//            domain.template getLastProp<0>() = x * x;
//            domain.template getLastProp<0>() = x;
            // Here fill the validation value for Df/Dx
            domain.template getLastProp<2>() = cos(domain.getLastPos()[0]) + cos(domain.getLastPos()[1]);
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

        Derivative_x Dx(domain, 2, rCut);
        Derivative_y Dy(domain, 2, rCut);
        Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);
        auto vv = getV<3>(domain);

        vv = Lap(P);
        v = Dx(P) + Dy(P);
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

        domain.deleteGhost();
        domain.write("v");

        //std::cout << demangle(typeid(decltype(v)).name()) << "\n";

        //Debug<decltype(expr)> a;

        //typedef decltype(expr)::blabla blabla;

        //auto err = Dx + Dx;
    }

    BOOST_AUTO_TEST_CASE(dcpse_test_diffusion) {
        const size_t sz[2] = {100, 100};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double[2]>> domain(0, box, bc, ghost);

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

        //Derivative_x Dx(domain, 2, rCut);
        //Derivative_y Dy(domain, 2, rCut);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut, 3);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);
        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        auto v = getV<0>(domain);
        auto dc = getV<1>(domain);
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
        //VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        //vtk_box.add(boxes);
        //vtk_box.write("vtk_box.vtk");

        domain.write("Diffusion");

        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = domain.getPos(p);

            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  0;
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  0;
            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = 0;
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  0;
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                if(xp[0]==0.5 && xp[1]==0.5)
                    domain.getProp<1>(p) =  0;

            }

            ++it2;
        }


/*        double t=0;
        while (t<2)
        {
            dc=Lap(v);
            t=t+0.1;
            v=dc;
            domain.deleteGhost();
            domain.write("Diffusion");
        }*/

        //auto flux = Dx(v) + v;

    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_op_test_lap) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        size_t edgeSemiSize = 160;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing[2];
        spacing[0] = 2 * M_PI / (sz[0] - 1);
        spacing[1] = 2 * M_PI / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 2.0 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        vector_dist<2, double, aggregate<double, double, double, double, VectorS<2, double>>> domain(0, box,
                                                                                                                 bc,
                                                                                                                 ghost);

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
            domain.template getLastProp<0>() = sin(domain.getLastPos()[0]) + sin(domain.getLastPos()[1]);
            // Here fill the validation value for Lap
            domain.template getLastProp<1>() = - sin(domain.getLastPos()[0]) - sin(domain.getLastPos()[1]);
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

        //Derivative_x Dx(domain, 2, rCut);
        //Derivative_y Dy(domain, 2, rCut);
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut,1.9);
        auto v = getV<1>(domain);
        auto P = getV<0>(domain);
        auto vv = getV<2>(domain);
        auto errv = getV<3>(domain);



        vv = Lap(P);
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
        errv=v-vv;
        domain.deleteGhost();
        domain.write("v");

        //std::cout << demangle(typeid(decltype(v)).name()) << "\n";

        //Debug<decltype(expr)> a;

        //typedef decltype(expr)::blabla blabla;

        //auto err = Dx + Dx;
    }

    BOOST_AUTO_TEST_CASE(dcpse_slice) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2.0 * spacing;

        vector_dist<2, double, aggregate<double,VectorS<2, double>,double,double[3][3]>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;

            Particles.getLastProp<1>()[0] = sin(x+y);
            Particles.getLastProp<1>()[1] = cos(x+y);

            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto S = getV<2>(Particles);
        auto Sig = getV<3>(Particles);


        Derivative_x Dx(Particles, 2, rCut,2);

        P = Dx(V[0]);
        S = V[0]*V[0] + V[1]*V[1];

        Sig[0][1] = V[0]*V[0] + V[1]*V[1];
        Sig[1][0] = P;
        Sig[2][2] = 5.0;

        auto it2 = Particles.getDomainIterator();

		double err = 0.0;

		while (it2.isNext())
		{
			auto p = it2.get();

			if (fabs(Particles.getProp<0>(p) - Particles.getProp<1>(p)[1]) >= err )
			{
				err = fabs(Particles.getProp<0>(p) - Particles.getProp<1>(p)[1]);
			}

			if (fabs(Particles.getProp<2>(p) - 1.0) >= err )
			{
				err = fabs(Particles.getProp<2>(p) - 1.0);
			}

			if (fabs(Particles.getProp<3>(p)[0][1] - 1.0) >= err )
			{
				err = fabs(Particles.getProp<3>(p)[0][1] - 1.0);
			}

			if (fabs(Particles.getProp<3>(p)[1][0] - Particles.getProp<1>(p)[1]) >= err )
			{
				err = fabs(Particles.getProp<3>(p)[0][1] - Particles.getProp<1>(p)[1]);
			}

			if (fabs(Particles.getProp<3>(p)[2][2] - 5.0) >= err )
			{
				err = fabs(Particles.getProp<3>(p)[2][2] - 5.0);
			}

			++it2;
		}

        Particles.write("test_out");
    }



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BOOST_AUTO_TEST_CASE(dcpse_slice_solver) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {1, 1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<2, double> ghost(spacing * 3);
        double rCut = 2 * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>,VectorS<2, double>>> domain(0, box, bc, ghost);


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
        //Gradient Grad(domain, 2, rCut);
        Laplacian Lap(domain, 2, rCut,1.9,support_options::RADIUS);
        //Advection Adv(domain, 3, rCut, 3);
        //Solver Sol_Lap(Lap),Sol_Dx(Dx);


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
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));


            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p)[0] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[0] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

                domain.getProp<1>(p)[1] = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p)[1] = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            }
            ++it2;
        }

        eq_id vx,vy;

        vx.setId(0);
        vy.setId(1);

        DCPSE_scheme<equations2d1,decltype(domain)> Solver( domain);
        auto Poisson0 = Lap(v[0]);
        auto Poisson1 = Lap(v[1]);
        //auto D_x = Dx(v[1]);
        //auto D_y = Dy(v[1]);
        Solver.impose(Poisson0, bulk, RHS[0],vx);
        Solver.impose(Poisson1, bulk, RHS[1],vy);
        Solver.impose(v[0], up_p, RHS[0],vx);
        Solver.impose(v[1], up_p, RHS[1],vy);
        Solver.impose(v[0], r_p,  RHS[0],vx);
        Solver.impose(v[1], r_p,  RHS[1],vy);
        Solver.impose(v[0], dw_p, RHS[0],vx);
        Solver.impose(v[1], dw_p, RHS[1],vy);
        Solver.impose(v[0], l_p,  RHS[0],vx);
        Solver.impose(v[1], l_p,  RHS[1],vy);
        Solver.solve(sol[0],sol[1]);
        DCPSE_sol=Lap(sol);
        double worst1 = 0.0;
        double worst2 = 0.0;


        v=sol-RHS;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p)[0] - domain.getProp<2>(p)[0]) >= worst1) {
                worst1 = fabs(domain.getProp<3>(p)[0] - domain.getProp<2>(p)[0]);
            }
            domain.getProp<4>(p)[0] = fabs(domain.getProp<3>(p)[0] - domain.getProp<2>(p)[0]);

        }
        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<3>(p)[1] - domain.getProp<2>(p)[1]) >= worst2) {
                worst2 = fabs(domain.getProp<3>(p)[1] - domain.getProp<2>(p)[1]);
            }
            domain.getProp<4>(p)[1] = fabs(domain.getProp<3>(p)[1] - domain.getProp<2>(p)[1]);

        }
        std::cout << "Maximum Analytic Error in slice x: " << worst1 << std::endl;
        std::cout << "Maximum Analytic Error in slice y: " << worst2 << std::endl;

        domain.write("Slice_anasol");
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE_END()


