//
// Created by Abhinav Singh on 27.10.20.
//

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

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
#include "Decomposition/Distribution/SpaceDistribution.hpp"
#include "../SphericalHarmonics.hpp"


BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests_paper)

    BOOST_AUTO_TEST_CASE(Lid_cavity_paper) {
        /*Data from Ghia et al.
         * y,u
        1.0000,1.00000
        0.9766,0.84123
        0.9688,0.78871
        0.9609,0.73722
        0.9531,0.68717
        0.8516,0.23151
        0.7344,0.00332
        0.6172,-0.13641
        0.5000,-0.20581
        0.4531,-0.21090
        0.2813,-0.15662
        0.1719,-0.10150
        0.1016,-0.06434
        0.0703,-0.04775
        0.0625,-0.04192
        0.0547,-0.03717
        0.0000,0.00000*/
        constexpr int x = 0;
        constexpr int y = 1;
        size_t gd_sz = 256;
        const size_t sz[2] = {gd_sz,gd_sz};
        Box<2, double> box({0, 0}, {1,1});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing;
        spacing = 1.0 / (sz[0] - 1);
        double rCut = 4.1 * spacing;
        int ord = 2;
        double sampling_factor = 4.0;
        double alpha=0.0125;
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");
        auto &v_cl = create_vcluster();
        vector_dist<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,double>> Particles(0, box,
                                                                                                                                                                           bc,
                                                                                                                                                                           ghost);
        //Init_DCPSE(Particles)
        BOOST_TEST_MESSAGE("Init Particles...");

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;
        openfpm::vector<aggregate<int>> corner_ul;
        openfpm::vector<aggregate<int>> corner_ur;
        openfpm::vector<aggregate<int>> corner_dl;
        openfpm::vector<aggregate<int>> corner_dr;


        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing;
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing;
            Particles.getLastPos()[1] = yp0;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync Particles across processors...");
        Particles.map();
        Particles.ghost_get<0>();

        double x0=box.getLow(0),x1=box.getHigh(0),y0=box.getLow(1),y1=box.getHigh(1);

        Box<2, double> up({x0 - spacing / 2.0, y1 - spacing / 2.0},
                          {x1 + spacing / 2.0, y1 + spacing / 2.0});

        Box<2, double> down({x0 - spacing / 2.0, y0 - spacing / 2.0},
                            {x1 + spacing / 2.0, y0 + spacing / 2.0});

        Box<2, double> left({x0 - spacing / 2.0, y0 - spacing / 2.0},
                            {x0 + spacing / 2.0, y1 + spacing / 2.0});

        Box<2, double> right({x1 - spacing / 2.0, y0 - spacing / 2.0},
                             {x1 + spacing / 2.0, y1 + spacing / 2.0});
        /*Box<2, double> mid({box.getHigh(0) / 2.0 - spacing, box.getHigh(1) / 2.0 - spacing / 2.0},
                           {box.getHigh(0) / 2.0, box.getHigh(1) / 2.0 + spacing / 2.0});
*/


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);


        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("vtk_box.vtk");

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            if (up.isInside(xp) == true) {
                if (left.isInside(xp) == true) {
                    corner_ul.add();
                    corner_ul.last().get<0>() = p.getKey();
                } else if (right.isInside(xp) == true) {
                    corner_ur.add();
                    corner_ur.last().get<0>() = p.getKey();
                } else {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                }
                up_p.add();
                up_p.last().get<0>() = p.getKey();

            } else if (down.isInside(xp) == true) {
                if (left.isInside(xp) == true) {
                    corner_dl.add();
                    corner_dl.last().get<0>() = p.getKey();
                } else if (right.isInside(xp) == true) {
                    corner_dr.add();
                    corner_dr.last().get<0>() = p.getKey();
                } else {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                }
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                if (up.isInside(xp) == false && down.isInside(xp) == false) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                }
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                if (up.isInside(xp) == false && down.isInside(xp) == false) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                }
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }

        vector_dist_subset<2, double, aggregate<double, VectorS<2, double>, VectorS<2, double>, VectorS<2, double>,double,VectorS<2, double>,double,double>> Particles_bulk(Particles,bulk);

        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto RHS = getV<2>(Particles);
        auto dV = getV<3>(Particles);
        auto div = getV<4>(Particles);
        auto V_star = getV<5>(Particles);
        auto H = getV<6>(Particles);
        auto dP = getV<7>(Particles);



        auto P_bulk = getV<0>(Particles_bulk);
        auto V_bulk = getV<1>(Particles_bulk);
        auto RHS_bulk =getV<2>(Particles_bulk);
        auto V_star_bulk = getV<5>(Particles_bulk);
        auto dP_bulk = getV<7>(Particles_bulk);


        P_bulk = 0;
        Advection Adv(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Dx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        //Laplacian Lap(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_xx Dxx(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);

        Derivative_y Dy(Particles, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_x Bulk_Dx(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);
        Derivative_y Bulk_Dy(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);
        Laplacian Bulk_Lap(Particles_bulk, 2, rCut,sampling_factor, support_options::RADIUS);

        int n = 0, nmax = 300, ctr = 0, errctr=1, Vreset = 0;
        double V_err=1;
        if (Vreset == 1) {
            P_bulk = 0;
            P = 0;
            Vreset = 0;
        }
        P=0;

        eq_id vx,vy;
        vx.setId(0);
        vy.setId(1);
/*          W.setVarId(0);
        Sf.setVarId(1);
*/

        double Re=100;
        double sum, sum1, sum_k,V_err_eps=1e-5,V_err_old;
        auto StokesX=V[x]*Dx(V_star[x])-(1.0/Re)*(Dxx(V_star[x])+Dyy(V_star[x]));
        auto StokesY=V[y]*Dy(V_star[y])-(1.0/Re)*(Dxx(V_star[y])+Dyy(V_star[y]));
        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc.setPreconditioner(PCJACOBI);
        petsc_solver<double> solverPetsc2;
        solverPetsc2.setSolver(KSPGMRES);
        //solverPetsc.setRestart(250);
        solverPetsc2.setPreconditioner(PCJACOBI);
        RHS[x] = V[x];
        RHS[y] = V[y];
        dV=0;
        timer tt;
        while (V_err >= V_err_eps && n <= nmax) {
            if (n%5==0){
            Particles.ghost_get<0,1,2,3,4,5,6,7>(SKIP_LABELLING);
            Particles.deleteGhost();
            Particles.write_frame("LID",n,BINARY);
            Particles.ghost_get<0>();
            }
            tt.start();
            Particles.ghost_get<0>(SKIP_LABELLING);
            RHS_bulk[x] = -Bulk_Dx(P);
            RHS_bulk[y] = -Bulk_Dy(P);
            DCPSE_scheme<equations2d2E, decltype(Particles)> Solver(Particles);
            Solver.impose(StokesX, bulk, RHS[0], vx);
            Solver.impose(StokesY, bulk, RHS[1], vy);
            Solver.impose(V_star[x], up_p1, 1.0, vx);
            Solver.impose(V_star[y], up_p1, 0.0, vy);
            Solver.impose(V_star[x], l_p, 0.0, vx);
            Solver.impose(V_star[y], l_p, 0.0, vy);
            Solver.impose(V_star[x], r_p, 0.0, vx);
            Solver.impose(V_star[y], r_p, 0.0, vy);
            Solver.impose(V_star[x], dw_p, 0.0, vx);
            Solver.impose(V_star[y], dw_p, 0.0, vy);
            Solver.impose(V_star[x], corner_ul, 0.0, vx);
            Solver.impose(V_star[y], corner_ul, 0.0, vy);
            Solver.impose(V_star[x], corner_ur, 0.0, vx);
            Solver.impose(V_star[y], corner_ur, 0.0, vy);
            //Solver.solve_with_solver(solverPetsc, V_star[x], V_star[y]);
            Solver.solve(V_star[x], V_star[y]);
            //Particles.write_frame("LIDV",n);
            Particles.ghost_get<5>(SKIP_LABELLING);
            div = (Dx(V_star[x]) + Dy(V_star[y]));
            DCPSE_scheme<equations2d1E,decltype(Particles)> SolverH(Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Dxx(H)+Dyy(H);
            SolverH.impose(Helmholtz,bulk,prop_id<4>());
/*            SolverH.impose(H, up_p,0);
            SolverH.impose(H, l_p,0);
            SolverH.impose(H, r_p,0);
            SolverH.impose(H, dw_p,0);*/
            SolverH.impose(Dy(H), up_p,0);
             SolverH.impose(Dx(H), l_p,0);
             SolverH.impose(Dx(H), r_p,0);
             SolverH.impose(Dy(H), dw_p,0);
             /*SolverH.impose(-Dx(H) + Dy(H), corner_ul, 0);
             SolverH.impose(Dx(H) + Dy(H), corner_ur, 0);
             SolverH.impose(-Dx(H) - Dy(H), corner_dl, 0);
             SolverH.impose(Dx(H) - Dy(H), corner_dr, 0);*/
            //SolverH.solve_with_solver(solverPetsc2,H);
            SolverH.solve(H);
            Particles.ghost_get<6>(SKIP_LABELLING);
            //dP_bulk=Bulk_Lap(H);
            P_bulk = P - alpha*div;
            //dV[x]=Dx(H);
            //dV[y]=Dy(H);
            V_star_bulk[0] = V_star[0] - Bulk_Dx(H);
            V_star_bulk[1] = V_star[1] - Bulk_Dy(H);
            sum = 0;
            sum1 = 0;

            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) *
                       (Particles.getProp<5>(p)[0] - Particles.getProp<1>(p)[0]) +
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]) *
                       (Particles.getProp<5>(p)[1] - Particles.getProp<1>(p)[1]);
                sum1 += Particles.getProp<5>(p)[0] * Particles.getProp<5>(p)[0] +
                        Particles.getProp<5>(p)[1] * Particles.getProp<5>(p)[1];
            }

            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V = V_star;
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            V_err_old = V_err;
            V_err = sum / sum1;
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
                errctr++;
                //alpha_P -= 0.1;
            } else {
                errctr = 0;
            }
            if (n > 3) {
                if (errctr > 1) {
                    //std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN ERROR" << std::endl;
                    std::cout << "Alpha Halfed DUE TO INCREASE/VERY SLOW DECREASE IN ERROR" << std::endl;
                    alpha=alpha/2;
                    Vreset = 1;
                    errctr=0;
                    //break;
                } else {
                    Vreset = 0;
                }
            }
            n++;
            tt.stop();
            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << " and took " <<tt.getwct() <<"("<<tt.getcputime()<<") seconds(CPU)." << std::endl;
            }
        }
        double worst1 = 0.0;
        double worst2 = 0.0;

      /*  for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(Particles.getProp<6>(p)[0] - Particles.getProp<1>(p)[0]) >= worst1) {
                worst1 = fabs(Particles.getProp<6>(p)[0] - Particles.getProp<1>(p)[0]);
            }
        }
        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(Particles.getProp<6>(p)[1] - Particles.getProp<1>(p)[1]) >= worst2) {
                worst2 = fabs(Particles.getProp<6>(p)[1] - Particles.getProp<1>(p)[1]);
            }
        }*/
        //Particles.deleteGhost();
        //Particles.write("PC_subset_lid");
/*        std::cout << "Maximum Analytic Error in Vx: " << worst1 << std::endl;
        std::cout << "Maximum Analytic Error in Vy: " << worst2 << std::endl;*/
        Particles.deleteGhost();
        Particles.write("LID");

    }
    BOOST_AUTO_TEST_CASE(Active2d_paper) {

    }
    BOOST_AUTO_TEST_CASE(Sphere3d_paper) {
        size_t grd_sz=17;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({-1.0, -1.0,-1.0}, {1.0,1.0,1.0});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = 1.0 / (sz[0] - 1);
        double rCut = 3.9 * spacing;
        double R=1.0;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_B           RHS            V_t   Helmholtz              RHS2            Polar_cord
        vector_dist<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles(0, box, bc, ghost);

        std::mt19937 generator(1);
        std::uniform_real_distribution<double> uniform01(0.0, 1.0);
        auto &v_cl = create_vcluster();


        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            auto key = it.get();
            double x = -1.0+key.get(0) * it.getSpacing(0);
            double y = -1.0+key.get(1) * it.getSpacing(1);
            double z = -1.0+key.get(2) * it.getSpacing(2);
            double r=sqrt(x*x+y*y+z*z);
            if (r<R) {
                Particles.add();
                Particles.getLastPos()[0] = x;
                Particles.getLastPos()[1] = y;
                Particles.getLastPos()[2] = z;
            }
            ++it;
        }
        for(int i=0;i<=2*int(grd_sz);i++)
        {
            for(int j=0;j<=1.5*int(grd_sz);j++)
            {
                double theta,phi;
                phi = (2 * M_PI - spacing) * i* spacing;
                theta = (M_PI - spacing) * j * spacing;
                double x = R * sin(theta) * cos(phi);
                double y = R * sin(theta) * sin(phi);
                double z = R * cos(theta);
                Particles.add();
                Particles.getLastPos()[0] = x;
                Particles.getLastProp<8>()[0] = R;
                Particles.getLastPos()[1] = y;
                Particles.getLastProp<8>()[1] = theta;
                Particles.getLastPos()[2] = z;
                Particles.getLastProp<8>()[2] = phi;
            }
        }
        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> bulk_without_center;
        openfpm::vector<aggregate<int>> Surface;
        openfpm::vector<aggregate<int>> temp_part;
        openfpm::vector<aggregate<int>> center;

        auto it2 = Particles.getDomainIterator();
        Particles.ghost_get<0,1,2,3,4,5>();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Point<3, double> xP = Particles.getProp<8>(p);
            Particles.getProp<0>(p) =0;
            if (sqrt(xp[0]*xp[0]+xp[1]*xp[1]+xp[2]*xp[2]) < R-spacing/2) {
                if (xp[0]==0 && xp[1]==0 && xp[1]==0)
                {   center.add();
                    center.last().get<0>() = p.getKey();
                }
                else{
                    bulk_without_center.add();
                    bulk_without_center.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) =  0;
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;

            }
            else {
                Surface.add();
                Surface.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) =  0;
                //if(xp[2]>0 && xp[2]<spacing) {
                //Particles.getProp<1>(p)[0] = boost::math::spherical_harmonic_r(0, 0,xP[1] ,xP[2]);
                Particles.getProp<1>(p)[0] =  -xp[1];//(pow(xp[0],2)+pow(xp[1],2));
                Particles.getProp<1>(p)[1] =  xp[0];//(pow(xp[0],2)+pow(xp[1],2));
                Particles.getProp<1>(p)[2] =  0;
                //}
            }
            ++it2;
        }

        vector_dist_subset<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles_bulk(Particles,bulk);
        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto V_B = getV<2>(Particles);
        V.setVarId(0);
        auto DIV = getV<3>(Particles);
        auto V_t = getV<4>(Particles);
        auto temp=getV<6>(Particles);
        auto RHS = getV<7>(Particles);
        auto P_bulk = getV<0>(Particles_bulk);
        auto RHS_bulk = getV<7>(Particles_bulk);

        std::cout << "Sphere Init" << std::endl;
        V_t=V;
        V_B=V;
        P=0;
        P_bulk=0;
        Particles.write("Sphere.vtk");
        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        double sampling=1.9;
        double sampling2=1.9;
        double rCut2=3.1*spacing;

        Derivative_x Dx(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dx(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dy(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dz(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        //Laplacian Lap(Particles, 2, rCut,sampling2);
        Derivative_xx Dxx(Particles, 2, rCut2,sampling2);
        Derivative_yy Dyy(Particles, 2, rCut2,sampling2);
        Derivative_zz Dzz(Particles, 2, rCut2,sampling2);

/*        std::cout << "Dx" << std::endl;
        Dx.checkMomenta(Particles);
        std::cout << "Dy" << std::endl;
        Dy.checkMomenta(Particles);
        std::cout << "Dz" << std::endl;
        Dz.checkMomenta(Particles);
        std::cout << "Lap" << std::endl;
        Lap.checkMomenta(Particles);
        Particles.write("Sphere.vtk");
        Dx.DrawKernel<0>(Particles,25828);*/
        //B_Dx.checkMomenta(Particles_subset);
        Particles.write("Sphere.vtk");
        std::cout << "DCPSE KERNELS DONE" << std::endl;
        petsc_solver<double> solverPetsc;
        //solverPetsc.setRestart(250);
        double nu=0.5;
        Particles.ghost_get<0,1,2,3,4,5,6,7>();
        timer tt;
        double sum=0,sum1=0;
        V_t=V;
        double V_err_eps = 1e-3;
        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 30;
        int ctr = 0, errctr, Vreset = 0;
        V_err = 1;
        Particles.write_frame("StokesSphere",0);
        n = 1;
        tt.start();
        while (V_err >= V_err_eps && n <= nmax) {
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
            std::cout << "Stokes Solved" << std::endl;
            Particles.ghost_get<1>();
            DIV = -(Dx(V[0])+Dy(V[1])+Dz(V[2]));
           /* DCPSE_scheme<equations3d1, decltype(Particles)> SolverH(Particles);//,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Lap(H);
            SolverH.impose(Helmholtz, bulk, prop_id<3>());
            SolverH.impose(H, Surface, 0);*/
/*            for(int j=0;j<Surface.size();j++)
            {   temp_part.clear();
                auto p=Surface.get<0>(j);
                temp_part.add();
                temp_part.last().get<0>() = Surface.template get<0>(j);
                Point<3, double> xp = Particles.getPos(p);
                SolverH.impose(xp[0]*Dx(H)+xp[1]*Dy(H)+xp[2]*Dz(H),temp_part,0);
            }*/
            //SolverH.solve(H);
            //SolverH.solve_with_solver(solverPetsc2, H);
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
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
                errctr++;
            } else {
                errctr = 0;
            }
            if (n > 3) {
                if (errctr > 3) {
                    std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN ERROR" << std::endl;
                    Vreset = 1;
                    break;
                } else {
                    Vreset = 0;
                }
            }
            n++;
            //Particles.write_frame("V_debug", n);
            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;
            }
        }
        tt.stop();
        if (v_cl.rank() == 0) {
            std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds."
                      << std::endl;
        }
        Particles.write_frame("StokesSphere",1);

    }

    BOOST_AUTO_TEST_CASE(Sph_harm) {
    BOOST_REQUIRE(boost::math::Y(2,1,0.5,0)+0.459674<0.001);
    BOOST_REQUIRE(boost::math::legendre_p(0,-1,1)==0.00);
    //TO BE TESTED...
    double a=boost::math::legendre_p(0,-1,1);
    std::cout<<a<<std::endl;
        double nu=0.5;
        size_t grd_sz=17;
        const size_t sz[3] = {grd_sz,grd_sz,grd_sz};
        Box<3, double> box({-1.0, -1.0,-1.0}, {1.0,1.0,1.0});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double spacing = 2.0 / (sz[0] - 1);
        double rCut = 3.9 * spacing;
        double R=1.0;
        Ghost<3, double> ghost(rCut);
        //                                  P        V                 v_B           RHS            V_t         P_anal              RHS2            Polar cord
        vector_dist<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles(0, box, bc, ghost);

        std::mt19937 generator(1);
        std::uniform_real_distribution<double> uniform01(0.0, 1.0);
        auto &v_cl = create_vcluster();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> Surface;

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            auto key = it.get();
            double x = -1.0+key.get(0) * it.getSpacing(0);
            double y = -1.0+key.get(1) * it.getSpacing(1);
            double z = -1.0+key.get(2) * it.getSpacing(2);
            double r=sqrt(x*x+y*y+z*z);
            if (r<R) {
                Particles.add();
                Particles.getLastPos()[0] = x;
                Particles.getLastPos()[1] = y;
                Particles.getLastPos()[2] = z;
                Particles.getLastProp<8>()[0] = r;
                if (r==0){
                Particles.getLastProp<8>()[1] = 0.0;
                }
                else{
                    Particles.getLastProp<8>()[1] = acos(z/sqrt(r));
                }
                Particles.getLastProp<8>()[2] = std::atan2(y,x);
            }
            ++it;
        }

        int n_sp=int(grd_sz)*int(grd_sz)*6.0;

        double Golden_angle=M_PI * (3.0 - sqrt(5.0));

        for(int i=1;i<=n_sp;i++)
        {
                double y = 1 - (i /(n_sp - 1.0)) * 2.0;
                double radius = sqrt(1 - y * y);
                double Golden_theta = Golden_angle * i;
                double x = cos(Golden_theta) * radius;
                double z = sin(Golden_theta) * radius;
                Particles.add();
                Particles.getLastPos()[0] = x;
                Particles.getLastPos()[1] = y;
                Particles.getLastPos()[2] = z;
                Particles.getLastProp<8>()[0] = 1.0 ;
                Particles.getLastProp<8>()[1] = acos(z);
                Particles.getLastProp<8>()[2] = std::atan2(y,x);
            }

        Particles.add();
        Particles.getLastPos()[0] = 0;
        Particles.getLastPos()[1] = 0;
        Particles.getLastPos()[2] = 1;
        Particles.getLastProp<8>()[0] = 1.0 ;
        Particles.getLastProp<8>()[1] = acos(1);
        Particles.getLastProp<8>()[2] = std::atan2(0.0,0.0);
        Particles.add();
        Particles.getLastPos()[0] = 0;
        Particles.getLastPos()[1] = 0;
        Particles.getLastPos()[2] = -1;
        Particles.getLastProp<8>()[0] = 1.0 ;
        Particles.getLastProp<8>()[1] = acos(-1);
        Particles.getLastProp<8>()[2] = std::atan2(0.0,0.0);
/*        int n1=34,n2=25;
        for(int i=0;i<=n1;i++)
        {
            for(int j=1;j<=n2;j++)
            {
                double spacing1=2 * M_PI/n1;
                double spacing2=M_PI/n2;
                double theta,phi;
                phi = (2 * M_PI - spacing2) * i* spacing2;
                theta = (M_PI - spacing1) * j * spacing1;
                double x = R * sin(theta) * cos(phi);
                double y = R * sin(theta) * sin(phi);
                double z = R * cos(theta);
                Particles.add();
                Particles.getLastPos()[0] = x;
                Particles.getLastProp<8>()[0] = R;
                Particles.getLastPos()[1] = y;
                Particles.getLastProp<8>()[1] = theta;
                Particles.getLastPos()[2] = z;
                Particles.getLastProp<8>()[2] = phi;
            }
        }*/
        Particles.map();
        Particles.ghost_get<0>();


        std::unordered_map<const lm,double,key_hash,key_equal> Vr;
        std::unordered_map<const lm,double,key_hash,key_equal> V1;
        std::unordered_map<const lm,double,key_hash,key_equal> V2;
        Vr[std::make_tuple(0,0)]=0.0;
        Vr[std::make_tuple(1,0)]=0.0;
        Vr[std::make_tuple(1,-1)]=0.0;
        Vr[std::make_tuple(1,1)]=0.0;

        V1[std::make_tuple(0,0)]=0.0;
        V1[std::make_tuple(1,0)]=0.0;
        V1[std::make_tuple(1,-1)]=0.0;
        V1[std::make_tuple(1,1)]=0.0;

        V2[std::make_tuple(0,0)]=0.0;
        V2[std::make_tuple(1,0)]=1.0;

        V2[std::make_tuple(1,-1)]=0.0;
        V2[std::make_tuple(1,1)]=0.0;

        //Particles.ghost_get<0,1,2,3,4,5>();
        auto it2 = Particles.getDomainIterator();
        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Point<3, double> xP = Particles.getProp<8>(p);
            Particles.getProp<0>(p) =0;
            if (xP[0]==1.0) {
                Surface.add();
                Surface.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) =  0;
                //if(xp[2]>0 && xp[2]<spacing) {
                //Particles.getProp<1>(p)[0] = boost::math::spherical_harmonic_r(0, 0,xP[1] ,xP[2]);
                /* double r,theta,phi;
                 r=Particles.getProp<8>(p)[0];
                 theta=Particles.getProp<8>(p)[1];
                 phi=Particles.getProp<8>(p)[2];*/
                //std::cout<<xP[0]<<","<<xP[1]<<","<<xP[2]<<std::endl;
                std::vector<double> SVel=boost::math::sumY<1>(xP[0],xP[1],xP[2],Vr,V1,V2);
                Particles.getProp<1>(p)[0] = SVel[0];//(pow(xp[0],2)+pow(xp[1],2));
                Particles.getProp<1>(p)[1] = SVel[1];//(pow(xp[0],2)+pow(xp[1],2));
                Particles.getProp<1>(p)[2] = SVel[2];
                Particles.getProp<9>(p)[0] = SVel[0];//(pow(xp[0],2)+pow(xp[1],2));
                Particles.getProp<9>(p)[1] = SVel[1];//(pow(xp[0],2)+pow(xp[1],2));
                Particles.getProp<9>(p)[2] = SVel[2];
                //Particles.getProp<8>(p)[0] =  -xp[1];//(pow(xp[0],2)+pow(xp[1],2));
                //Particles.getProp<8>(p)[1] =  xp[0];//(pow(xp[0],2)+pow(xp[1],2));
                //Particles.getProp<8>(p)[2] =  0;
                //}


            }
            else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                Particles.getProp<0>(p) =  0;
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            ++it2;
        }
        Particles.write("SphereInit.vtk");

        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            Point<3, double> xp = Particles.getPos(p);
            Point<3, double> xP = Particles.getProp<8>(p);

/*            if(xp[0]==0 &&xp[1]==0 &&xp[2]==0){
                std::cout<<"Debug"<<std::endl;
            }*/

            std::unordered_map<const lm,double,key_hash,key_equal> Ur;
            std::unordered_map<const lm,double,key_hash,key_equal> U2;
            std::unordered_map<const lm,double,key_hash,key_equal> U1;
            std::unordered_map<const lm,double,key_hash,key_equal> Plm;

            for (int l = 0; l <= 1; l++) {
                for (int m = -l; m <= l; m++) {
                    auto Er= Vr.find(std::make_tuple(l,m));
                    auto E1= V1.find(std::make_tuple(l,m));
                    auto E2= V2.find(std::make_tuple(l,m));
                    std::vector<double> Sol=boost::math::sph_anasol_u(nu,l,m,Er->second,E1->second,E2->second,xP[0]);
                    Ur[std::make_tuple(l,m)]=Sol[0];
                    U1[std::make_tuple(l,m)]=Sol[1];
                    U2[std::make_tuple(l,m)]=Sol[2];
                    Plm[std::make_tuple(l,m)]=Sol[3];
/*                    if(xp[0]==0 &&xp[1]==0 &&xp[2]==0){
                        std::cout<<l<<","<<m<<":"<<Sol[0]<<","<<Sol[1]<<","<<Sol[2]<<","<<Sol[3]<<","<<std::endl;
                    }*/

                }
            }
            std::vector<double> SVel=boost::math::sumY<1>(xP[0],xP[1],xP[2],Ur,U1,U2);
            Particles.getProp<9>(p)[0] = SVel[0];
            Particles.getProp<9>(p)[1] = SVel[1];
            Particles.getProp<9>(p)[2] = SVel[2];
            Particles.getProp<5>(p) =boost::math::sumY_Scalar<1>(xP[0],xP[1],xP[2],Plm);
        }

        vector_dist_subset<3, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>,double,double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>>> Particles_bulk(Particles,bulk);
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

        std::cout << "Sphere Init" << std::endl;
        V_t=V;
        V_B=V_anal;
        P=0;
        P_bulk=0;
        Particles.write("Sphere.vtk");
        eq_id vx,vy,vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);

        double sampling=1.9;
        double sampling2=1.9;
        double rCut2=3.1*spacing;

        Derivative_x Dx(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dx(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_y Dy(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dy(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,sampling, support_options::RADIUS),B_Dz(Particles_bulk, 2, rCut,sampling, support_options::RADIUS);
/*        Derivative_x Dx(Particles, 2, rCut,sampling),B_Dx(Particles_bulk, 2, rCut,sampling);
        Derivative_y Dy(Particles, 2, rCut,sampling),B_Dy(Particles_bulk, 2, rCut,sampling);
        Derivative_z Dz(Particles, 2, rCut,sampling),B_Dz(Particles_bulk, 2, rCut,sampling);*/
        //Laplacian Lap(Particles, 2, rCut2,sampling2);//, support_options::RADIUS);
        Derivative_xx Dxx(Particles, 2, rCut2,sampling2);//,support_options::RADIUS);
        Derivative_yy Dyy(Particles, 2, rCut2,sampling2);//,support_options::RADIUS);
        Derivative_zz Dzz(Particles, 2, rCut2,sampling2);//,support_options::RADIUS);

/*        std::cout << "Dx" << std::endl;
        Dx.checkMomenta(Particles);
        std::cout << "Dy" << std::endl;
        Dy.checkMomenta(Particles);
        std::cout << "Dz" << std::endl;
        Dz.checkMomenta(Particles);
        std::cout << "Lap" << std::endl;
        Lap.checkMomenta(Particles);
        Particles.write("Sphere.vtk");*/
        //Dx.DrawKernel<0>(Particles,25828);
        //B_Dx.checkMomenta(Particles_subset);
        //Particles.write("Sphere.vtk");
        std::cout << "DCPSE KERNELS DONE" << std::endl;
        petsc_solver<double> solverPetsc;
        //solverPetsc.setRestart(250);
        //Particles.ghost_get<0,1,2,3,4,5,6,7>();
        timer tt;
        double sum=0,sum1=0;
        V_t=V;
        double V_err_eps = 1e-5;
        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 60;
        int ctr = 0, errctr, Vreset = 0;
        V_err = 1;
        n = 0;
        tt.start();
        while (V_err >= V_err_eps && n <= nmax) {
            Particles.write_frame("StokesSphere",n);
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
            std::cout << "Stokes Solved" << std::endl;
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
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
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
            //Particles.write_frame("V_debug", n);
            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;
            }
            double worst=0;
            double L2=0;
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                double dx=Particles.getProp<2>(p)[0] - Particles.getProp<1>(p)[0];
                double dy=Particles.getProp<2>(p)[1] - Particles.getProp<1>(p)[1];
                double dz=Particles.getProp<2>(p)[2] - Particles.getProp<1>(p)[2];
                L2 += dx*dx+dy*dy+dz*dz;
                if (fabs(dx*dx+dy*dy+dz*dz) > worst) {
                    worst = fabs(dx*dx+dy*dy+dz*dz);
                }
            }
            v_cl.sum(worst);
            v_cl.sum(L2);
            v_cl.execute();

                std::cout << "L2: " << sqrt(L2 / bulk.size())
                          << std::endl;
                std::cout << "L_inf: " << worst
                          << std::endl;

        }
        tt.stop();
        if (v_cl.rank() == 0) {
            std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds."
                      << std::endl;
        }
        Particles.write("StokesSphere");

        double worst=0;
        double L2=0;
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            double dx=Particles.getProp<2>(p)[0] - Particles.getProp<1>(p)[0];
            double dy=Particles.getProp<2>(p)[1] - Particles.getProp<1>(p)[1];
            double dz=Particles.getProp<2>(p)[2] - Particles.getProp<1>(p)[2];
            L2 += dx*dx+dy*dy+dz*dz;
            if (fabs(dx*dx+dy*dy+dz*dz) > worst) {
                worst = fabs(dx*dx+dy*dy+dz*dz);
            }
        }

        v_cl.sum(worst);
        v_cl.sum(L2);
        v_cl.execute();

        std::cout.precision(17);

        if (v_cl.rank() == 0) {
            std::cout<<"Bulk Size: "<<grd_sz<<","<<bulk.size()<<std::endl;
        std::cout << "L2_Final: " << sqrt(L2/bulk.size())
                  << std::endl;
        std::cout << "L_inf_Final: " << worst
                  << std::endl;
            }

    }



BOOST_AUTO_TEST_SUITE_END()

