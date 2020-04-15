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
#include "DCPSE_op.hpp"
#include "DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"

//int vector_dist_expression_op<void,void,VECT_COPY_N_TO_N>::i = 0;
//int vector_dist_expression_op<void,void,VECT_COPY_1_TO_N>::i = 0;

//! Specify the general characteristic of system to solve
struct equations3d3 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};

struct equations3d {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};

const bool equations3d::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations3d3::boundary[] = {NON_PERIODIC, NON_PERIODIC};


//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests3)

    BOOST_AUTO_TEST_CASE(stokes_3d) {
        size_t grd_sz=8;
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
            Particles.getProp<3>(p) =0;
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
        vy.setId(2);

        Derivative_x Dx(Particles, 2, rCut,1.9,support_options::RADIUS );
        Derivative_y Dy(Particles, 2, rCut,1.9,support_options::RADIUS);
        Derivative_z Dz(Particles, 2, rCut,1.9,support_options::RADIUS);
        Gradient Grad(Particles, 2, rCut,1.9,support_options::RADIUS );
        Laplacian Lap(Particles, 2, rCut,1.9,support_options::RADIUS);
        Advection Adv(Particles, 2, rCut,1.9,support_options::RADIUS);
        Divergence Div(Particles, 2, rCut,1.9,support_options::RADIUS);

        double nu=1e-2;

        double sum=0,sum2=0;
        int n=15;
        Particles.write_frame("Stokes3d",0);
        V_t=V;
        for(int i=1; i<=n ;i++)
        {   RHS2=-Grad(P);
            DCPSE_scheme<equations3d3,decltype(Particles)> Solver( Particles);
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
            //std::cout << "Stokes Solved" << std::endl;
            RHS=Div(V_star);
            DCPSE_scheme<equations3d,decltype(Particles)> SolverH(Particles,options_solver::LAGRANGE_MULTIPLIER);
            auto Helmholtz = Lap(H);
            auto D_x=Dx(H);
            auto D_y=Dy(H);
            auto D_z=Dz(H);
            SolverH.impose(Helmholtz,bulk,prop_id<3>());
            SolverH.impose(D_y, up_p,0);
            SolverH.impose(D_x, r_p, 0);
            SolverH.impose(-D_y, dw_p,0);
            SolverH.impose(-D_x, l_p,0);
            SolverH.impose(D_z, f_p,0);
            SolverH.impose(-D_z, b_p,0);
            SolverH.solve(H);
            //std::cout << "Helmholtz Solved" << std::endl;
            V=V_star-Grad(H);
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
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            for(int j=0;j<b_p.size();j++)
            {   auto p=dw_p.get<0>(j);
                Particles.getProp<1>(p)[0] =  0;
                Particles.getProp<1>(p)[1] =  0;
                Particles.getProp<1>(p)[2] =  0;
            }
            P=P-Lap(H)+0.5*Adv(V_t,H);
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

/*    BOOST_AUTO_TEST_CASE(dcpse_sphere) {
        const size_t sz[3] = {31,31,31};
        Box<3, double> box({0, 0}, {1,1});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.1 * spacing;
        Ghost<2, double> ghost(rCut);
        //                                  P        V                 v_star           RHS            V_t       Helmholtz
        vector_dist<2, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>, double,    double>> Particles(0, box, bc, ghost);
        double r=5,theta=0,phi=-M_PI/2;
        double dtheta=2*M_PI/sz[0],dphi=M_PI/sz[0];
        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = r*cos(theta)*cos(phi);
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = r*sin(theta)*cos(phi);
            double z = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = r*sin(phi);
            theta+=dtheta;
            phi+=dphi;
            ++it;
        }

        Sphere<3, double> bulk({0,0,0},5);
        openfpm::vector<Sphere<3, double>> spheres;
        spheres.add(bulk);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_sph;
        vtk_sph.add(bulk);
        vtk_sph.write("vtk_sph.vtk");

        Particles.map();
        Particles.ghost_get<0>();

        Particles.write_frame("Sphere",i);\
    }
*/
BOOST_AUTO_TEST_SUITE_END()


