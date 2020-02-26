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


//! Specify the general characteristic of system to solve
struct equations2 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
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

const bool equations2::boundary[] = {NON_PERIODIC, NON_PERIODIC};

//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests2)

    BOOST_AUTO_TEST_CASE(dcpse_Lid2) {
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
        openfpm::vector<aggregate<int>> ref_p;
        openfpm::vector<aggregate<int>> ref_p2;
        openfpm::vector<aggregate<int>> FullBulk;


        auto P = getV<0>(Particles);
        auto V = getV<1>(Particles);
        auto dV = getV<2>(Particles);
        auto RHS = getV<3>(Particles);
        auto Vtemp = getV<4>(Particles);
        auto divV = getV<5>(Particles);


        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) - spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
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
            if (xp[0]!=0 || xp[1]!=0) {
                if (xp[0]!=0+spacing || xp[1]!=0) {
                    FullBulk.add();
                    FullBulk.last().get<0>() = p.getKey();
                }
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
                    Particles.getProp<1>(p)[0] =  0;
                    Particles.getProp<1>(p)[1] =  0;
                }
                else if (xp[0]==0+spacing && xp[1]==0) {
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

        Derivative_x Dx(Particles, 2, rCut,2);
        Derivative_y Dy(Particles, 2, rCut,2);
        Gradient Grad(Particles, 2, rCut,2);
        Laplacian Lap(Particles, 2, rCut, 2);
        Advection Adv(Particles, 2, rCut, 2);
        Divergence Div(Particles, 2, rCut, 2);
        double dt=5e-4;
        int n=30;
        double nu=1e-2;
/*        divV=Div(V);
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
        return;*/
        Particles.write_frame("Re1000-1e-4-Lid",0);
        for(int i =1; i<=n ;i++)
        {   dV=dt*(nu*Lap(V)-Adv(V,V));
            RHS=1/dt*Div(dV);
            std::cout<<"RHS Done"<<std::endl;
            DCPSE_scheme<equations2,decltype(Particles)> Solver(2*rCut, Particles);
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
            Vtemp=V+dV-dt*Grad(P);
            V=Vtemp;
            divV=Div(V);
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


BOOST_AUTO_TEST_SUITE_END()


