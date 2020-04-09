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
struct equations3d2 {
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
const bool equations3d2::boundary[] = {NON_PERIODIC, NON_PERIODIC};


//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests3)

/*    BOOST_AUTO_TEST_CASE(dcpse_sphere) {
        const size_t sz[3] = {31,31,31};
        Sphere<3, double> sphere(0, 5);
        double spacing = 0.1;
        Ghost<3, double> ghost(spacing * 3);
        double rCut = 2.5 * spacing;
        //                                  P        V                 v_star           RHS            V_t       Helmholtz
        vector_dist<2, double, aggregate<double,VectorS<3, double>,VectorS<3, double>,double,VectorS<3, double>, double,    double>> Particles(0, Sphere, bc, ghost);
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
            ++it;
        }

        Sphere<3, double> bulk(0,5);

        openfpm::vector<Sphere<3, double>> spheres;
        spheres.add(bulk);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_sph;
        vtk_sph.add(spheres);
        vtk_sph.write("vtk_sph.vtk");

        Particles.map();
        Particles.ghost_get<0>();

        Particles.write_frame("Sphere",i);\
    }*/

BOOST_AUTO_TEST_SUITE_END()


