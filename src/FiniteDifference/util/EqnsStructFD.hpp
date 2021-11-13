//
// Created by Abhinav Singh on 17.04.20.
//

#ifndef OPENFPM_PDATA_EQNSSTRUCT_HPP
#define OPENFPM_PDATA_EQNSSTRUCT_HPP

#include "Solvers/umfpack_solver.hpp"
#include "Solvers/petsc_solver.hpp"

//! Specify the general characteristic of system to solve
struct equations2d1 {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double,double,double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

const bool equations2d1::boundary[] = {NON_PERIODIC,NON_PERIODIC};

struct equations2d2 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations2d2::boundary[]={NON_PERIODIC,NON_PERIODIC};

struct equations2d1p {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations2d1p::boundary[]={PERIODIC,PERIODIC};

struct equations2d2p {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations2d2p::boundary[]={PERIODIC,PERIODIC};

struct equations2d3p {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations2d3p::boundary[]={PERIODIC,PERIODIC};


struct equations2d3 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations2d3::boundary[]={NON_PERIODIC,NON_PERIODIC};

struct equations2d4 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 4;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations2d4::boundary[]={NON_PERIODIC,NON_PERIODIC};


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
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations3d3::boundary[]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

struct equations3d1 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
const bool equations3d1::boundary[]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};


//! Specify the general characteristic of system to solve
struct equations2d1E {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations2d1E::boundary[]={NON_PERIODIC,NON_PERIODIC};

struct equations2d2E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations2d2E::boundary[]={NON_PERIODIC,NON_PERIODIC};

struct equations2d3E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations2d3E::boundary[]={NON_PERIODIC,NON_PERIODIC};

struct equations2d4E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 4;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations2d4E::boundary[]={NON_PERIODIC,NON_PERIODIC};



struct equations2d1pE {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations2d1pE::boundary[]={PERIODIC,PERIODIC};


struct equations2d2pE {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

const bool equations2d2pE::boundary[]={PERIODIC,PERIODIC};


struct equations2d3pE {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations2d3pE::boundary[]={PERIODIC,PERIODIC,PERIODIC};


struct equations3d3E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

const bool equations3d3E::boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

struct equations3d1E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef grid_dist_id<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
const bool equations3d1E::boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};


//! Specify the general characteristic of system to solve
struct equations2d1_stag {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static const bool boundary[];

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef staggered_grid_dist<dims, double, aggregate<double,double,double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

const bool equations2d1_stag::boundary[] = {NON_PERIODIC,NON_PERIODIC};


#endif //OPENFPM_PDATA_EQNSSTRUCT_HPP
