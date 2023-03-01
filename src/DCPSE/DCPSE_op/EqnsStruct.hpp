//
// Created by Abhinav Singh on 17.04.20.
//
//Example file to be used only when DCPSE_Solver is defined before.
//Create your own structures for your equations after including the solver.
#ifndef OPENFPM_PDATA_EQNSSTRUCT_HPP
#define OPENFPM_PDATA_EQNSSTRUCT_HPP

#include "Solvers/umfpack_solver.hpp"
#include "Solvers/petsc_solver.hpp"

#ifdef HAVE_PETSC
//! Specify the general characteristic of system to solve
struct equations2d1 {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d2 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d1p {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d2p {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};


struct equations2d3p {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d3 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d4 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 4;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d1 {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3Pz {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3Pyz {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, PERIODIC,PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3Pxz {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d1Pz {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};


#ifdef __NVCC__
struct equations2d1_gpu {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d2_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d1p_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d2p_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d3p_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d3_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations2d4_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 4;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d1_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3Pz_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3Pyz_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, PERIODIC,PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d3Pxz_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};

struct equations3d1Pz_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double, PETSC_BASE> Vector_type;

    typedef petsc_solver<double> solver_type;
};
#endif //__NVCC__

#endif //HAVE_PETSC


//! Specify the general characteristic of system to solve
struct equations2d1E {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d2E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d3E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d4E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 4;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};


struct equations2d1pE {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d2pE {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d3pE {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d3E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d1E {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};


struct equations3d3EPxz {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d3EPz {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC,PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};


#ifdef __NVCC__
struct equations2d1E_gpu {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims=2;
    //! number of fields in the system
    static const unsigned int nvar=1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d2E_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d3E_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d4E_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 4;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};


struct equations2d1pE_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d2pE_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 2;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations2d3pE_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 2;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d3E_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ..
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d1E_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 1;

    //! boundary at X and Y
    static constexpr bool boundary[]={NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d3EPxz_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, NON_PERIODIC,PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};

struct equations3d3EPz_gpu {
    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims = 3;
    //! number of fields in the system
    static const unsigned int nvar = 3;

    //! boundary at X and Y
    static constexpr bool boundary[]={PERIODIC, PERIODIC,PERIODIC};

    //! type of space float, double, ...
    typedef double stype;

    //! type of base particles
    typedef vector_dist_gpu<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;

    typedef umfpack_solver<double> solver_type;
};
#endif //__NVCC__

#endif //OPENFPM_PDATA_EQNSSTRUCT_HPP
