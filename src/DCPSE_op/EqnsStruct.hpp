//
// Created by Abhinav Singh on 17.04.20.
//

#ifndef OPENFPM_PDATA_EQNSSTRUCT_HPP
#define OPENFPM_PDATA_EQNSSTRUCT_HPP

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
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};

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
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};


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
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};

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
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};

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
    typedef vector_dist<dims, double, aggregate<double>> b_part;

    //! type of SparseMatrix for the linear solver
    typedef SparseMatrix<double, int, EIGEN_BASE> SparseMatrix_type;

    //! type of Vector for the linear solver
    typedef Vector<double> Vector_type;
};


#endif //OPENFPM_PDATA_EQNSSTRUCT_HPP
