//
// Created by Abhinav Singh on 17.04.20.
//

#ifndef OPENFPM_PDATA_EQNSSTRUCT_HPP
#define OPENFPM_PDATA_EQNSSTRUCT_HPP

//! Specify the general characteristic of system to solve
struct equations {

    //! dimensionaly of the equation ( 3D problem ...)
    static const unsigned int dims;
    //! number of fields in the system
    static const unsigned int nvar;

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

struct equations* create_equations( int dims, int nvar, bool boundary[])
{
    struct equations* newEquations = (struct equations*)malloc( sizeof( struct equations ) );
    if( newEquations )
    {
        newEquations->dims = a;
        newEquations->nvar = b;
        newEquations->boundary[] = boundary[];
    }
    return newEquations;
}

void destroy_equations( struct equations* obj )
{
    if( obj )
        free( obj );
}

void print_equations( struct equations* obj )
{
    if( obj )
    {
        std::cout<<("Equation dimension->dims = %d\n"<<obj->dims)<<std::endl;
        printf("Equation variables->vars = %d\n"<<obj->vars)<<std::endl;
        printf("Equation boundary conditions->vars = %d\n"<<obj->boundary[])<<std::endl;
    }
}


/*
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

const bool equations2d1p::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations2d::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations2d1p::boundary[] = {PERIODIC, NON_PERIODIC};
const bool equations2d2p::boundary[] = {PERIODIC, NON_PERIODIC};


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

const bool equations3d1::boundary[] = {NON_PERIODIC, NON_PERIODIC};
const bool equations3d3::boundary[] = {NON_PERIODIC, NON_PERIODIC};
*/



#endif //OPENFPM_PDATA_EQNSSTRUCT_HPP
