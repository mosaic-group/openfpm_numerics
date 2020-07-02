

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "FD_Solver.hpp"
#include "Solvers/petsc_solver.hpp"
#include "FD_expressions.hpp"
#include "FD_op.hpp"


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

BOOST_AUTO_TEST_SUITE( FD_Solver_test )


BOOST_AUTO_TEST_CASE(solver_check_diagonal)
{
	const size_t sz[2] = {81,81};
    Box<2, double> box({0, 0}, {1, 1});
    periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};
    Ghost<2,long int> ghost(1);

    grid_dist_id<2, double, aggregate<double,double,double>> domain(sz, box, ghost, bc);


    auto it = domain.getDomainIterator();
    while (it.isNext())
    {
    	auto key = it.get();
    	auto gkey = it.getGKey(key);
        double x = gkey.get(0) * domain.spacing(0);
        double y = gkey.get(1) * domain.spacing(1);
        domain.get<1>(key) = x+y;

        ++it;
    }

    domain.ghost_get<0>();

    auto v =  FD::getV<0>(domain);
    auto RHS= FD::getV<1>(domain);
    auto sol= FD::getV<2>(domain);

    FD_scheme<equations2d1,decltype(domain)> Solver(ghost,domain);

    Solver.impose(5.0*v,{0,0},{80,80}, prop_id<1>());
    Solver.solve(sol);


    domain.write("basic_test");
}


    BOOST_AUTO_TEST_CASE(solver_Lap)
    {
        const size_t sz[2] = {82,82};
        Box<2, double> box({0, 0}, {1, 1});
        periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};
        Ghost<2,long int> ghost(1);

        grid_dist_id<2, double, aggregate<double,double,double>> domain(sz, box, ghost, bc);


        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            auto gkey = it.getGKey(key);
            double x = gkey.get(0) * domain.spacing(0);
            double y = gkey.get(1) * domain.spacing(1);
            domain.get<0>(key) = sin(M_PI*x)*sin(M_PI*y);
            domain.get<1>(key) = -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
            ++it;
        }

        domain.ghost_get<0>();
        FD::Derivative_x Dx;
        FD::Derivative_y Dy;
        auto v =  FD::getV<0>(domain);
        auto RHS= FD::getV<1>(domain);
        auto sol= FD::getV<2>(domain);


        FD_scheme<equations2d1,decltype(domain)> Solver(ghost,domain);
        FD::Laplacian_xy Lap;
        FD::L2Error L2Error;
        FD::LInfError LInfError;

/*        Solver.impose(Lap(v),{1,1},{79,79}, prop_id<1>());
        Solver.impose(v,{0,0},{80,0}, prop_id<0>());
        Solver.impose(v,{0,1},{0,79}, prop_id<0>());
        Solver.impose(v,{0,80},{80,80}, prop_id<0>());
        Solver.impose(v,{80,1},{80,79}, prop_id<0>());
        Solver.solve(sol);*/

        Solver.impose(Lap(v),{1,1},{80,80}, prop_id<1>());
        Solver.impose(v,{0,0},{81,0}, prop_id<0>());
        Solver.impose(v,{0,1},{0,80}, prop_id<0>());
        Solver.impose(v,{0,81},{81,81}, prop_id<0>());
        Solver.impose(v,{81,1},{81,80}, prop_id<0>());
        /*auto A=Solver.getA();
        A.write("Lap_Matrix");*/

        Solver.solve(sol);
        auto l2error = L2Error(v, sol);
        auto linferror = LInfError(v, sol);

        std::cout << "L2 error: " << l2error << std::endl;
        std::cout << "L∞ Error: " << linferror << std::endl;

        domain.write("FDSOLVER_Lap_test");
    }


    /*
In 3D we use exact solution:

Vx = x^2 + y^2
Vy = y^2 + z^2
Vz = x^2 + y^2 - 2(x+y)z
p = x + y + z - 3/2
f_x = f_y = f_z = 3

-\Delta u + \nabla p + f = <-4, -4, -4> + <1, 1, 1> + <3, 3, 3> = 0
\nabla \cdot u           = 2x + 2y - 2(x + y)                   = 0
*/

BOOST_AUTO_TEST_SUITE_END()
