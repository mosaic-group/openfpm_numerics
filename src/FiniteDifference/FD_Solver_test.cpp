

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "FD_Solver.hpp"
#include "Solvers/petsc_solver.hpp"
#include "Solvers/umfpack_solver.hpp"
#include "FD_expressions.hpp"
#include "FD_op.hpp"
#include "Grid/staggered_dist_grid.hpp"


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
        FD::Lap Lap;
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

    BOOST_AUTO_TEST_CASE(solver_Lap_stag)
    {
        const size_t sz[2] = {82,82};
        Box<2, double> box({0, 0}, {1, 1});
        periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};
        Ghost<2,long int> ghost(1);

        staggered_grid_dist<2, double, aggregate<double,double,double>> domain(sz, box, ghost, bc);


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


        FD_scheme<equations2d1_stag,decltype(domain)> Solver(ghost,domain);
        FD::Lap Lap;

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

        auto it2 = domain.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();

            if (fabs(domain.getProp<0>(p) - domain.getProp<2>(p)) > worst) {
                worst = fabs(domain.getProp<0>(p) - domain.getProp<2>(p));
            }

            ++it2;
        }

        std::cout << "Maximum Error: " << worst << std::endl;
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

    struct lid_nn_stag
    {
    	// dimensionaly of the equation (2D problem 3D problem ...)
    	static const unsigned int dims = 2;

    	// number of fields in the system v_x, v_y, P so a total of 3
    	static const unsigned int nvar = 3;

    	// boundary conditions PERIODIC OR NON_PERIODIC
    	static const bool boundary[];

    	// type of space float, double, ...
    	typedef float stype;

    	// type of base grid, it is the distributed grid that will store the result
    	// Note the first property is a 2D vector (velocity), the second is a scalar (Pressure)
    	typedef staggered_grid_dist<2,float,aggregate<float[2],float>,CartDecomposition<2,float>> b_grid;

    	// type of SparseMatrix, for the linear system, this parameter is bounded by the solver
    	// that you are using, in case of umfpack it is the only possible choice
    	typedef SparseMatrix<double,int> SparseMatrix_type;

    	// type of Vector for the linear system, this parameter is bounded by the solver
    	// that you are using, in case of umfpack it is the only possible choice
    	typedef Vector<double> Vector_type;

        typedef umfpack_solver<double> solver_type;
    };

    const bool lid_nn_stag::boundary[] = {NON_PERIODIC,NON_PERIODIC};

    BOOST_AUTO_TEST_CASE(solver_Lid_driven_cavity_stag)
    {
    	Vcluster<> & v_cl = create_vcluster();

    	if (v_cl.getProcessingUnits() > 3)
    		return;

    	//! [lid-driven cavity 2D]

    	// velocity in the grid is the property 0, pressure is the property 1
    	constexpr int velocity = 0;
    	constexpr int pressure = 1;

    	constexpr int x = 0;
    	constexpr int y = 1;

    	// Domain, a rectangle
    	Box<2,float> domain({0.0,0.0},{3.0,1.0});

    	// Ghost (Not important in this case but required)
    	Ghost<2,float> g(0.011); // suited for particles
        //Ghost<2,long int> g(1); //better for grids
    	// Grid points on x=256 and y=64
    	long int sz[] = {256,64};
    	size_t szu[2];
    	szu[0] = (size_t)sz[0];
    	szu[1] = (size_t)sz[1];

    	// We need one more point on the left and down part of the domain
    	// This is given by the boundary conditions that we impose, the
    	// reason is mathematical in order to have a well defined system
    	// and cannot be discussed here
    	Padding<2> pd({1,1},{0,0});

    	// Distributed grid that store the solution
    	staggered_grid_dist<2,float,aggregate<Point<2,float>,float,Point<2,float>,float>> g_dist(szu,domain,g);
    	grid_dist_id<2,float,aggregate<Point<2,float>,float,Point<2,float>,float>> g_dist_normal(g_dist.getDecomposition(),szu,g);

    	openfpm::vector<comb<2>> cmb_v;
    	cmb_v.add({0,-1});
    	cmb_v.add({-1,0});

    	g_dist.setDefaultStagPosition();
    	g_dist.setStagPosition<0>(cmb_v);

    	// It is the maximum extension of the stencil
    	Ghost<2,long int> stencil_max(1);

    	// Finite difference scheme
    	FD_scheme<lid_nn_stag,decltype(g_dist)> fd(pd,stencil_max,g_dist);

        auto P =  FD::getV_stag<1>(g_dist);
        auto v = FD::getV_stag<0>(g_dist);
        auto RHS = FD::getV_stag<2>(g_dist);
        //auto RHSy = FD::getV_stag<3>(g_dist);


        v.setVarId(0);
        P.setVarId(2);

        FD::Derivative_x_stag Dx;
        FD::Derivative_y_stag Dy;
        FD::Lap Lap;

        double nu = 1.0;

        auto Stokes_vx = nu * Lap(v[x]) - Dx(P);
        auto Stokes_vy = nu * Lap(v[y]) - Dy(P);

        auto incompressibility = Dx(v[x]) + Dy(v[y]);

        eq_id ic,vx,vy;

        ic.setId(2);
        vx.setId(0);
        vy.setId(1);

        comb<2> center_cell({0,0});
        comb<2> left_cell({0,-1});
        comb<2> bottom_cell({-1,0});
        comb<2> corner_right({-1,1});
        comb<2> corner_dw({-1,-1});
        comb<2> corner_up({1,-1});

        auto it = g_dist.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            auto gkey = it.getGKey(key);
            double xg = gkey.get(0) * g_dist.spacing(0);
            double yg = gkey.get(1) * g_dist.spacing(1);
            if (xg==3.0){
                //v[y].value(key,bottom_cell)=1.0;
                g_dist.getProp<2>(key)[1] = 1.0;
                g_dist.getProp<3>(key) = 1.0;
                //std::cout<<it.getGKey(key).get(0)<<","<<it.getGKey(key).get(1)<<":"<<g_dist.getProp<3>(key)<<std::endl;
            }
            ++it;
        }

    	// Here we impose the equation, we start from the incompressibility Eq imposed in the bulk with the
    	// exception of the first point {0,0} and than we set P = 0 in {0,0}, why we are doing this is again
    	// mathematical to have a well defined system, an intuitive explanation is that P and P + c are both
    	// solution for the incompressibility equation, this produce an ill-posed problem to make it well posed
    	// we set one point in this case {0,0} the pressure to a fixed constant for convenience P = 0

        fd.impose(incompressibility, {0,0},{sz[0]-2,sz[1]-2}, 0.0,ic,true);
    	fd.impose(P, {0,0},{0,0},0.0,ic);

    	// Here we impose the Eq1 and Eq2
    	fd.impose(Stokes_vx, {1,0},{sz[0]-2,sz[1]-2},0.0,vx,left_cell);
    	fd.impose(Stokes_vy, {0,1},{sz[0]-2,sz[1]-2},0.0,vy,bottom_cell);

    	// Staggering pattern in the domain
    	// 		+--Vy-+
        //		|     |
        //	   Vx  P  Vx
        //		|     |
        //		0--Vy-+
        //
    	// Imposing B1
    	// Left Wall Vx without touching top wall and with specified location in the cell ("left_cell").
    	fd.impose(v[x], {0,0},{0,sz[1]-2},0.0,vx,left_cell);
        // Left Wall Vy without touching top wall
        // (Starts at -1 index since Vy is on top of the cell and needs special boundary treatment "corner_dw")
        //  cells in the Domain near 0,0 look like:
        //
        //		       |     |
        //	           Vx  P(0,0)
        //		       |     |
        // 	           Vy---Vy--+
        //             :
        //         "corner_right of -1" specified by starting at -1,0
    	fd.impose(v[y], {-1,0},{-1,sz[1]-1},0.0,vy,corner_right);

        //Imposing B2
    	// Similarly Right Wall (Also need "corner_dw" treatment for Vy, hence +1 in the y index)
    	fd.impose(v[x],{sz[0]-1,0},{sz[0]-1,sz[1]-2},0.0,vx,left_cell);
    	fd.impose(v[y],{sz[0]-1,0},{sz[0]-1,sz[1]-1},prop_id<3>(),vy,corner_dw);

    	// Imposing B3
        // Similarly Bottom Wall (needs "corner_up" treatment for Vx, hence -1 in the x index)
    	fd.impose(v[x], {0,-1},{sz[0]-1,-1},0.0,vx,corner_up);
    	fd.impose(v[y], {0,0},{sz[0]-2,0},0.0,vy,bottom_cell);
    	// Imposing B4
        // Similarly Top Wall (needs "corner_dw" treatment for Vx, hence +1 in the x index)
        fd.impose(v[x],{0,sz[1]-1},{sz[0]-1,sz[1]-1},0.0,vx,corner_dw);
    	fd.impose(v[y],{0,sz[1]-1},{sz[0]-2,sz[1]-1},0.0,vy,bottom_cell);

    	// When we pad the grid, there are points of the grid that are not
    	// touched by the previous condition. Mathematically this lead
    	// to have too many variables for the conditions that we are imposing.
    	// Here we are imposing variables that we do not touch to zero
    	//

    	// Padding pressure
    	fd.impose(P, {-1,-1},{sz[0]-1,-1},0.0,ic);
    	fd.impose(P, {-1,sz[1]-1},{sz[0]-1,sz[1]-1},0.0,ic);
    	fd.impose(P, {-1,0},{-1,sz[1]-2},0.0,ic);
    	fd.impose(P, {sz[0]-1,0},{sz[0]-1,sz[1]-2},0.0,ic);

    	// Impose v_x Padding Impose v_y padding
    	fd.impose(v[x], {-1,-1},{-1,sz[1]-1},0.0,vx,left_cell);
    	fd.impose(v[y], {-1,-1},{sz[0]-1,-1},0.0,vy,bottom_cell);

    	fd.solve(v[x],v[y],P);


    	auto it2 = g_dist_normal.getDomainIterator();

    	while (it2.isNext())
    	{
    		auto key = it2.get();

    		g_dist_normal.template getProp<0>(key)[0] = v[x].value(key,left_cell);
    		g_dist_normal.template getProp<0>(key)[1] = v[y].value(key,bottom_cell);

    		g_dist_normal.template getProp<1>(key) = P.value(key,center_cell);
            g_dist_normal.template getProp<2>(key)[x] = RHS[x].value(key,center_cell);
            //SOME PROBLEM IN THIS LINE:
            return;
            g_dist_normal.template getProp<2>(key)[y] = RHS[y].value(key,center_cell);


            ++it2;
    	}

    	g_dist_normal.write("out_test");

    	//! [Copy the solution to grid]

    	//g_dist.write("out_test")
    }


    struct ana_nn_stag
    {
        // dimensionaly of the equation (2D problem 3D problem ...)
        static const unsigned int dims = 2;

        // number of fields in the system v_x, v_y, P so a total of 3
        static const unsigned int nvar = 3;

        // boundary conditions PERIODIC OR NON_PERIODIC
        static const bool boundary[];

        // type of space float, double, ...
        typedef float stype;

        // type of base grid, it is the distributed grid that will store the result
        // Note the first property is a 2D vector (velocity), the second is a scalar (Pressure)
        typedef staggered_grid_dist<2,float,aggregate<float[2],float,float[2],float>,CartDecomposition<2,float>> b_grid;

        // type of SparseMatrix, for the linear system, this parameter is bounded by the solver
        // that you are using, in case of umfpack it is the only possible choice
        // typedef SparseMatrix<double,int> SparseMatrix_type;
        // // typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

        // // type of Vector for the linear system, this parameter is bounded by the solver
        // // that you are using, in case of umfpack it is the only possible choice
        // typedef Vector<double> Vector_type;
        // // typedef Vector<double, PETSC_BASE> Vector_type;

        // typedef umfpack_solver<double> solver_type;
        // // typedef petsc_solver<double> solver_type;
        //! type of SparseMatrix for the linear solver
        typedef SparseMatrix<double, int, PETSC_BASE> SparseMatrix_type;

        //! type of Vector for the linear solver
        typedef Vector<double, PETSC_BASE> Vector_type;

        typedef petsc_solver<double> solver_type;
    };

    const bool ana_nn_stag::boundary[] = {NON_PERIODIC,NON_PERIODIC};

    /*
    In 2D we use exact solution:

     u = x^2 + y^2
     v = 2 x^2 - 2xy
     p = x + y - 1
     f_x = f_y = 3

    so that

     -\Delta u + \nabla p + f = <-4, -4> + <1, 1> + <3, 3> = 0
     \Delta u - \nabla p  = +f
     \nabla \cdot u           = 2x - 2x                    = 0
    */

    BOOST_AUTO_TEST_CASE(solver_ana_stag)
    {
        Vcluster<> & v_cl = create_vcluster();

        if (v_cl.getProcessingUnits() > 3)
            return;

        // velocity in the grid is the property 0, pressure is the property 1
        constexpr int velocity = 0;
        constexpr int pressure = 1;

        constexpr int x = 0;
        constexpr int y = 1;

        // Domain, a rectangle
        Box<2,float> domain({0.0,0.0},{1.0,1.0});

        // Grid points on x=256 and y=256
        long int sz[] = {256,256};
        // Ghost (Not important in this case but required)
        Ghost<2,float> g(domain.getHigh(0)/sz[0]);

        size_t szu[2];
        szu[0] = (size_t)sz[0];
        szu[1] = (size_t)sz[1];

        Padding<2> pd({1,1},{0,0});

        // Distributed grid that store the solution
        staggered_grid_dist<2,float,aggregate<Point<2,float>,float,Point<2,float>,float,float,float>> g_dist(szu,domain,g);
        grid_dist_id<2,float,aggregate<Point<2,float>,float,Point<2,float>,float,float,float>> g_dist_normal(g_dist.getDecomposition(),szu,g);

        double hx = g_dist.spacing(0);
        double hy = g_dist.spacing(1);

        openfpm::vector<comb<2>> cmb_v;
        cmb_v.add({0,-1});
        cmb_v.add({-1,0});

        g_dist.setDefaultStagPosition();
        g_dist.setStagPosition<0>(cmb_v);

        // It is the maximum extension of the stencil
        Ghost<2,long int> stencil_max(1);

        // Finite difference scheme
        FD_scheme<lid_nn_stag,decltype(g_dist)> fd(pd,stencil_max,g_dist);

        auto P =  FD::getV_stag<1>(g_dist);
        auto v = FD::getV_stag<0>(g_dist);
        auto Ana_v = FD::getV_stag<2>(g_dist);
        auto Ana_P = FD::getV_stag<3>(g_dist);

        v.setVarId(0);
        P.setVarId(2);

        FD::Derivative_x_stag Dx;
        FD::Derivative_y_stag Dy;
        FD::Lap Lap;

        auto Stokes_vx = Lap(v[x]) - Dx(P);
        auto Stokes_vy = Lap(v[y]) - Dy(P);

        auto incompressibility = Dx(v[x]) + Dy(v[y]);

        eq_id ic,vx,vy;

        ic.setId(2);
        vx.setId(0);
        vy.setId(1);

        comb<2> center_cell({0,0});
        comb<2> left_cell({0,-1});
        comb<2> bottom_cell({-1,0});
        comb<2> corner_right({-1,1});
        comb<2> corner_dw({-1,-1});
        comb<2> corner_up({1,-1});

        auto fu = [](double x, double y) { return x*x + y*y; };
        auto fv = [](double x, double y) { return 2.0*x*x - 2.0*x*y; };

        // Staggering pattern in the domain
        //      +--Vy-+
        //      |     |
        //     Vx  P  Vx
        //      |     |
        //      0--Vy-+
        //

        fd.impose(incompressibility, {0,0},{sz[0]-2,sz[1]-2}, 0.0,ic,true);
        fd.impose(P, {0,0},{0,0},(hx+hy)/2-1.0,ic);

        fd.impose(Stokes_vx, {1,0},{sz[0]-2,sz[1]-2},3.0,vx,left_cell);
        fd.impose(Stokes_vy, {0,1},{sz[0]-2,sz[1]-2},3.0,vy,bottom_cell);

        // Imposing B1 Left Wall
        fd.impose(v[x], {0,0},{0,sz[1]-2},fu,vx,left_cell);
        fd.impose(v[y], {-1,0},{-1,sz[1]-1},fv,vy,corner_right);

        // //Imposing B2 Right Wall
        fd.impose(v[x],{sz[0]-1,0},{sz[0]-1,sz[1]-2},fu,vx,left_cell);
        fd.impose(v[y],{sz[0]-1,0},{sz[0]-1,sz[1]-1},fv,vy,corner_dw);

        // Imposing B3 Bottom Wall
        fd.impose(v[x],{0,-1},{sz[0]-1,-1},fu,vx,corner_up);
        fd.impose(v[y], {0,0},{sz[0]-2,0},fv,vy,bottom_cell);

        // Imposing B4 Top Wall
        fd.impose(v[x],{0,sz[1]-1},{sz[0]-1,sz[1]-1},fu,vx,corner_dw);
        fd.impose(v[y],{0,sz[1]-1},{sz[0]-2,sz[1]-1},fv,vy,bottom_cell);

        // Padding pressure
        fd.impose(P, {-1,-1},{sz[0]-1,-1},0.0,ic);
        fd.impose(P, {-1,sz[1]-1},{sz[0]-1,sz[1]-1},0.0,ic);
        fd.impose(P, {-1,0},{-1,sz[1]-2},0.0,ic);
        fd.impose(P, {sz[0]-1,0},{sz[0]-1,sz[1]-2},0.0,ic);

        // Impose v_x Padding Impose v_y padding
        fd.impose(v[x], {-1,-1},{-1,sz[1]-1},0.0,vx,left_cell);
        fd.impose(v[y], {-1,-1},{sz[0]-1,-1},0.0,vy,bottom_cell);

        fd.solve(v[x],v[y],P);

        auto it2 = g_dist_normal.getDomainIterator();

        while (it2.isNext())
        {
            auto key = it2.get();
            auto gkey = it2.getGKey(key);
            double x = gkey.get(0) * g_dist_normal.spacing(0);
            double y = gkey.get(1) * g_dist_normal.spacing(1);

            g_dist_normal.template getProp<0>(key)[0] = v[0].value(key,left_cell);
            g_dist_normal.template getProp<0>(key)[1] = v[1].value(key,bottom_cell);

            g_dist_normal.template getProp<1>(key) = P.value(key,center_cell);

            g_dist_normal.template getProp<2>(key)[0] = x*x + y*y;
            g_dist_normal.template getProp<2>(key)[1] = 2*x*x - 2*x*y;
            g_dist_normal.template getProp<3>(key) = x + y - 1;
            ++it2;
        }

        g_dist_normal.write("ana_stokes");

        //! [Copy the solution to grid]

        //g_dist.write("ana_stokes")
    }


BOOST_AUTO_TEST_SUITE_END()
