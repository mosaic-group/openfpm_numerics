/*
 * eq_unit_test.hpp
 *
 *  Created on: Oct 13, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_

#include "Laplacian.hpp"
#include "FiniteDifference/eq.hpp"
#include "FiniteDifference/sum.hpp"
#include "FiniteDifference/mul.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/scalar.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "Vector/Vector.hpp"
#include "Solvers/umfpack_solver.hpp"
#include "data_type/aggregate.hpp"
#include "FiniteDifference/FDScheme.hpp"

BOOST_AUTO_TEST_SUITE( eq_test_suite )

//! [Definition of the system]

struct lid_nn
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
	typedef grid_dist_id<2,float,aggregate<float[2],float>,CartDecomposition<2,float>> b_grid;

	// type of SparseMatrix, for the linear system, this parameter is bounded by the solver
	// that you are using, in case of umfpack it is the only possible choice
	typedef SparseMatrix<double,int> SparseMatrix_type;

	// type of Vector for the linear system, this parameter is bounded by the solver
	// that you are using, in case of umfpack it is the only possible choice
	typedef Vector<double> Vector_type;

	// Define that the underline grid where we discretize the system of equation is staggered
	static const int grid_type = STAGGERED_GRID;
};

const bool lid_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

//! [Definition of the system]

//! [Definition of the equation of the system in the bulk and at the boundary]

// Constant Field
struct eta
{
	typedef void const_field;

	static float val()	{return 1.0;}
};

// Convenient constants
constexpr unsigned int v[] = {0,1};
constexpr unsigned int P = 2;
constexpr unsigned int ic = 2;

// Create field that we have v_x, v_y, P
typedef Field<v[x],lid_nn> v_x;
typedef Field<v[y],lid_nn> v_y;
typedef Field<P,lid_nn> Prs;

// Eq1 V_x

typedef mul<eta,Lap<v_x,lid_nn>,lid_nn> eta_lap_vx;
typedef D<x,Prs,lid_nn> p_x;
typedef minus<p_x,lid_nn> m_p_x;
typedef sum<eta_lap_vx,m_p_x,lid_nn> vx_eq;

// Eq2 V_y

typedef mul<eta,Lap<v_y,lid_nn>,lid_nn> eta_lap_vy;
typedef D<y,Prs,lid_nn> p_y;
typedef minus<p_y,lid_nn> m_p_y;
typedef sum<eta_lap_vy,m_p_y,lid_nn> vy_eq;

// Eq3 Incompressibility

typedef D<x,v_x,lid_nn,FORWARD> dx_vx;
typedef D<y,v_y,lid_nn,FORWARD> dy_vy;
typedef sum<dx_vx,dy_vy,lid_nn> ic_eq;


// Equation for boundary conditions

/* Consider the staggered cell
 *
 	 	 \verbatim

		+--$--+
		|     |
		#  *  #
		|     |
		0--$--+

	  # = velocity(y)
	  $ = velocity(x)
	  * = pressure

		\endverbatim
 *
 *
 * If we want to impose v_y = 0 on 0 we have to interpolate between # of this cell
 * and # of the previous cell on y, (Average) or Avg operator
 *
 */

// Directional Avg
typedef Avg<x,v_y,lid_nn> avg_vy;
typedef Avg<y,v_x,lid_nn> avg_vx;

typedef Avg<x,v_y,lid_nn,FORWARD> avg_vy_f;
typedef Avg<y,v_x,lid_nn,FORWARD> avg_vx_f;

#define EQ_1 0
#define EQ_2 1
#define EQ_3 2

//! [Definition of the equation of the system in the bulk and at the boundary]

template<typename solver_type,typename lid_nn> void lid_driven_cavity_2d()
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 3)
		return;

	//! [lid-driven cavity 2D]

	// velocity in the grid is the property 0, pressure is the property 1
	constexpr int velocity = 0;
	constexpr int pressure = 1;

	// Domain, a rectangle
	Box<2,float> domain({0.0,0.0},{3.0,1.0});

	// Ghost (Not important in this case but required)
	Ghost<2,float> g(0.01);

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
	grid_dist_id<2,float,aggregate<float[2],float>,CartDecomposition<2,float>> g_dist(szu,domain,g);

	// It is the maximum extension of the stencil
	Ghost<2,long int> stencil_max(1);

	// Finite difference scheme
	FDScheme<lid_nn> fd(pd, stencil_max, domain,g_dist);

	// Here we impose the equation, we start from the incompressibility Eq imposed in the bulk with the
	// exception of the first point {0,0} and than we set P = 0 in {0,0}, why we are doing this is again
	// mathematical to have a well defined system, an intuitive explanation is that P and P + c are both
	// solution for the incompressibility equation, this produce an ill-posed problem to make it well posed
	// we set one point in this case {0,0} the pressure to a fixed constant for convenience P = 0
	fd.impose(ic_eq(),0.0, EQ_3, {0,0},{sz[0]-2,sz[1]-2},true);
	fd.impose(Prs(),  0.0, EQ_3, {0,0},{0,0});

	// Here we impose the Eq1 and Eq2
	fd.impose(vx_eq(),0.0, EQ_1, {1,0},{sz[0]-2,sz[1]-2});
	fd.impose(vy_eq(),0.0, EQ_2, {0,1},{sz[0]-2,sz[1]-2});

	// v_x and v_y
	// Imposing B1
	fd.impose(v_x(),0.0, EQ_1, {0,0},{0,sz[1]-2});
	fd.impose(avg_vy_f(),0.0, EQ_2 , {-1,0},{-1,sz[1]-1});
	// Imposing B2
	fd.impose(v_x(),0.0, EQ_1, {sz[0]-1,0},{sz[0]-1,sz[1]-2});
	fd.impose(avg_vy(),1.0, EQ_2,    {sz[0]-1,0},{sz[0]-1,sz[1]-1});

	// Imposing B3
	fd.impose(avg_vx_f(),0.0, EQ_1, {0,-1},{sz[0]-1,-1});
	fd.impose(v_y(), 0.0, EQ_2, {0,0},{sz[0]-2,0});
	// Imposing B4
	fd.impose(avg_vx(),0.0, EQ_1,   {0,sz[1]-1},{sz[0]-1,sz[1]-1});
	fd.impose(v_y(), 0.0, EQ_2, {0,sz[1]-1},{sz[0]-2,sz[1]-1});

	// When we pad the grid, there are points of the grid that are not
	// touched by the previous condition. Mathematically this lead
	// to have too many variables for the conditions that we are imposing.
	// Here we are imposing variables that we do not touch to zero
	//

	// Padding pressure
	fd.impose(Prs(), 0.0, EQ_3, {-1,-1},{sz[0]-1,-1});
	fd.impose(Prs(), 0.0, EQ_3, {-1,sz[1]-1},{sz[0]-1,sz[1]-1});
	fd.impose(Prs(), 0.0, EQ_3, {-1,0},{-1,sz[1]-2});
	fd.impose(Prs(), 0.0, EQ_3, {sz[0]-1,0},{sz[0]-1,sz[1]-2});

	// Impose v_x Padding Impose v_y padding
	fd.impose(v_x(), 0.0, EQ_1, {-1,-1},{-1,sz[1]-1});
	fd.impose(v_y(), 0.0, EQ_2, {-1,-1},{sz[0]-1,-1});

	solver_type solver;
	auto x = solver.solve(fd.getA(),fd.getB());

	//! [lid-driven cavity 2D]

	//! [Copy the solution to grid]

	fd.template copy<velocity,pressure>(x,{0,0},{sz[0]-1,sz[1]-1},g_dist);

	std::string s = std::string(demangle(typeid(solver_type).name()));
	s += "_";

	//! [Copy the solution to grid]

	g_dist.write(s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid");

#ifdef HAVE_OSX

    std::string file1 = std::string("test/") + s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid_" + std::to_string(v_cl.getProcessUnitID()) + "_test_osx.vtk";
    std::string file2 = s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk";

#else

	#if __GNUC__ == 5

		std::string file1 = std::string("test/") + s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid_" + std::to_string(v_cl.getProcessUnitID()) + "_test_GCC5.vtk";
		std::string file2 = s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk";

	#else

		std::string file1 = std::string("test/") + s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid_" + std::to_string(v_cl.getProcessUnitID()) + "_test_GCC4.vtk";
		std::string file2 = s + "lid_driven_cavity_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk";

	#endif


#endif

    std::cout << "File1: " << file1 << std::endl;
    std::cout << "File2: " << file2 << std::endl;

	// Check that match
	bool test = compare(file1,file2);
	BOOST_REQUIRE_EQUAL(test,true);

}

// Lid driven cavity, incompressible fluid

BOOST_AUTO_TEST_CASE(lid_driven_cavity)
{
	lid_driven_cavity_2d<umfpack_solver<double>,lid_nn>();
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_ */
