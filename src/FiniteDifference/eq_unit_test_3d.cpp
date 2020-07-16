/*
 * eq_unit_test_3d.hpp
 *
 *  Created on: Jan 4, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_3D_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_3D_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config.h"
#include "Laplacian.hpp"
#include "FiniteDifference/eq.hpp"
#include "FiniteDifference/sum.hpp"
#include "FiniteDifference/mul.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "Vector/Vector.hpp"
#include "Solvers/umfpack_solver.hpp"
#include "data_type/aggregate.hpp"
#include "Solvers/petsc_solver.hpp"
#include "FiniteDifference/FDScheme.hpp"

BOOST_AUTO_TEST_SUITE( eq_test_suite_3d )

//! Specify the general caratteristic of system to solve
struct lid_nn_3d_eigen
{
	//! dimensionaly of the equation ( 3D problem ...)
	static const unsigned int dims = 3;
	//! number of fields in the system
	static const unsigned int nvar = 4;

	//! boundary at X and Y
	static const bool boundary[];

	//! type of space float, double, ...
	typedef float stype;

	//! type of base grid
	typedef grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> b_grid;

	//! type of SparseMatrix for the linear solver
	typedef SparseMatrix<double,int,EIGEN_BASE> SparseMatrix_type;

	//! type of Vector for the linear solver
	typedef Vector<double> Vector_type;

	//! Define the underline grid is staggered
	static const int grid_type = STAGGERED_GRID;
};

struct lid_nn_3d_petsc
{
	//! dimensionaly of the equation ( 3D problem ...)
	static const unsigned int dims = 3;
	//! number of fields in the system
	static const unsigned int nvar = 4;

	//! boundary at X and Y
	static const bool boundary[];

	//! type of space float, double, ...
	typedef float stype;

	//! type of base grid
	typedef grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> b_grid;

	//! type of SparseMatrix for the linear solver
	typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;

	//! type of Vector for the linear solver
	typedef Vector<double,PETSC_BASE> Vector_type;

	//! Define the underline grid is staggered
	static const int grid_type = STAGGERED_GRID;
};

//typedef lid_nn_3d_eigen lid_nn_3d;

const bool lid_nn_3d_eigen::boundary[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
const bool lid_nn_3d_petsc::boundary[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

// Constant Field
struct eta
{
	//! define that eta is a constant field
	typedef void const_field;

	//! therutn the value of the constant
	static float val()	{return 1.0;}
};

template<typename solver_type,typename lid_nn_3d> void lid_driven_cavity_3d()
{
	#include "Equations/stoke_flow_eq_3d.hpp"

	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 3)
		return;

	// Domain
	Box<3,float> domain({0.0,0.0,0.0},{3.0,1.0,1.0});

	// Ghost
	Ghost<3,float> g(0.01);

	long int sz[] = {36,12,12};
	size_t szu[3];
	szu[0] = (size_t)sz[0];
	szu[1] = (size_t)sz[1];
	szu[2] = (size_t)sz[2];

	Padding<3> pd({1,1,1},{0,0,0});

	// velocity in the grid is the property 0, pressure is the property 1
	constexpr int velocity = 0;
	constexpr int pressure = 1;

	// Initialize openfpm
	grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> g_dist(szu,domain,g);

	// Ghost stencil
	Ghost<3,long int> stencil_max(1);

	// Distributed grid
	FDScheme<lid_nn_3d> fd(pd,stencil_max,domain,g_dist);

	// start and end of the bulk

	fd.impose(ic_eq(),0.0, EQ_4, {0,0,0},{sz[0]-2,sz[1]-2,sz[2]-2},true);
	fd.impose(Prs(),  0.0, EQ_4, {0,0,0},{0,0,0});
	fd.impose(vx_eq(),0.0, EQ_1, {1,0},{sz[0]-2,sz[1]-2,sz[2]-2});
	fd.impose(vy_eq(),0.0, EQ_2, {0,1},{sz[0]-2,sz[1]-2,sz[2]-2});
	fd.impose(vz_eq(),0.0, EQ_3, {0,0,1},{sz[0]-2,sz[1]-2,sz[2]-2});

	// v_x
	// R L
	fd.impose(v_x(),0.0, EQ_1, {0,0,0},      {0,sz[1]-2,sz[2]-2});
	fd.impose(v_x(),0.0, EQ_1, {sz[0]-1,0,0},{sz[0]-1,sz[1]-2,sz[2]-2});

	// T B
	fd.impose(avg_y_vx_f(),0.0, EQ_1, {0,-1,0},     {sz[0]-1,-1,sz[2]-2});
	fd.impose(avg_y_vx(),0.0, EQ_1,   {0,sz[1]-1,0},{sz[0]-1,sz[1]-1,sz[2]-2});

	// A F
	fd.impose(avg_z_vx_f(),0.0, EQ_1, {0,-1,-1},     {sz[0]-1,sz[1]-1,-1});
	fd.impose(avg_z_vx(),0.0, EQ_1, {0,-1,sz[2]-1},{sz[0]-1,sz[1]-1,sz[2]-1});

	// v_y
	// R L
	fd.impose(avg_x_vy_f(),0.0, EQ_2,  {-1,0,0},     {-1,sz[1]-1,sz[2]-2});
	fd.impose(avg_x_vy(),1.0, EQ_2,    {sz[0]-1,0,0},{sz[0]-1,sz[1]-1,sz[2]-2});

	// T B
	fd.impose(v_y(), 0.0, EQ_2, {0,0,0},      {sz[0]-2,0,sz[2]-2});
	fd.impose(v_y(), 0.0, EQ_2, {0,sz[1]-1,0},{sz[0]-2,sz[1]-1,sz[2]-2});

	// F A
	fd.impose(avg_z_vy(),0.0, EQ_2,   {-1,0,sz[2]-1}, {sz[0]-1,sz[1]-1,sz[2]-1});
	fd.impose(avg_z_vy_f(),0.0, EQ_2, {-1,0,-1},      {sz[0]-1,sz[1]-1,-1});

	// v_z
	// R L
	fd.impose(avg_x_vz_f(),0.0, EQ_3, {-1,0,0},     {-1,sz[1]-2,sz[2]-1});
	fd.impose(avg_x_vz(),1.0, EQ_3,   {sz[0]-1,0,0},{sz[0]-1,sz[1]-2,sz[2]-1});

	// T B
	fd.impose(avg_y_vz(),0.0, EQ_3, {-1,sz[1]-1,0},{sz[0]-1,sz[1]-1,sz[2]-1});
	fd.impose(avg_y_vz_f(),0.0, EQ_3, {-1,-1,0},   {sz[0]-1,-1,sz[2]-1});

	// F A
	fd.impose(v_z(),0.0, EQ_3, {0,0,0},      {sz[0]-2,sz[1]-2,0});
	fd.impose(v_z(),0.0, EQ_3, {0,0,sz[2]-1},{sz[0]-2,sz[1]-2,sz[2]-1});

	// Padding pressure

	// L R
	fd.impose(Prs(), 0.0, EQ_4, {-1,-1,-1},{-1,sz[1]-1,sz[2]-1});
	fd.impose(Prs(), 0.0, EQ_4, {sz[0]-1,-1,-1},{sz[0]-1,sz[1]-1,sz[2]-1});

	// T B
	fd.impose(Prs(), 0.0, EQ_4, {0,sz[1]-1,-1}, {sz[0]-2,sz[1]-1,sz[2]-1});
	fd.impose(Prs(), 0.0, EQ_4, {0,-1     ,-1}, {sz[0]-2,-1,     sz[2]-1});

	// F A
	fd.impose(Prs(), 0.0, EQ_4, {0,0,sz[2]-1}, {sz[0]-2,sz[1]-2,sz[2]-1});
	fd.impose(Prs(), 0.0, EQ_4, {0,0,-1},      {sz[0]-2,sz[1]-2,-1});

	// Impose v_x  v_y v_z padding
	fd.impose(v_x(), 0.0, EQ_1, {-1,-1,-1},{-1,sz[1]-1,sz[2]-1});
	fd.impose(v_y(), 0.0, EQ_2, {-1,-1,-1},{sz[0]-1,-1,sz[2]-1});
	fd.impose(v_z(), 0.0, EQ_3, {-1,-1,-1},{sz[0]-1,sz[1]-1,-1});

	solver_type solver;
	auto x_ = solver.try_solve(fd.getA(),fd.getB());

	// Bring the solution to grid
	fd.template copy<velocity,pressure>(x_,{0,0},{sz[0]-1,sz[1]-1,sz[2]-1},g_dist);

	std::string s = std::string(demangle(typeid(solver_type).name()));
	s += "_";

	g_dist.write(s + "lid_driven_cavity_3d_p" + std::to_string(v_cl.getProcessingUnits()) + "_grid");

#if !(defined(SE_CLASS3) || defined(TEST_COVERAGE_MODE))

	// Initialize openfpm
	grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> g_dist2(g_dist.getDecomposition(),szu,g);
	g_dist2.load("test/lid_driven_cavity_3d_reference.hdf5");

	auto it2 = g_dist2.getDomainIterator();

	bool test = true;
	while (it2.isNext())
	{
		auto p = it2.get();

		test &= fabs(g_dist2.template getProp<velocity>(p)[0] - g_dist.template getProp<velocity>(p)[0]) < 3.5e-5;
		test &= fabs(g_dist2.template getProp<velocity>(p)[1] - g_dist.template getProp<velocity>(p)[1]) < 3.5e-5;
		test &= fabs(g_dist2.template getProp<velocity>(p)[2] - g_dist.template getProp<velocity>(p)[2]) < 3.5e-5;

		test &= fabs(g_dist2.template getProp<pressure>(p) - g_dist.template getProp<pressure>(p)) < 3.0e-4;

		if (test == false)
		{
			std::cout << g_dist2.template getProp<velocity>(p)[0] << "   " << g_dist.template getProp<velocity>(p)[0] << std::endl;
			std::cout << g_dist2.template getProp<velocity>(p)[1] << "   " << g_dist.template getProp<velocity>(p)[1] << std::endl;
			std::cout << g_dist2.template getProp<velocity>(p)[2] << "   " << g_dist.template getProp<velocity>(p)[2] << std::endl;

			std::cout << g_dist2.template getProp<pressure>(p) << "   " << g_dist.template getProp<pressure>(p) << std::endl;

			break;
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(test,true);

#endif

}

// Lid driven cavity, uncompressible fluid

BOOST_AUTO_TEST_CASE(lid_driven_cavity)
{
#if defined(HAVE_EIGEN) && defined(HAVE_SUITESPARSE)
	lid_driven_cavity_3d<umfpack_solver<double>,lid_nn_3d_eigen>();
#endif
#ifdef HAVE_PETSC
	lid_driven_cavity_3d<petsc_solver<double>,lid_nn_3d_petsc>();
#endif
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_3D_HPP_ */
