/*
 * eq_unit_test_3d.hpp
 *
 *  Created on: Jan 4, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_3D_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_3D_HPP_


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

BOOST_AUTO_TEST_SUITE( eq_test_suite_3d )

//!

struct lid_nn_3d_eigen
{
	// dimensionaly of the equation ( 3D problem ...)
	static const unsigned int dims = 3;
	// number of fields in the system
	static const unsigned int nvar = 4;

	// boundary at X and Y
	static const bool boundary[];

	// type of space float, double, ...
	typedef float stype;

	// type of base grid
	typedef grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> b_grid;

	// type of SparseMatrix for the linear solver
	typedef SparseMatrix<double,int> SparseMatrix_type;

	// type of Vector for the linear solver
	typedef Vector<double> Vector_type;

	// Define the underline grid is staggered
	static const int grid_type = STAGGERED_GRID;
};

struct lid_nn_3d_petsc
{
	// dimensionaly of the equation ( 3D problem ...)
	static const unsigned int dims = 3;
	// number of fields in the system
	static const unsigned int nvar = 4;

	// boundary at X and Y
	static const bool boundary[];

	// type of space float, double, ...
	typedef float stype;

	// type of base grid
	typedef grid_dist_id<3,float,aggregate<float[3],float>,CartDecomposition<3,float>> b_grid;

	// type of SparseMatrix for the linear solver
	typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;

	// type of Vector for the linear solver
	typedef Vector<double,PETSC_BASE> Vector_type;

	// Define the underline grid is staggered
	static const int grid_type = STAGGERED_GRID;
};

//typedef lid_nn_3d_eigen lid_nn_3d;

const bool lid_nn_3d_eigen::boundary[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
const bool lid_nn_3d_petsc::boundary[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

// Constant Field
struct eta
{
	typedef void const_field;

	static float val()	{return 1.0;}
};

//#define SYSEQ_TYPE lid_nn_3d_eigen;
//#include "Equations/stoke_flow_eq_3d.hpp"

template<typename solver_type,typename lid_nn_3d> void lid_driven_cavity_3d()
{
	#include "Equations/stoke_flow_eq_3d.hpp"

	Vcluster & v_cl = create_vcluster();

//	if (v_cl.getProcessingUnits() > 3)
//		return;

	// Domain
	Box<3,float> domain({0.0,0.0,0.0},{3.0,1.0,1.0});

	// Ghost
	Ghost<3,float> g(0.01);

	long int sz[] = {96,64,64};
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
	FDScheme<lid_nn_3d> fd(pd,stencil_max,domain,g_dist.getGridInfo(),g_dist);

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

	std::cout << "Solving " << "\n";

	solver_type solver;
	solver.best_solve();
	auto x = solver.solve(fd.getA(),fd.getB());

	std::cout << "Solved" << "\n";

	// Bring the solution to grid
	fd.template copy<velocity,pressure>(x,{0,0},{sz[0]-1,sz[1]-1,sz[2]-1},g_dist);

	g_dist.write("lid_driven_cavity_3d_p" + std::to_string(v_cl.getProcessingUnits()));

/*	if (v_cl.getProcessUnitID() == 0)
	{
		if (v_cl.getProcessingUnits() == 1)
		{
			// Check that match
			bool test = compare("lid_driven_cavity_3d_p1_grid_0_test.vtk","lid_driven_cavity_3d_p1_grid_0.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
		}
		else if (v_cl.getProcessingUnits() == 2)
		{
			// Check that match
			bool test = compare("lid_driven_cavity_p2_grid_0_test.vtk","lid_driven_cavity_p2_grid_0.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
			test = compare("lid_driven_cavity_p2_grid_1_test.vtk","lid_driven_cavity_p2_grid_1.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
		}
		else if (v_cl.getProcessingUnits() == 3)
		{
			// Check that match
			bool test = compare("lid_driven_cavity_p3_grid_0_test.vtk","lid_driven_cavity_p3_grid_0.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
			test = compare("lid_driven_cavity_p3_grid_1_test.vtk","lid_driven_cavity_p3_grid_1.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
			test = compare("lid_driven_cavity_p3_grid_2_test.vtk","lid_driven_cavity_p3_grid_2.vtk");
			BOOST_REQUIRE_EQUAL(test,true);
		}
	}*/
}

// Lid driven cavity, uncompressible fluid

BOOST_AUTO_TEST_CASE(lid_driven_cavity)
{
//	lid_driven_cavity_3d<umfpack_solver<double>,lid_nn_3d_eigen>();
	lid_driven_cavity_3d<petsc_solver<double>,lid_nn_3d_petsc>();
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_3D_HPP_ */
