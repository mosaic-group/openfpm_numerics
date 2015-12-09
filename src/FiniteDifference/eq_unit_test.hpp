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

// Stokes flow

struct lid_nn
{
	// dimensionaly of the equation (2D problem 3D problem ...)
	static const unsigned int dims = 2;
	// number of fields in the system
	static const unsigned int nvar = 3;
	static const unsigned int ord = EQS_FIELD;

	// boundary at X and Y
	static const bool boundary[];

	//
	static constexpr unsigned int num_cfields = 0;

	// type of space float, double, ...
	typedef float stype;

	// type of base grid
	typedef grid_dist_id<2,float,aggregate<float[2],float>,CartDecomposition<2,float>> b_grid;

	// type of SparseMatrix
	typedef SparseMatrix<double,int> SparseMatrix_type;

	// type of Vector
	typedef Vector<double> Vector_type;

	// Define the the underline grid is staggered
	static const int grid_type = STAGGERED_GRID;
};

const bool lid_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

// Constant Field
struct eta
{
	typedef void const_field;

	static float val()	{return 1.0;}
};

// Model the equations

constexpr unsigned int v[] = {0,1};
constexpr unsigned int P = 2;
constexpr unsigned int ic = 2;

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


// Directional Avg
typedef Avg<x,v_y,lid_nn> avg_vy;
typedef Avg<y,v_x,lid_nn> avg_vx;

typedef Avg<x,v_y,lid_nn,FORWARD> avg_vy_f;
typedef Avg<y,v_x,lid_nn,FORWARD> avg_vx_f;

#define EQ_1 0
#define EQ_2 1
#define EQ_3 2

BOOST_AUTO_TEST_SUITE( eq_test_suite )

// Lid driven cavity, uncompressible fluid

BOOST_AUTO_TEST_CASE(lid_driven_cavity)
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Ghost
	Ghost<2,float> g(0.01);

	long int sz[] = {8,8};
	size_t szu[2];
	szu[0] = (size_t)sz[0];
	szu[1] = (size_t)sz[1];

	Padding<2> pd({1,1},{0,0});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Initialize openfpm
	grid_dist_id<2,float,aggregate<float[2],float>,CartDecomposition<2,float>> g_dist(szu,domain,g);

	// Distributed grid
	FDScheme<lid_nn> fd(pd,domain,g_dist.getGridInfo(),g_dist.getDecomposition(),g_dist.getVC());

	// start and end of the bulk

	fd.impose(ic_eq(),0.0, EQ_3, {0,0},{sz[0]-2,sz[1]-2},true);
	fd.impose(Prs(),  0.0, EQ_3, {0,0},{0,0});
	fd.impose(vx_eq(),0.0, EQ_1, {1,0},{sz[0]-2,sz[1]-2});
	fd.impose(vy_eq(),0.0, EQ_2, {0,1},{sz[0]-2,sz[1]-2});

	// v_x
	// R L
	fd.impose(v_x(),0.0, EQ_1, {0,0},{0,sz[1]-2});
	fd.impose(v_x(),0.0, EQ_1, {sz[0]-1,0},{sz[0]-1,sz[1]-2});

	// T B
	fd.impose(avg_vx_f(),0.0, EQ_1, {0,-1},{sz[0]-1,-1});
	fd.impose(avg_vx(),0.0, EQ_1,   {0,sz[1]-1},{sz[0]-1,sz[1]-1});

	// v_y
	// R L
	fd.impose(avg_vy_f(),0.0, EQ_2 , {-1,0},{-1,sz[1]-1});
	fd.impose(avg_vy(),1.0, EQ_2,    {sz[0]-1,0},{sz[0]-1,sz[1]-1});

	// T B
	fd.impose(v_y(), 0.0, EQ_2, {0,0},{sz[0]-2,0});
	fd.impose(v_y(), 0.0, EQ_2, {0,sz[1]-1},{sz[0]-2,sz[1]-1});

	// Padding pressure
	fd.impose(Prs(), 0.0, EQ_3, {-1,-1},{sz[0]-1,-1});
	fd.impose(Prs(), 0.0, EQ_3, {-1,sz[1]-1},{sz[0]-1,sz[1]-1});
	fd.impose(Prs(), 0.0, EQ_3, {-1,0},{-1,sz[1]-2});
	fd.impose(Prs(), 0.0, EQ_3, {sz[0]-1,0},{sz[0]-1,sz[1]-2});

	// Impose v_x Padding Impose v_y padding
	fd.impose(v_x(), 0.0, EQ_1, {-1,-1},{-1,sz[1]-1});
	fd.impose(v_y(), 0.0, EQ_2, {-1,-1},{sz[0]-1,-1});

	auto x = umfpack_solver<double>::solve(fd.getA(),fd.getB());

	// Bring the solution to grid
	x.copy<FDScheme<lid_nn>,decltype(g_dist),0,1>(fd,{0,0},{sz[0]-1,sz[1]-1},g_dist);

	g_dist.write("lid_driven_cavity");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_ */
