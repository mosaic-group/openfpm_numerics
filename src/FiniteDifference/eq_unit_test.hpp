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
	typedef grid_dist_id<2,float,scalar<float>,CartDecomposition<2,float>> b_grid;
};

const bool lid_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

// Constant Field
struct eta
{
	float val()	{return 1.0;}
};

// Model the equations

constexpr unsigned int v[] = {0,1};
constexpr unsigned int P = 2;

typedef Field<v[x],lid_nn> v_x;
typedef Field<v[y],lid_nn> v_y;
typedef Field<P,lid_nn> Prs;

// Eq1 V_x

typedef mul<eta,Lap<v_x,lid_nn>,lid_nn> eta_lap_vx;
typedef D<x,Prs,lid_nn> p_x;
typedef sum<eta_lap_vx,p_x,lid_nn> vx_eq;

// Eq2 V_y

typedef mul<eta,Lap<v_y,lid_nn>,lid_nn> eta_lap_vy;
typedef D<x,Prs,lid_nn> p_x;
typedef sum<eta_lap_vy,p_x,lid_nn> vy_eq;

// Eq3 Incompressibility

typedef D<x,v_x,lid_nn> dx_vx;
typedef D<y,v_y,lid_nn> dy_vy;
typedef sum<dx_vx,dy_vy> incompressibility;

BOOST_AUTO_TEST_SUITE( eq_test_suite )

// Lid driven cavity, uncompressible fluid

BOOST_AUTO_TEST_CASE( lid_driven_cavity )
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Ghost
	Ghost<2,float> g(0.01);

	size_t sz[] = {32,32};

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Initialize openfpm
	grid_dist_id<2,float,scalar<float>,CartDecomposition<2,float>> g_dist(sz,domain,g);

	// Distributed grid
	FDScheme<lid_nn> fd;

	// start and end of the bulk
	grid_key_dx<2> bulk_start(1,1);
	grid_key_dx<2> bulk_end(sz[0]-1,sz[1]-1);

	// Impose the vx equation
	vx_eq vx;
	fd.impose(vx, g_dist.getGridInfo() , g_dist.getSubDomainIterator(bulk_start,bulk_end));

	// Impose the vy equation
	vy_eq vy;
	fd.impose(vy, g_dist.getGridInfo(), g_dist.getSubDomainIterator(bulk_start,bulk_end));
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_ */
