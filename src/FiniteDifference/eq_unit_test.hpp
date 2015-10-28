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
typedef sum<eta_lap_vx,p_x,lid_nn> vx_eq;

// Eq2 V_y

typedef mul<eta,Lap<v_y,lid_nn>,lid_nn> eta_lap_vy;
typedef D<x,Prs,lid_nn> p_x;
typedef sum<eta_lap_vy,p_x,lid_nn> vy_eq;

// Eq3 Incompressibility

typedef D<x,v_x,lid_nn> dx_vx;
typedef D<y,v_y,lid_nn> dy_vy;
typedef sum<dx_vx,dy_vy,lid_nn> ic_eq;


BOOST_AUTO_TEST_SUITE( eq_test_suite )

// Lid driven cavity, uncompressible fluid

BOOST_AUTO_TEST_CASE( lid_driven_cavity )
{
	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Ghost
	Ghost<2,float> g(0.01);

	size_t sz[] = {32,32};
	Padding<2> pd({1,2},{1,2});

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Initialize openfpm
	grid_dist_id<2,float,scalar<float>,CartDecomposition<2,float>> g_dist(sz,domain,g);

	// Distributed grid
	FDScheme<lid_nn> fd(pd,domain,g_dist.getGridInfo(),g_dist.getDecomposition(),g_dist.getVC());

	// start and end of the bulk
	grid_key_dx<2> bulk_start(0,0);
	grid_key_dx<2> bulk_end(sz[0],sz[1]);

	// Impose the vx and vy equation in the bulk
	fd.imposeA(vx_eq(), g_dist.getSubDomainIterator(bulk_start,bulk_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bulk_start,bulk_end));
	fd.imposeA(vy_eq(), g_dist.getSubDomainIterator(bulk_start,bulk_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bulk_start,bulk_end));
	fd.imposeA(ic_eq(), g_dist.getSubDomainIterator(bulk_start,bulk_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bulk_start,bulk_end));

	// Boundary condition on the left
	grid_key_dx<2> bleft_start(0,-1);
	grid_key_dx<2> bleft_end(0,33);
	fd.imposeA(v_x(), g_dist.getSubDomainIterator(bleft_start,bleft_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bleft_start,bleft_end));

	// Boundary condition on the right
	grid_key_dx<2> bright_start(32,-1);
	grid_key_dx<2> bright_end(32,33);
	fd.imposeA(v_x(), g_dist.getSubDomainIterator(bright_start,bright_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bright_start,bright_end));

	// Padding condition Pressure Left
	grid_key_dx<2> pleft_start(-1,0);
	grid_key_dx<2> pleft_end(-1,31);
	fd.imposeA(Prs(), g_dist.getSubDomainIterator(bleft_start,bleft_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bleft_start,bleft_end));

	// Padding condition Pressure Right
	grid_key_dx<2> pright_start(32,-1);
	grid_key_dx<2> pright_end(33,33);
	fd.imposeA(Prs(), g_dist.getSubDomainIterator(pright_start,pright_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(pright_start,pright_end));

	// Padding condition vy Right
	grid_key_dx<2> pvright_start(33,-1);
	grid_key_dx<2> pvright_end(33,33);
	fd.imposeA(v_y(), g_dist.getSubDomainIterator(pvright_start,pright_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(pvright_start,pright_end));

	// Padding Bottom pressure
	grid_key_dx<2> pbot_start(0,-1);
	grid_key_dx<2> pbot_end(31,-1);
	fd.imposeA(Prs(), g_dist.getSubDomainIterator(pbot_start,pbot_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(pbot_start,pbot_end));

	// Bottom boundary condition vy
	grid_key_dx<2> bbot_start(-1,0);
	grid_key_dx<2> bbot_end(33,0);
	fd.imposeA(v_y(), g_dist.getSubDomainIterator(bbot_start,bbot_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(bbot_start,bbot_end));

	// Padding top pressure
	grid_key_dx<2> ptop_start(0,32);
	grid_key_dx<2> ptop_end(31,33);
	fd.imposeA(Prs(), g_dist.getSubDomainIterator(ptop_start,ptop_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(ptop_start,ptop_end));

	// Top boundary condition v_y
	grid_key_dx<2> btop_start(-1,32);
	grid_key_dx<2> btop_end(33,32);
	fd.imposeA(v_y(), g_dist.getSubDomainIterator(btop_start,btop_end));
	fd.imposeB(1.0, g_dist.getSubDomainIterator(btop_start,btop_end));

	// Padding top v_x
	grid_key_dx<2> pvtop_start(0,33);
	grid_key_dx<2> pvtop_end(31,33);
	fd.imposeA(v_x(), g_dist.getSubDomainIterator(pvtop_start,pvtop_end));
	fd.imposeB(0.0, g_dist.getSubDomainIterator(pvtop_start,pvtop_end));
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_ */
