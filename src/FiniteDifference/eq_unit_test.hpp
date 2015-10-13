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

	// type of space float, double, ...
	typedef float stype;
};

const bool lid_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

// Constant Field

struct eta
{
};

// Model the equation

constexpr unsigned int v[] = {0,1};
constexpr unsigned int P = 2;

typedef Field<v[x],lid_nn> v_x;
typedef Field<v[y],lid_nn> v_y;
typedef Field<P,lid_nn> Prs;

// Eq1 V_x

typedef mul<eta,Lap<v_x,lid_nn>,lid_nn> eta_lap_vx;
typedef D<x,Prs,lid_nn> p_x;
typedef sum<eta_lap_vx,p_x> vx_eq;

// Eq2 V_y

typedef mul<eta,Lap<v_y,lid_nn>,lid_nn> eta_lap_vy;
typedef D<x,Prs,lid_nn> p_x;
typedef sum<eta_lap_vy,p_x,lid_nn> vy_eq;

// Eq3 Incompressibility

typedef D<x,v_x,lid_nn> dx_vx;
typedef D<y,v_y,lid_nn> dy_vy;
typedef sum<dx_vx,dy_vy> incompressibility;

// Lid driven cavity, uncompressible fluid

BOOST_AUTO_TEST_CASE( lid_driven_cavity )
{

}


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_UNIT_TEST_HPP_ */
