/*
 * eq_unit_tests.hpp
 *
 *  Created on: Sep 18, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_UNIT_TESTS_HPP_

#include "FiniteDifference/FDScheme.hpp"

BOOST_AUTO_TEST_SUITE( eq_test )

struct Sys_p
{
	static const unsigned int dims = 3;
	typedef float stype;
};

// define coordinates id

constexpr unsigned int x = 0;
constexpr unsigned int y = 1;

// define constants

constexpr unsigned int eta = 0;
constexpr unsigned int zeta = 1;
constexpr unsigned int gamma = 2;
constexpr unsigned int nu = 3;
constexpr unsigned int lambda = 4;
constexpr unsigned int dmu = 5;
constexpr unsigned int sigma[2][2] = {{6,7},{8,9}};
constexpr unsigned int p[] = {10,11};
constexpr unsigned int del_nu = 12;
constexpr unsigned int g_ext = 13;

// define variables fields

// velocity
constexpr unsigned int V[2] = {0,1};

// pressure
constexpr unsigned int P = 2;

// model Eq1 active Gel

//// First level expression



/*
 *
 * p_x^2
 *
 */
//typedef mul<p[x],p[x],Sys_p> px2;


/*
 *
 * p_y^2
 *
 */
//typedef mul<p[y],p[y],Sys_p> py2;


/*
 * p_x * p_y
 *
 */
//typedef mul<p[x],p[y],Sys_p> px_py;



BOOST_AUTO_TEST_CASE( eq_test_use)
{
	// create an equation object


}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_UNIT_TESTS_HPP_ */
