/*
 * operators_unit_tests.cpp
 *
 *  Created on: Dec 09, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_UNIT_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_UNIT_TEST_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "FiniteDifference/Derivative.hpp"
#include "FiniteDifference/Laplacian.hpp"
#include "FiniteDifference/Average.hpp"
#include "FiniteDifference/sum.hpp"
#include "FiniteDifference/mul.hpp"
#include "eq.hpp"
#include "FiniteDifference/operators.hpp"

struct  op_sys_nn
{
	//! dimensionaly of the equation (2D problem 3D problem ...)
	static const unsigned int dims = 2;
	//! number of degree of freedoms
	static const unsigned int nvar = 1;

	static const unsigned int ord = EQS_FIELD;

	//! boundary at X and Y
	static const bool boundary[];

	//! type of space float, double, ...
	typedef float stype;

	//! Base grid
	typedef grid_dist_id<dims,stype,aggregate<float>,CartDecomposition<2,stype> > b_grid;

	//! specify that we are on testing
	typedef void testing;
};

const bool op_sys_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

constexpr unsigned int x = 0;
constexpr unsigned int y = 1;
constexpr unsigned int z = 2;
constexpr unsigned int V = 0;

template<typename T>
struct debug;

BOOST_AUTO_TEST_SUITE( operators_test )

BOOST_AUTO_TEST_CASE( operator_plus )
{
  Field<x,op_sys_nn> f1;
  Field<y,op_sys_nn> f2;
  Field<z,op_sys_nn> f3;

  auto sum2 = f1+f2;
  auto sum3 = f1+f2+f3;
  
  auto mul2 = f1*f2;
  auto mul3 = f2*f1*f3;

  auto test1 = f2*(f3+f1);

  Der<x,op_sys_nn,CENTRAL> dx;
  auto der1 = dx(f1);
  auto der2 = dx(sum2);

  coeff<double,op_sys_nn> c{3.0};
  auto coeff_test = c*f1;

  // debug<decltype(f1*f1)> dbg;
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_UNIT_TEST_HPP_ */

