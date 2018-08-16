/*
 * util_num_unit_tests.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_UTIL_UTIL_NUM_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_UTIL_UTIL_NUM_UNIT_TESTS_HPP_

#include "util_num.hpp"

//! [Constant fields struct definition]

//! define a constant field
struct anyname_field
{
	//! define that is a constant field
	typedef void const_field;

	//! Evaluate the constant field (in this case always return 1.0)
	float val()	{return 1.0;}
};

//! define a non-constant (in space) field
struct anyname_field_with_pos
{
	//! It define that is a constant field
	typedef void const_field;

	//! It define that is not constant in space
	typedef void with_position;

	//! Evaluate the field in one point
	float val_pos(grid_key_dx<2> & pos)	{return 1.0;}
};

//! stub field
struct no_field
{
};

//! [Constant fields struct definition]

//! [Is on test struct definition]

//! define that we are on testing mode
struct on_test
{
	//! specify testing mode
	typedef void testing;
};

//! Not on testing mode
struct not_on_test
{
};

//! [Is on test struct definition]

BOOST_AUTO_TEST_SUITE( util_num_suite )

BOOST_AUTO_TEST_CASE( util_num )
{
	//! [Usage of is_const_field]

	bool ret = is_const_field<anyname_field>::value;
	BOOST_REQUIRE_EQUAL(ret,true);

	ret = is_const_field<no_field>::value;
	BOOST_REQUIRE_EQUAL(ret,false);

	//! [Usage of is_const_field]

	//! [Usage of has_val_pos]

	ret = is_const_field<anyname_field_with_pos>::value;
	BOOST_REQUIRE_EQUAL(ret,true);

	ret = is_const_field<no_field>::value;
	BOOST_REQUIRE_EQUAL(ret,false);

	//! [Usage of has_val_pos]

	//! [Usage of is_testing]

	ret = is_testing<on_test>::value;
	BOOST_REQUIRE_EQUAL(ret,true);

	ret = is_testing<not_on_test>::value;
	BOOST_REQUIRE_EQUAL(ret,false);

	//! [Usage of is_testing]

	//! [Usage of stub_or_real]

	ret = std::is_same<stub_or_real<on_test,2,float,CartDecomposition<2,float>>::type, grid_dist_testing<2>>::value ;
	BOOST_REQUIRE_EQUAL(ret,true);

	ret = std::is_same<stub_or_real<not_on_test,2,float,CartDecomposition<2,float>>::type, grid_dist_id<2,float,aggregate<size_t>,CartDecomposition<2,float>> >::value;
	BOOST_REQUIRE_EQUAL(ret,true);

	//! [Usage of stub_or_real]
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_UTIL_UTIL_NUM_UNIT_TESTS_HPP_ */
