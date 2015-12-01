/*
 * common_test.hpp
 *
 *  Created on: Oct 11, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_UTIL_COMMON_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_UTIL_COMMON_TEST_HPP_

#include "FiniteDifference/util/common.hpp"

//! [Define structures]

struct test_grid_type_staggered
{
	static const int grid_type = STAGGERED_GRID;
};

struct test_grid_type_normal
{
	static const int grid_type = NORMAL_GRID;
};

struct test_grid_type_no_def
{
};

//! [Define structures]

BOOST_AUTO_TEST_SUITE( util_pdata_test )

BOOST_AUTO_TEST_CASE( check_pdata_templates_util_function )
{
	{
	//! [Check has grid_type]

	BOOST_REQUIRE_EQUAL(has_grid_type<test_grid_type_staggered>::type::value,true);
	BOOST_REQUIRE_EQUAL(has_grid_type<test_grid_type_normal>::type::value,true);
	BOOST_REQUIRE_EQUAL(has_grid_type<test_grid_type_no_def>::type::value,false);

	//! [Check has grid_type]
	}

	{
	//! [Check grid_type staggered]

	BOOST_REQUIRE_EQUAL(is_grid_staggered<test_grid_type_staggered>::value(),true);
	BOOST_REQUIRE_EQUAL(is_grid_staggered<test_grid_type_normal>::value(),false);
	BOOST_REQUIRE_EQUAL(is_grid_staggered<test_grid_type_no_def>::value(),false);

	//! [Check grid_type staggered]
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_UTIL_COMMON_TEST_HPP_ */
