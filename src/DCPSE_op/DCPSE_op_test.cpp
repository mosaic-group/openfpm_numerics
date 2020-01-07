/*
 * DCPSE_op_test.cpp
 *
 *  Created on: Jan 7, 2020
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK
#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "DCPSE_op.hpp"

//template<typename T>
//struct Debug;

BOOST_AUTO_TEST_SUITE( dcpse_op_suite_tests )

BOOST_AUTO_TEST_CASE( dcpse_op_tests )
{
	Lap L;
	Field F;

	auto expr = L(L(F)) + L(F) + F;

	std::cout << demangle(typeid(decltype(expr)).name()) << "\n";

	//Debug<decltype(expr)> a;

	typedef decltype(expr)::blabla blabla;

	auto err = L + L;
}

BOOST_AUTO_TEST_SUITE_END()


