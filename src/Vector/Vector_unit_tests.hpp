/*
 * Vector_unit_tests.hpp
 *
 *  Created on: Apr 6, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_UNIT_TESTS_HPP_

#include "Vector/Vector.hpp"

BOOST_AUTO_TEST_SUITE( vector_test_suite )

BOOST_AUTO_TEST_CASE(vector_eigen_parallel)
{
	Vcluster & vcl = create_vcluster();

	if (vcl.getProcessingUnits() != 3)
		return;

	// 3 Processors 9x9 Matrix to invert

	Vector<double> v(9);

	if (vcl.getProcessUnitID() == 0)
	{
		// v row1
		v.insert(0,0);
		v.insert(1,1);
		v.insert(2,2);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		v.insert(3,3);
		v.insert(4,4);
		v.insert(5,5);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		v.insert(6,6);
		v.insert(7,7);
		v.insert(8,8);
	}

	Vector<double> v3;
	v3 = v;

	// force to sync
	v.getVec();

	Vector<double> v2;
	v2 = v;

	// Master has the full vector

	if (vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(v(0),0);
		BOOST_REQUIRE_EQUAL(v(1),1);
		BOOST_REQUIRE_EQUAL(v(2),2);

		BOOST_REQUIRE_EQUAL(v(3),3);
		BOOST_REQUIRE_EQUAL(v(4),4);
		BOOST_REQUIRE_EQUAL(v(5),5);

		BOOST_REQUIRE_EQUAL(v(6),6);
		BOOST_REQUIRE_EQUAL(v(7),7);
		BOOST_REQUIRE_EQUAL(v(8),8);

		// Change the vector on Master

		v(0) = 8;
		v(1) = 7;
		v(2) = 6;
		v(3) = 5;
		v(4) = 4;
		v(5) = 3;
		v(6) = 2;
		v(7) = 1;
		v(8) = 0;
	}

	v.scatter();
	v2.scatter();
	v3.scatter();

	if (vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(v(0),8);
		BOOST_REQUIRE_EQUAL(v(1),7);
		BOOST_REQUIRE_EQUAL(v(2),6);

		BOOST_REQUIRE_EQUAL(v2(0),0);
		BOOST_REQUIRE_EQUAL(v2(1),1);
		BOOST_REQUIRE_EQUAL(v2(2),2);

		BOOST_REQUIRE_EQUAL(v3(0),0);
		BOOST_REQUIRE_EQUAL(v3(1),1);
		BOOST_REQUIRE_EQUAL(v3(2),2);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_EQUAL(v(3),5);
		BOOST_REQUIRE_EQUAL(v(4),4);
		BOOST_REQUIRE_EQUAL(v(5),3);

		BOOST_REQUIRE_EQUAL(v2(3),3);
		BOOST_REQUIRE_EQUAL(v2(4),4);
		BOOST_REQUIRE_EQUAL(v2(5),5);

		BOOST_REQUIRE_EQUAL(v3(0),3);
		BOOST_REQUIRE_EQUAL(v3(1),4);
		BOOST_REQUIRE_EQUAL(v3(2),5);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(v(6),2);
		BOOST_REQUIRE_EQUAL(v(7),1);
		BOOST_REQUIRE_EQUAL(v(8),0);

		BOOST_REQUIRE_EQUAL(v2(6),6);
		BOOST_REQUIRE_EQUAL(v2(7),7);
		BOOST_REQUIRE_EQUAL(v2(8),8);

		BOOST_REQUIRE_EQUAL(v3(0),6);
		BOOST_REQUIRE_EQUAL(v3(1),7);
		BOOST_REQUIRE_EQUAL(v3(2),8);
	}
}


BOOST_AUTO_TEST_CASE(vector_petsc_parallel)
{
	Vcluster & vcl = create_vcluster();

	if (vcl.getProcessingUnits() != 3)
		return;

	// 3 Processors 9x9 Matrix to invert

	Vector<double,OFPM_PETSC_VEC> v(9);

	if (vcl.getProcessUnitID() == 0)
	{
		// v row1
		v.insert(0,0);
		v.insert(1,1);
		v.insert(2,2);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		v.insert(3,3);
		v.insert(4,4);
		v.insert(5,5);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		v.insert(6,6);
		v.insert(7,7);
		v.insert(8,8);
	}

	Vector<double,OFPM_PETSC_VEC> v3;
	v3 = v;

	if (vcl.getProcessUnitID() == 0)
	{
		// v row1
		v(0) = 8;
		v(1) = 7;
		v(2) = 6;
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		v(3) = 5;
		v(4) = 4;
		v(5) = 3;
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		v(6) = 2;
		v(7) = 1;
		v(8) = 0;
	}

	// Master has the full vector

	if (vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(v(0),8);
		BOOST_REQUIRE_EQUAL(v(1),7);
		BOOST_REQUIRE_EQUAL(v(2),6);

		BOOST_REQUIRE_EQUAL(v3(0),0);
		BOOST_REQUIRE_EQUAL(v3(1),1);
		BOOST_REQUIRE_EQUAL(v3(2),2);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_EQUAL(v(3),5);
		BOOST_REQUIRE_EQUAL(v(4),4);
		BOOST_REQUIRE_EQUAL(v(5),3);

		BOOST_REQUIRE_EQUAL(v3(0),3);
		BOOST_REQUIRE_EQUAL(v3(1),4);
		BOOST_REQUIRE_EQUAL(v3(2),5);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(v(6),2);
		BOOST_REQUIRE_EQUAL(v(7),1);
		BOOST_REQUIRE_EQUAL(v(8),0);

		BOOST_REQUIRE_EQUAL(v3(0),6);
		BOOST_REQUIRE_EQUAL(v3(1),7);
		BOOST_REQUIRE_EQUAL(v3(2),9);
	}

	auto & v2 = v.getVec();
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_UNIT_TESTS_HPP_ */
