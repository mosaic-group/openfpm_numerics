/*
 * Vector_unit_tests.hpp
 *
 *  Created on: Apr 6, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_UNIT_TESTS_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <VCluster/VCluster.hpp>

#include <iostream>
#include "Vector/Vector.hpp"

BOOST_AUTO_TEST_SUITE( vector_test_suite )

#ifdef HAVE_EIGEN

BOOST_AUTO_TEST_CASE(vector_eigen_parallel)
{
	Vcluster<> & vcl = create_vcluster();

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

		BOOST_REQUIRE_EQUAL(v3(3),3);
		BOOST_REQUIRE_EQUAL(v3(4),4);
		BOOST_REQUIRE_EQUAL(v3(5),5);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(v(6),2);
		BOOST_REQUIRE_EQUAL(v(7),1);
		BOOST_REQUIRE_EQUAL(v(8),0);

		BOOST_REQUIRE_EQUAL(v2(6),6);
		BOOST_REQUIRE_EQUAL(v2(7),7);
		BOOST_REQUIRE_EQUAL(v2(8),8);

		BOOST_REQUIRE_EQUAL(v3(6),6);
		BOOST_REQUIRE_EQUAL(v3(7),7);
		BOOST_REQUIRE_EQUAL(v3(8),8);
	}
}

#endif

#ifdef HAVE_PETSC

BOOST_AUTO_TEST_CASE(vector_petsc_parallel)
{
	Vcluster<> & vcl = create_vcluster();

	if (vcl.getProcessingUnits() != 3)
		return;

	// 3 Processors 9x9 Matrix to invert

	Vector<double,PETSC_BASE> v(9,3);

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

	Vector<double,PETSC_BASE> v3;
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

		BOOST_REQUIRE_EQUAL(v3(3),3);
		BOOST_REQUIRE_EQUAL(v3(4),4);
		BOOST_REQUIRE_EQUAL(v3(5),5);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(v(6),2);
		BOOST_REQUIRE_EQUAL(v(7),1);
		BOOST_REQUIRE_EQUAL(v(8),0);

		BOOST_REQUIRE_EQUAL(v3(6),6);
		BOOST_REQUIRE_EQUAL(v3(7),7);
		BOOST_REQUIRE_EQUAL(v3(8),8);
	}

	// Here we get the petsc vector
	auto & vp = v.getVec();

	// We check the correctness of the PETSC vector

	if (vcl.getProcessUnitID() == 0)
	{
		PetscInt ix[] = {0,1,2};
		PetscScalar y[3];

		VecGetValues(vp,3,ix,y);

		BOOST_REQUIRE_EQUAL(y[0],8);
		BOOST_REQUIRE_EQUAL(y[1],7);
		BOOST_REQUIRE_EQUAL(y[2],6);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		PetscInt ix[] = {3,4,5};
		PetscScalar y[3];

		VecGetValues(vp,3,ix,y);

		BOOST_REQUIRE_EQUAL(y[0],5);
		BOOST_REQUIRE_EQUAL(y[1],4);
		BOOST_REQUIRE_EQUAL(y[2],3);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		PetscInt ix[] = {6,7,8};
		PetscScalar y[3];

		VecGetValues(vp,3,ix,y);

		BOOST_REQUIRE_EQUAL(y[0],2);
		BOOST_REQUIRE_EQUAL(y[1],1);
		BOOST_REQUIRE_EQUAL(y[2],0);
	}

}

#endif

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_UNIT_TESTS_HPP_ */
