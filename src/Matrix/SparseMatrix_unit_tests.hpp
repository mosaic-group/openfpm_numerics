/*
 * SparseMatrix_unit_tests.hpp
 *
 *  Created on: Apr 4, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_UNIT_TESTS_HPP_

#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "Solvers/umfpack_solver.hpp"

BOOST_AUTO_TEST_SUITE( sparse_matrix_test_suite )



BOOST_AUTO_TEST_CASE(sparse_matrix_eigen_parallel)
{
	Vcluster & vcl = create_vcluster();

	if (vcl.getProcessingUnits() != 3)
		return;

	// 3 Processors 9x9 Matrix to invert

	SparseMatrix<double,int> sm(9,9);
	Vector<double> v(9);

	typedef SparseMatrix<double,int>::triplet_type triplet;

	auto & triplets = sm.getMatrixTriplets();

	if (vcl.getProcessUnitID() == 0)
	{
		// row 1
		triplets.add(triplet(0,0,-2));
		triplets.add(triplet(0,1,1));

		// row 2
		triplets.add(triplet(1,0,1));
		triplets.add(triplet(1,1,-2));
		triplets.add(triplet(1,2,1));

		// row 3
		triplets.add(triplet(2,1,1));
		triplets.add(triplet(2,2,-2));
		triplets.add(triplet(2,3,1));

		// v row1
		v.insert(0,1);
		v.insert(1,1);
		v.insert(2,1);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		// row 4
		triplets.add(triplet(3,2,1));
		triplets.add(triplet(3,3,-2));
		triplets.add(triplet(3,4,1));

		// row 5
		triplets.add(triplet(4,3,1));
		triplets.add(triplet(4,4,-2));
		triplets.add(triplet(4,5,1));

		// row 6
		triplets.add(triplet(5,4,1));
		triplets.add(triplet(5,5,-2));
		triplets.add(triplet(5,6,1));

		v.insert(3,1);
		v.insert(4,1);
		v.insert(5,1);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		// row 7
		triplets.add(triplet(6,5,1));
		triplets.add(triplet(6,6,-2));
		triplets.add(triplet(6,7,1));

		// row 8
		triplets.add(triplet(7,6,1));
		triplets.add(triplet(7,7,-2));
		triplets.add(triplet(7,8,1));

		// row 9
		triplets.add(triplet(8,7,1));
		triplets.add(triplet(8,8,-2));

		v.insert(6,1);
		v.insert(7,1);
		v.insert(8,1);
	}

	// force to sync
	sm.getMat();

	// Master has the full Matrix

	if (vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(sm(0,0),-2);
		BOOST_REQUIRE_EQUAL(sm(0,1),1);

		BOOST_REQUIRE_EQUAL(sm(1,0),1);
		BOOST_REQUIRE_EQUAL(sm(1,1),-2);
		BOOST_REQUIRE_EQUAL(sm(1,2),1);

		BOOST_REQUIRE_EQUAL(sm(2,1),1);
		BOOST_REQUIRE_EQUAL(sm(2,2),-2);
		BOOST_REQUIRE_EQUAL(sm(2,3),1);

		BOOST_REQUIRE_EQUAL(sm(3,2),1);
		BOOST_REQUIRE_EQUAL(sm(3,3),-2);
		BOOST_REQUIRE_EQUAL(sm(3,4),1);

		BOOST_REQUIRE_EQUAL(sm(4,3),1);
		BOOST_REQUIRE_EQUAL(sm(4,4),-2);
		BOOST_REQUIRE_EQUAL(sm(4,5),1);

		BOOST_REQUIRE_EQUAL(sm(5,4),1);
		BOOST_REQUIRE_EQUAL(sm(5,5),-2);
		BOOST_REQUIRE_EQUAL(sm(5,6),1);

		BOOST_REQUIRE_EQUAL(sm(6,5),1);
		BOOST_REQUIRE_EQUAL(sm(6,6),-2);
		BOOST_REQUIRE_EQUAL(sm(6,7),1);

		BOOST_REQUIRE_EQUAL(sm(7,6),1);
		BOOST_REQUIRE_EQUAL(sm(7,7),-2);
		BOOST_REQUIRE_EQUAL(sm(7,8),1);

		BOOST_REQUIRE_EQUAL(sm(8,7),1);
		BOOST_REQUIRE_EQUAL(sm(8,8),-2);
	}

	// try to invert the Matrix with umfpack

	auto x = umfpack_solver<double>::solve(sm,v);

	// we control the solution

	if (vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_CLOSE(x(0), -4.5, 0.001);
		BOOST_REQUIRE_CLOSE(x(1), -8, 0.001);
		BOOST_REQUIRE_CLOSE(x(2), -10.5, 0.001);
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_CLOSE(x(3), -12.0, 0.001);
		BOOST_REQUIRE_CLOSE(x(4), -12.5, 0.001);
		BOOST_REQUIRE_CLOSE(x(5), -12.0, 0.001);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_CLOSE(x(6), -4.5, 0.001);
		BOOST_REQUIRE_CLOSE(x(7), -8, 0.001);
		BOOST_REQUIRE_CLOSE(x(8), -10.5, 0.001);
	}
}



BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_UNIT_TESTS_HPP_ */
