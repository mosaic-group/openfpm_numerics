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
#include "Solvers/petsc_solver.hpp"

#ifdef HAVE_PETSC
#include <petscmat.h>
#endif

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
		BOOST_REQUIRE_CLOSE(x(6), -10.5, 0.001);
		BOOST_REQUIRE_CLOSE(x(7), -8, 0.001);
		BOOST_REQUIRE_CLOSE(x(8), -4.5, 0.001);
	}
}

BOOST_AUTO_TEST_CASE(sparse_matrix_eigen_petsc)
{
	Vcluster & vcl = create_vcluster();

	if (vcl.getProcessingUnits() != 3)
		return;

	// 3 Processors 9x9 Matrix to invert

	SparseMatrix<double,int,PETSC_BASE> sm(9,9,3);
	Vector<double,PETSC_BASE> v(9,3);

	typedef SparseMatrix<double,int,PETSC_BASE>::triplet_type triplet;

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

	// Get the petsc Matrix
	Mat & matp = sm.getMat();

	if (vcl.getProcessUnitID() == 0)
	{
		PetscInt r0[1] = {0};
		PetscInt r1[1] = {1};
		PetscInt r2[1] = {2};
		PetscInt r0cols[3] = {0,1};
		PetscInt r1cols[3] = {0,1,2};
		PetscInt r2cols[3] = {1,2,3};
		PetscScalar y[3];

		MatGetValues(matp,1,r0,2,r0cols,y);

		BOOST_REQUIRE_EQUAL(y[0],-2);
		BOOST_REQUIRE_EQUAL(y[1],1);

		MatGetValues(matp,1,r1,3,r1cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);

		MatGetValues(matp,1,r2,3,r2cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);

	}
	else if (vcl.getProcessUnitID() == 1)
	{
		PetscInt r0[1] = {3};
		PetscInt r1[1] = {4};
		PetscInt r2[1] = {5};
		PetscInt r0cols[3] = {2,3,4};
		PetscInt r1cols[3] = {3,4,5};
		PetscInt r2cols[3] = {4,5,6};
		PetscScalar y[3];

		MatGetValues(matp,1,r0,3,r0cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);

		MatGetValues(matp,1,r1,3,r1cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);

		MatGetValues(matp,1,r2,3,r2cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		PetscInt r0[1] = {6};
		PetscInt r1[1] = {7};
		PetscInt r2[1] = {8};
		PetscInt r0cols[3] = {5,6,7};
		PetscInt r1cols[3] = {6,7,8};
		PetscInt r2cols[3] = {7,8};
		PetscScalar y[3];

		MatGetValues(matp,1,r0,3,r0cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);

		MatGetValues(matp,1,r1,3,r1cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);
		BOOST_REQUIRE_EQUAL(y[2],1);

		MatGetValues(matp,1,r2,2,r2cols,y);

		BOOST_REQUIRE_EQUAL(y[0],1);
		BOOST_REQUIRE_EQUAL(y[1],-2);

	}

	petsc_solver<double> solver;
	solver.solve(sm,v);
}

BOOST_AUTO_TEST_CASE(sparse_matrix_eigen_petsc_solve)
{
	Vcluster & vcl = create_vcluster();

	if (vcl.getProcessingUnits() != 3)
		return;

	const int loc = 200;

	// 3 Processors 9x9 Matrix to invert

	SparseMatrix<double,int,PETSC_BASE> sm(3*loc,3*loc,loc);
	Vector<double,PETSC_BASE> v(3*loc,loc);

	typedef SparseMatrix<double,int,PETSC_BASE>::triplet_type triplet;

	auto & triplets = sm.getMatrixTriplets();

	if (vcl.getProcessUnitID() == 0)
	{
		// row 1
		triplets.add(triplet(0,0,-2));
		triplets.add(triplet(0,1,1));
		v.insert(0,1);

		// row 2
		for (size_t i = 1 ; i < loc ; i++)
		{
			triplets.add(triplet(i,i-1,1));
			triplets.add(triplet(i,i,-2));
			triplets.add(triplet(i,i+1,1));

			v.insert(i,1);
		}
	}
	else if (vcl.getProcessUnitID() == 1)
	{
		// row 2
		for (size_t i = loc ; i < 2*loc ; i++)
		{
			triplets.add(triplet(i,i-1,1));
			triplets.add(triplet(i,i,-2));
			triplets.add(triplet(i,i+1,1));

			v.insert(i,1);
		}
	}
	else if (vcl.getProcessUnitID() == 2)
	{
		// row 2
		for (size_t i = 2*loc ; i < 3*loc-1 ; i++)
		{
			triplets.add(triplet(i,i-1,1));
			triplets.add(triplet(i,i,-2));
			triplets.add(triplet(i,i+1,1));

			v.insert(i,1);
		}

		// row 9
		triplets.add(triplet(3*loc-1,3*loc-2,1));
		triplets.add(triplet(3*loc-1,3*loc-1,-2));
		v.insert(3*loc-1,1);
	}

	// Get the petsc Matrix
	Mat & matp = sm.getMat();


	petsc_solver<double> solver;
	solver.solve(sm,v);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_UNIT_TESTS_HPP_ */
