/*
 * SparseMatrix_petsc.hpp
 *
 *  Created on: Apr 26, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_PETSC_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_PETSC_HPP_


#include "Vector/map_vector.hpp"
#include <boost/mpl/int.hpp>
#include <petscmat.h>

#define PETSC_BASE 2

/*! \brief It store one non-zero element in the sparse matrix
 *
 * Given a row, and a column, store a value
 *
 *
 */
template<typename T>
class triplet<T,PETSC_BASE>
{
	PetscInt row_;
	PetscInt col_;
	PetscScalar val_;

public:

	PetscInt & row()
	{
		return row_;
	}

	PetscInt & col()
	{
		return col_;
	}

	PetscScalar & value()
	{
		return val_;
	}

	/*! \brief Constructor from row, colum and value
	 *
	 * \param i row
	 * \param j colum
	 * \param val value
	 *
	 */
	triplet(long int i, long int j, T val)
	{
		row_ = i;
		col_ = j;
		val_ = val;
	}

	// Default constructor
	triplet()	{};
};

/* ! \brief Sparse Matrix implementation, that map over Eigen
 *
 * \tparam T Type of the sparse Matrix store on each row,colums
 * \tparam id_t type of id
 * \tparam impl implementation
 *
 */
template<typename T, typename id_t>
class SparseMatrix<T,id_t, PETSC_BASE>
{
public:

	//! Triplet implementation id
	typedef boost::mpl::int_<PETSC_BASE> triplet_impl;

	//! Triplet type
	typedef triplet<T,PETSC_BASE> triplet_type;

private:

	// Size of the matrix
	size_t g_row;
	size_t g_col;

	// Local size of the matrix
	size_t l_row;
	size_t l_col;

	// starting row for this processor
	size_t start_row;

	// indicate if the matrix has been created
	bool m_created = false;

	// PETSC Matrix
	Mat mat;
	openfpm::vector<triplet_type> trpl;
	openfpm::vector<triplet_type> trpl_recv;

	// temporary list of values
	mutable openfpm::vector<PetscScalar> vals;
	// temporary list of colums
	mutable openfpm::vector<PetscInt> cols;
	// PETSC d_nnz and o_nnz
	mutable openfpm::vector<PetscInt> d_nnz;
	mutable openfpm::vector<PetscInt> o_nnz;

	/*! \brief Fill the petsc Matrix
	 *
	 *
	 */
	void fill_petsc()
	{
		d_nnz.resize(l_row);
		o_nnz.resize(l_row);

		d_nnz.fill(0);
		o_nnz.fill(0);

		// Here we explore every row to count how many non zero we have in the diagonal matrix part,
		// and the non diagonal matrix part, needed by MatMPIAIJSetPreallocation

		size_t i = 0;

		// Set the Matrix from triplet
		while (i < trpl.size())
		{
			PetscInt row = trpl.get(i).row();

			while(row == trpl.get(i).row() && i < trpl.size())
			{
				if ((size_t)trpl.get(i).col() >= start_row && (size_t)trpl.get(i).col() < start_row + l_row)
					d_nnz.get(row - start_row)++;
				else
					o_nnz.get(row - start_row)++;
				i++;
			}
		}

		PETSC_SAFE_CALL(MatMPIAIJSetPreallocation(mat,0,static_cast<const PetscInt*>(d_nnz.getPointer()),0,
														static_cast<const PetscInt*>(o_nnz.getPointer())));

		// Counter i is zero
		i = 0;

		// Set the Matrix from triplet
		while (i < trpl.size())
		{
			vals.clear();
			cols.clear();

			PetscInt row = trpl.get(i).row();

			while(row == trpl.get(i).row() && i < trpl.size())
			{
				vals.add(trpl.get(i).value());
				cols.add(trpl.get(i).col());
				i++;
			}
			PETSC_SAFE_CALL(MatSetValues(mat,1,&row,cols.size(),static_cast<const PetscInt*>(cols.getPointer()),
					                                            static_cast<const PetscScalar *>(vals.getPointer()),
					                                            INSERT_VALUES));
		}

		PETSC_SAFE_CALL(MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY));
		PETSC_SAFE_CALL(MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY));
	}


public:

	/*! \brief Create an empty Matrix
	 *
	 * \param N1 number of row
	 * \param N2 number of colums
	 * \param N1_loc number of local row
	 *
	 */
	SparseMatrix(size_t N1, size_t N2, size_t n_row_local)
	:l_row(n_row_local),l_col(n_row_local)
	{
		PETSC_SAFE_CALL(MatCreate(PETSC_COMM_WORLD,&mat));
		PETSC_SAFE_CALL(MatSetType(mat,MATMPIAIJ));
		PETSC_SAFE_CALL(MatSetSizes(mat,n_row_local,n_row_local,N1,N2));

		Vcluster & v_cl = create_vcluster();

		openfpm::vector<size_t> vn_row_local;
		v_cl.allGather(l_row,vn_row_local);
		v_cl.execute();

		// Calculate the starting row for this processor

		start_row = 0;
		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
			start_row += vn_row_local.get(i);
	}

	/*! \brief Create an empty Matrix
	 *
	 */
	SparseMatrix()
	{
		PETSC_SAFE_CALL(MatCreate(PETSC_COMM_WORLD,&mat));
		PETSC_SAFE_CALL(MatSetType(mat,MATMPIAIJ));
	}

	~SparseMatrix()
	{
		// Destroy the matrix
		PETSC_SAFE_CALL(MatDestroy(&mat));
	}

	/*! \brief Get the Matrix triplets bugger
	 *
	 * It return a buffer that can be filled with triplets
	 *
	 * \return Petsc Matrix
	 *
	 */
	openfpm::vector<triplet_type> & getMatrixTriplets()
	{
		return this->trpl;
	}

	/*! \brief Get the Patsc Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	const Mat & getMat() const
	{
		fill_petsc();

		return mat;
	}

	/*! \brief Get the Petsc Matrix object
	 *
	 * \return the Petsc Matrix
	 *
	 */
	Mat & getMat()
	{
		fill_petsc();

		return mat;
	}

	/*! \brief Resize the Sparse Matrix
	 *
	 * \param row number for row
	 * \param col number of colums
	 * \param local number of row
	 * \param local number of colums
	 *
	 */
	void resize(size_t row, size_t col, size_t l_row, size_t l_col)
	{
		this->g_row = row;
		this->g_col = col;

		this->l_row = l_row;
		this->l_col = l_col;

		PETSC_SAFE_CALL(MatSetSizes(mat,l_row,l_col,g_row,g_col));

		Vcluster & v_cl = create_vcluster();

		openfpm::vector<size_t> vn_row_local;
		v_cl.allGather(l_row,vn_row_local);
		v_cl.execute();

		// Calculate the starting row for this processor

		start_row = 0;
		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
			start_row += vn_row_local.get(i);
	}

	/*! \brief Get the row i and the colum j of the Matrix
	 *
	 * \warning it is slow, consider to get blocks of the matrix
	 *
	 * \return the value of the matrix at row i colum j
	 *
	 */
	T operator()(id_t i, id_t j)
	{
		T v;

		MatGetValues(mat,1,&i,1,&j,&v);

		return v;
	}

	/*! \brief Get the value from triplet
	 *
	 * \warning It is extremly slow because it do a full search across the triplets elements
	 *
	 * \param r row
	 * \param c colum
	 *
	 */
	T getValue(size_t r, size_t c)
	{
		for (size_t i = 0 ; i < trpl.size() ; i++)
		{
			if (r == (size_t)trpl.get(i).row() && c == (size_t)trpl.get(i).col())
				return trpl.get(i).value();
		}

		return 0;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_PETSC_HPP_ */
