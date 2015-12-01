/*
 * Matrix.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_HPP_

#include <Eigen/Sparse>

/*! \brief It store the non zero elements of the matrix
 *
 *
 */
template<typename T> struct cval
{
	size_t j;
	T value;
};

/*! \brief It store the non zero elements of the matrix
 *
 *
 */
template<typename T, unsigned int impl> struct triplet
{
	long int i;
	long int j;
	T val;

	long int & row()
	{
		return i;
	}

	long int & col()
	{
		return j;
	}

	T & value()
	{
		return val;
	}
};

/* ! \brief Sparse Matrix implementation
 *
 * \tparam T Type of the sparse Matrix
 * \tparam id_t type of id
 * \tparam impl implementation
 *
 */
template<typename T,typename id_t ,typename Mi = Eigen::SparseMatrix<T,0,id_t>>
class SparseMatrix
{
};

#include "SparseMatrix_Eigen.hpp"

#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_HPP_ */
