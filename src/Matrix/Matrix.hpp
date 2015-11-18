/*
 * Matrix.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_MATRIX_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_MATRIX_HPP_

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
template<typename T> struct triplet
{
	long int i;
	long int j;
	T value;
};

template<typename T>
class SparseMatrix
{

};

#endif /* OPENFPM_NUMERICS_SRC_MATRIX_MATRIX_HPP_ */
