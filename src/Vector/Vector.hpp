/*
 * Vector.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_

#include <Eigen/Sparse>

/*! \brief It store one row value of a vector
 *
 * Given a row, store a value
 *
 *
 */
template<typename T, unsigned int impl>
class rval
{
};


/* ! \brief Sparse Matrix implementation
 *
 * \tparam T Type of the sparse Matrix
 * \tparam impl implementation
 *
 */
template<typename T,typename impl = Eigen::Matrix<T, Eigen::Dynamic, 1> >
class Vector
{
};

#include "Vector_eigen.hpp"
#include "Vector_petsc.hpp"

#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_ */
