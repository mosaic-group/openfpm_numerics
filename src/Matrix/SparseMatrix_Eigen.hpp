/*
 * SparseMatrix_Eigen.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_

#include "Vector/map_vector.hpp"
#include <boost/mpl/int.hpp>

#define EIGEN_TRIPLET 1

/*! \brief It store one non-zero element in the sparse matrix
 *
 * Given a row, and a column, store a value
 *
 *
 */
template<typename T>
struct triplet<T,EIGEN_TRIPLET>
{
	Eigen::Triplet<T,long int> triplet;

	long int & row()
	{
		size_t ptr = (size_t)&triplet.row();

		return (long int &)*(long int *)ptr;
	}

	long int & col()
	{
		size_t ptr = (size_t)&triplet.col();

		return (long int &)*(long int *)ptr;
	}

	T & value()
	{
		size_t ptr = (size_t)&triplet.value();

		return (T &)*(T *)ptr;
	}
};

/* ! \brief Sparse Matrix implementation, that map over Eigen
 *
 * \tparam T Type of the sparse Matrix store on each row,colums
 * \tparam id_t type of id
 * \tparam impl implementation
 *
 */
template<typename T, typename id_t>
class SparseMatrix<T,id_t,Eigen::SparseMatrix<T,0,id_t>>
{
	Eigen::SparseMatrix<T,0,id_t> mat;

public:

	// Triplet implementation id
	typedef boost::mpl::int_<EIGEN_TRIPLET> triplet_impl;

	// Triplet type
	typedef triplet<T,EIGEN_TRIPLET> triplet_type;

	/*! \brief Create an empty Matrix
	 *
	 * \param N1 number of row
	 * \param N2 number of colums
	 *
	 */
	SparseMatrix(size_t N1, size_t N2)
	:mat(N1,N2)
	{
	}

	/*! \brief Create a Matrix from a set of triplet
	 *
	 * \param N1 number of row
	 * \param N2 number of colums
	 * \param trpl triplet set
	 *
	 */
	SparseMatrix(size_t N1, size_t N2, openfpm::vector<Eigen::Triplet<T,id_t>> & trpl)
	:mat(N1,N2)
	{
		mat.setFromTriplets(trpl.begin(),trpl.end());
	}

	/*! \brief Create an empty Matrix
	 *
	 */
	SparseMatrix()	{}

	/*! \brief Set the Matrix from triplets
	 *
	 * \param trpl Triplet array from where to set
	 *
	 */
	void setFromTriplets(openfpm::vector<triplet_type> & trpl)
	{
		mat.setFromTriplets(trpl.begin(),trpl.end());
	}

	/*! \brief Get the Eigen Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	const Eigen::SparseMatrix<T,0,id_t> & getMat() const
	{
		return mat;
	}

	/*! \brief Get the Eigen Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	Eigen::SparseMatrix<T,0,id_t> & getMat()
	{
		return mat;
	}

	/*! \brief Resize the Sparse Matrix
	 *
	 * \param row number for row
	 * \param col number of colums
	 *
	 */
	void resize(size_t row, size_t col)
	{
		mat.resize(row,col);
	}
};



#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_ */
