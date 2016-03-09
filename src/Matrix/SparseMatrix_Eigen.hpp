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
#include "VCluster.hpp"

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
public:

	//! Triplet implementation id
	typedef boost::mpl::int_<EIGEN_TRIPLET> triplet_impl;

	//! Triplet type
	typedef triplet<T,EIGEN_TRIPLET> triplet_type;

private:

	Eigen::SparseMatrix<T,0,id_t> mat;
	openfpm::vector<triplet_type> trpl;
	openfpm::vector<triplet_type> trpl_recv;

	/*! \brief Assemble the matrix
	 *
	 * \param trpl Matrix in form of triplets
	 * \param mat Matrix to assemble
	 *
	 */
	void assemble()
	{
		Vcluster & vcl = *global_v_cluster;

		////// On Master and only here
		// we assemble the Matrix from the collected data
		if (vcl.getProcessingUnits() != 1)
		{
			collect();
			// only master assemble the Matrix
			if (vcl.getProcessUnitID() == 0)
				mat.setFromTriplets(trpl_recv.begin(),trpl_recv.end());
		}
		else
			mat.setFromTriplets(trpl.begin(),trpl.end());
	}

	/*! \brief Here we collect the full matrix on master
	 *
	 */
	void collect()
	{
		Vcluster & vcl = *global_v_cluster;

		trpl_recv.clear();

		// here we collect all the triplet in one array on the root node
		vcl.SGather(trpl,trpl_recv,0);
	}

public:



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
	SparseMatrix(size_t N1, size_t N2, const openfpm::vector<Eigen::Triplet<T,id_t>> & trpl)
	:mat(N1,N2)
	{
		this->trpl = trpl;
	}

	/*! \brief Create an empty Matrix
	 *
	 */
	SparseMatrix()	{}

	/*! \brief Get the Matrix triplets bugger
	 *
	 * It return a buffer that can be filled with triplets
	 *
	 * \return triplet buffer
	 *
	 */
	openfpm::vector<triplet_type> & getMatrixTriplets()
	{
		return this->trpl;
	}

	/*! \brief Get the Eigen Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	const Eigen::SparseMatrix<T,0,id_t> & getMat() const
	{
		// Here we collect the information on master
		assemble();

		return mat;
	}

	/*! \brief Get the Eigen Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	Eigen::SparseMatrix<T,0,id_t> & getMat()
	{
		assemble();

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
