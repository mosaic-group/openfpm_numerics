/*
 * Matrix.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_HPP_

#include "config/config.h"
#include "util/linalgebra_lib.hpp"

#ifdef HAVE_EIGEN
#include <Eigen/Sparse>
#define DEFAULT_MATRIX  = EIGEN_BASE
#else
#define DEFAULT_MATRIX = 0
#endif

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
template<typename T, int impl> struct triplet
{
	long int i;
	long int j;
	T val;

	triplet(long int i, long int j, T val)
	{
		row() = i;
		col() = j;
		value() = val;
	}

	triplet()
	{
	}

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

/*! \brief Sparse Matrix implementation
 *
 * \tparam T Type of the sparse Matrix store on each row,colums
 * \tparam id_t type of id
 * \tparam Mi implementation
 *
 */
template<typename T,typename id_t ,unsigned int Mi DEFAULT_MATRIX>
class SparseMatrix
{
public:

	//! Triplet implementation id
	typedef boost::mpl::int_<-1> triplet_impl;

	//! Triplet type
	typedef triplet<T,-1> triplet_type;

	openfpm::vector<triplet_type> stub_vt;
	int stub_i;
	T stub_t;

public:


	SparseMatrix(size_t N1, size_t N2) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	SparseMatrix(size_t N1, size_t N2, size_t loc) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	SparseMatrix()	{std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	openfpm::vector<triplet_type> & getMatrixTriplets() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_vt;}
	const int & getMat() const {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_i;}
	int & getMat() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_i;}
	void resize(size_t row, size_t col, size_t row_n, size_t col_n) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	T operator()(id_t i, id_t j) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_t;}
	bool save(const std::string & file) const {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return true;}
	bool load(const std::string & file) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return false;}
	T getValue(size_t r, size_t c) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_i;}
};

#ifdef HAVE_EIGEN
#include "SparseMatrix_Eigen.hpp"
#endif

#ifdef HAVE_PETSC
#include "SparseMatrix_petsc.hpp"
#endif

#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_HPP_ */
