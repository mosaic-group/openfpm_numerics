/*
 * Vector.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_

#include "util/linalgebra_lib.hpp"

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


// Eigen::Matrix<T, Eigen::Dynamic, 1>

#ifdef HAVE_EIGEN
#include <Eigen/Sparse>
#define DEFAULT_VECTOR  = EIGEN_BASE
#else
#define DEFAULT_VECTOR = 0
#endif

/* ! \brief Sparse Matrix implementation stub object when OpenFPM is compiled with no linear algebra support
 *
 */
template<typename T,unsigned int impl DEFAULT_VECTOR >
class Vector
{
	T stub;
	int stub_i;

public:

	Vector(const Vector<T> & v) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	Vector(const Vector<T> && v) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	Vector(size_t n) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	Vector() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	Vector(size_t n, size_t n_row_local) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	void resize(size_t row) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	void insert(size_t i, T val) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	inline T & insert(size_t i) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return stub;}
	inline const T & insert(size_t i) const {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub;}
	const T & operator()(size_t i) const {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return stub;}
	T & operator()(size_t i) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return stub;}
	void scatter() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	void fromFile(std::string file) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}
	Vector<T> & operator=(const Vector<T> & v) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return *this;}
	Vector<T> & operator=(const Vector<T> && v) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return *this;}
	int & getVec() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_i;}
};

#ifdef HAVE_EIGEN
#include "Vector_eigen.hpp"
#endif

#ifdef HAVE_PETSC
#include "Vector_petsc.hpp"
#endif

#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_ */
