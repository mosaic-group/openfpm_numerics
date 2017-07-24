/*
 * Vector.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_

#include "config/config.h"
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

/*! \brief Sparse Matrix implementation stub object when OpenFPM is compiled with no linear algebra support
 *
 */
template<typename T,unsigned int impl DEFAULT_VECTOR >
class Vector
{
	//! stub
	T stub;

	//! stub
	int stub_i;

public:

	/*! \brief stub copy constructor
	 *
	 * \param v stub
	 *
	 */
	Vector(const Vector<T> & v)
	:stub(0),stub_i(0)
	{std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief stub copy constructor
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(Vector<T> && v)
	:stub(0),stub_i(0)
	{std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief stub constructor from number of rows
	 *
	 * \param n stub
	 *
	 */
	Vector(size_t n)
	:stub(0),stub_i(0)
	{std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	//! stub default constructor
	Vector()
	:stub(0),stub_i(0)
	{std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief stub constructor
	 *
	 * \param n global number of row
	 * \param n_row_local local number of rows
	 *
	 */
	Vector(size_t n, size_t n_row_local)
	:stub(0),stub_i(0)
	{std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief stub resize
	 *
	 * \param row stub
	 * \param row_n stub
	 *
	 */
	void resize(size_t row, size_t row_n) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief stub insert
	 *
	 * \param i stub
	 * \param val stub
	 *
	 */
	void insert(size_t i, T val) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief stub insert
	 *
	 * \param i stub
	 *
	 * \return stub
	 *
	 */
	inline T & insert(size_t i) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return stub;}

	/*! \brief stub insert
	 *
	 * \param i stub
	 *
	 */
	inline const T & insert(size_t i) const {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub;}

	/*! \brief stub
	 *
	 * \param i stub
	 *
	 * \return stub
	 *
	 */
	const T & operator()(size_t i) const {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return stub;}

	/*! \brief stub
	 *
	 * \param i stub
	 *
	 * \return stub
	 *
	 */
	T & operator()(size_t i) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;return stub;}

	/*! \brief scatter
	 *
	 */
	void scatter() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}

	/*! \brief fromFile
	 *
	 * \param file stub
	 *
	 */
	void fromFile(std::string file) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl;}


	/*! \brief stub operator=
	 *
	 * \param v stub
	 *
	 * \return itself
	 *
	 */
	Vector<T> & operator=(const Vector<T> & v) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return *this;}

	/*! \brief stub operator=
	 *
	 * \param v stub
	 *
	 * \return itself
	 *
	 */
	Vector<T> & operator=(Vector<T> && v) {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return *this;}

	/*! \brief stub getVec
	 *
	 * \return stub
	 *
	 */
	int & getVec() {std::cerr << __FILE__ << ":" << __LINE__ << " Error in order to use this class you must compile OpenFPM with linear algebra support" << std::endl; return stub_i;}
};

#ifdef HAVE_EIGEN
#include "Vector_eigen.hpp"
#endif

#ifdef HAVE_PETSC
#include "Vector_petsc.hpp"
#endif

#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_HPP_ */
