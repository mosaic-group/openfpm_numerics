/*
 * Vector_petsc.hpp
 *
 *  Created on: Apr 29, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_PETSC_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_PETSC_HPP_

#include "Vector/map_vector.hpp"
#include <boost/mpl/int.hpp>
#include <petscvec.h>
#include "util/petsc_util.hpp"

#define PETSC_RVAL 2

/*! \brief It store one row value of a vector
 *
 * Given a row, store a value
 *
 *
 */
template<typename T>
class rval<T,PETSC_RVAL>
{
public:

	//! boost fusion that store the point
	typedef boost::fusion::vector<PetscInt,T> type;

	//! structure that store the data of the point
	type data;

	//! Property id of the point
	static const unsigned int row = 0;
	static const unsigned int value = 1;
	static const unsigned int max_prop = 2;

	// Get the row
	long int & rw()
	{
		return boost::fusion::at_c<row>(data);
	}

	// Get the value
	T & val()
	{
		return boost::fusion::at_c<value>(data);
	}

	/*! \brief Default constructor
	 *
	 */
	rval()	{}

	/*! \brief Constructor from row, colum and value
	 *
	 * \param i row
	 * \param val value
	 *
	 */
	rval(long int i, T val)
	{
		rw() = i;
		val() = val;
	}
};

constexpr unsigned int row_id = 0;
constexpr unsigned int val_id = 1;

template<typename T>
class Vector<T,Vec>
{
	// n_row
	size_t n_row;

	// n_row_local
	size_t n_row_local;

	// Mutable vector
	mutable Vec v;

	// Mutable row value vector
	openfpm::vector<rval<T,PETSC_RVAL>,HeapMemory,typename memory_traits_inte<rval<T,PETSC_RVAL>>::type > row_val;

	// Global to local map
	mutable std::unordered_map<size_t,size_t> map;

	// invalid
	T invalid;

	/*! \brief Set the Eigen internal vector
	 *
	 *
	 */
	void setPetsc() const
	{
		// Create the vector
		PETSC_SAFE_CALL(VecCreate(PETSC_COMM_WORLD,&v));
		PETSC_SAFE_CALL(VecSetSizes(v,n_row_local,n_row));

		// set the vector

		PETSC_SAFE_CALL(VecSetValues(v,n_row_local,&row_val.template get<row_id>(0),&row_val.template get<val_id>(0),INSERT_VALUES))

//		for (size_t i = 0 ; i < row_val.size() ; i++)
//			v[row_val.get(i).row()] = row_val.get(i).value();
	}

public:

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(const Vector<T> & v)
	{
		this->operator=(v);
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(const Vector<T> && v)
	{
		this->operator=(v);
	}

	/*! \brief Create a vector with n elements
	 *
	 * \param n number of elements in the vector
	 *
	 */
	Vector(size_t n)
	{
		resize(n);
	}

	/*! \brief Create a vector with 0 elements
	 *
	 */
	Vector()
	{
	}

	/*! \brief Resize the Vector
	 *
	 * \param row numbers of row
	 *
	 */
	void resize(size_t row)
	{
		n_row = row;
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 * \param val value
	 *
	 */
	void insert(size_t i, T val)
	{
		row_val.add();

		// Map
		map[i] = row_val.size()-1;

		row_val.last().template get<row_id>() = i;
		row_val.last().template get<val_id>() = val;
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	inline T & insert(size_t i)
	{
		row_val.add();

		// Map
		map[i] = row_val.size()-1;

		row_val.last().template get<row_id>() = i;
		return row_val.last().template get<val_id>();
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	inline const T & insert(size_t i) const
	{
		row_val.add();

		// Map
		map[i] = row_val.size()-1;

		row_val.last().row() = i;
		return row_val.last().value();
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \warning The element must exist
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	const T & operator()(size_t i) const
	{
		// Search if exist

		std::unordered_map<size_t,size_t>::iterator it = map.find(i);

		if ( it != map.end() )
			return row_val.template get<val_id>(it->second);

		return insert(i);
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \warning The element must exist
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	T & operator()(size_t i)
	{
		// Search if exist

		std::unordered_map<size_t,size_t>::iterator it = map.find(i);

		if ( it != map.end() )
			return row_val.template get<val_id>(it->second);

		return insert(i);
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	const Vec & getVec() const
	{
		setPetsc();

		return v;
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	Vec & getVec()
	{
		setPetsc();

		return v;
	}


	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector<T,Vec> & operator=(const Vector<T,Vec> & v)
	{
		map = v.map;
		row_val = v.row_val;

		return *this;
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector<T,Vec> & operator=(const Vector<T,Vec> && v)
	{
		map = v.map;
		row_val = v.row_val;

		return *this;
	}

	/*! \brief Copy the vector (it is used for special purpose)
	 *
	 * \warning v MUST contain at least all the elements of the vector
	 *
	 * \param v base eigen vector to copy
	 *
	 */
	Vector<T> & operator=(Eigen::Matrix<T, Eigen::Dynamic, 1> & v)
	{
		for (size_t i = 0 ; i < row_val.size() ; i++)
			row_val.get(i).value() = v(row_val.get(i).row());

		return *this;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_ */

