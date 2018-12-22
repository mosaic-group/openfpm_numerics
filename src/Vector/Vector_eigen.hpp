/*
 * Vector_eigen.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_

#include <type_traits>
#include "util/mul_array_extents.hpp"
#include <fstream>
#include "Grid/staggered_dist_grid.hpp"
#include "Space/Ghost.hpp"
#include "FiniteDifference/util/common.hpp"
#include <boost/mpl/vector_c.hpp>
#include <unordered_map>
#include "Vector_util.hpp"

#define EIGEN_RVAL 1

/*! \brief It store one row value of a vector
 *
 * Given a row, store a value
 *
 *
 */
template<typename T>
class rval<T,EIGEN_RVAL>
{
	//! row
	long int r;

	//! value
	T val;

public:

	/*! \brief Return the row index
	 *
	 * \return a reference to the row index
	 *
	 */
	long int & row()
	{
		return r;
	}

	/*! \brief Return the value
	 *
	 * \return a reference to the row value
	 *
	 */
	T & value()
	{
		return val;
	}

	/*! \brief Default constructor
	 *
	 */
	rval()
	:r(0)
	{}

	/*! \brief Constructor from row, colum and value
	 *
	 * \param i row
	 * \param val value
	 *
	 */
	rval(long int i, T val)
	{
		row() = i;
		value() = val;
	}

	/*! \brief Indicate that the structure has no pointer
	 *
	 * \return true
	 *
	 */
	static inline bool noPointers()
	{
		return true;
	}
};

template<typename T>
class Vector<T, EIGEN_BASE>
{
	//! Eigen vector
	mutable Eigen::Matrix<T, Eigen::Dynamic, 1> v;

	//! row val vector
	mutable openfpm::vector<rval<T,EIGEN_RVAL>> row_val;

	//! row val vector received
	mutable openfpm::vector<rval<T,EIGEN_RVAL>> row_val_recv;

	//! global to local map
	mutable std::unordered_map<size_t,size_t> map;

	//! invalid
	T invalid;

	//! Processors from where we gather
	mutable openfpm::vector<size_t> prc;

	//size of each chunk
	mutable openfpm::vector<size_t> sz;


	/*! \brief Here we collect the full vector on master
	 *
	 */
	void collect() const
	{
		Vcluster<> & vcl = create_vcluster();

		row_val_recv.clear();

		// here we collect all the triplet in one array on the root node
		vcl.SGather(row_val,row_val_recv,prc,sz,0);

		if (vcl.getProcessUnitID() != 0)
			row_val.resize(0);
		else
			row_val.swap(row_val_recv);

		build_map();
	}

	/*! \brief Set the Eigen internal vector
	 *
	 *
	 */
	void setEigen() const
	{
		// set the vector

		for (size_t i = 0 ; i < row_val.size() ; i++)
			v[row_val.get(i).row()] = row_val.get(i).value();
	}

	/*! \brief Build the map
	 *
	 *
	 */
	void build_map() const
	{
		map.clear();

		for (size_t i = 0 ; i < row_val.size() ; i++)
			map[row_val.get(i).row()] = i;
	}

public:

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(const Vector<T> & v)
	:invalid(0)
	{
		this->operator=(v);
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(const Vector<T> && v)
	:invalid(0)
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
		resize(n,0);
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
	 * \param l_row unused
	 *
	 */
	void resize(size_t row, size_t l_row)
	{
		v.resize(row);
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

		row_val.last().row() = i;
		row_val.last().value() = val;
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

		row_val.last().row() = i;
		return row_val.last().value();
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
			return row_val.get(it->second).value();

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
			return row_val.get(it->second).value();

		return insert(i);
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	const Eigen::Matrix<T, Eigen::Dynamic, 1> & getVec() const
	{
		collect();
		setEigen();

		return v;
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	Eigen::Matrix<T, Eigen::Dynamic, 1> & getVec()
	{
		collect();
		setEigen();

		return v;
	}

	/*! \brief Scatter the vector information to the other processors
	 *
	 * Eigen does not have a real parallel vector, so in order to work we have to scatter
	 * the vector from one processor to the other
	 *
	 *
	 */
	void scatter()
	{
		row_val_recv.clear();
		Vcluster<> & vcl = create_vcluster();

		vcl.SScatter(row_val,row_val_recv,prc,sz,0);

		// if we do not receive anything a previous collect has not been performed
		// and so nothing is scattered
		if (row_val_recv.size() != 0)
		{
			row_val.clear();
			row_val.add(row_val_recv);
			build_map();
		}
	}

	/*! \brief Load from file
	 *
	 *
	 *
	 */
	void fromFile(std::string file)
	{
		std::ifstream inputf;
		inputf.open(file);

		for (size_t i = 0 ; i < v.size() ; i++)
			inputf >> v(i);

		inputf.close();

	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector<T> & operator=(const Vector<T> & v)
	{
		prc = v.prc;
		sz = v.sz;
		map = v.map;
		row_val = v.row_val;

		return *this;
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector<T> & operator=(const Vector<T> && v)
	{
		prc = v.prc;
		sz = v.sz;
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
