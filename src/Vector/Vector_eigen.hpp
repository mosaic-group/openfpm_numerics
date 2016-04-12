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
#include "Vector_eigen_util.hpp"
#include "Grid/staggered_dist_grid.hpp"
#include "Space/Ghost.hpp"
#include "FiniteDifference/util/common.hpp"
#include <boost/mpl/vector_c.hpp>
#include <unordered_map>

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
	// row
	long int r;

	// value
	T val;

public:

	// Get the row
	long int & row()
	{
		return r;
	}

	// Get the value
	T & value()
	{
		return val;
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
		row() = i;
		value() = val;
	}
};

template<typename T>
class Vector<T,Eigen::Matrix<T, Eigen::Dynamic, 1>>
{
	mutable Eigen::Matrix<T, Eigen::Dynamic, 1> v;

	// row value vector
	mutable openfpm::vector<rval<T,EIGEN_RVAL>> row_val;
	mutable openfpm::vector<rval<T,EIGEN_RVAL>> row_val_recv;

	// global to local map
	mutable std::unordered_map<size_t,size_t> map;

	// invalid
	T invalid;

	// Processors from where we gather
	mutable openfpm::vector<size_t> prc;

	//size of each chunk
	mutable openfpm::vector<size_t> sz;

	/*! \brief Copy in a staggered grid
	 *
	 *
	 */
	template<typename scheme, typename Grid_dst ,unsigned int ... pos> void copy_staggered_to_staggered(scheme & sc, Grid_dst & g_dst)
	{
		// get the map
		const auto & g_map = sc.getMap();

		// check that g_dst is staggered
		if (is_grid_staggered<Grid_dst>::value() == false)
			std::cerr << __FILE__ << ":" << __LINE__ << " The destination grid must be staggered " << std::endl;

#ifdef SE_CLASS1

		if (g_map.getLocalDomainSize() != g_dst.getLocalDomainSize())
			std::cerr << __FILE__ << ":" << __LINE__ << " The staggered and destination grid in size does not match " << std::endl;
#endif

		// sub-grid iterator over the grid map
		auto g_map_it = g_map.getDomainIterator();

		// Iterator over the destination grid
		auto g_dst_it = g_dst.getDomainIterator();

		while (g_map_it.isNext() == true)
		{
			typedef typename to_boost_vmpl<pos...>::type vid;
			typedef boost::mpl::size<vid> v_size;

			auto key_src = g_map_it.get();
			size_t lin_id = g_map.template get<0>(key_src);

			// destination point
			auto key_dst = g_dst_it.get();

			// Transform this id into an id for the Eigen vector

			copy_ele<typename scheme::Sys_eqs_typ,Grid_dst,typename std::remove_reference<decltype(*this)>::type> cp(key_dst,g_dst,*this,lin_id,g_map.size());

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

			++g_map_it;
			++g_dst_it;
		}
	}


	/*! \brief Here we collect the full matrix on master
	 *
	 */
	void collect() const
	{
		Vcluster & vcl = *global_v_cluster;

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

	/*! \brief Copy the vector into the grid
	 *
	 * ## Copy from the vector to the destination grid
	 * \snippet eq_unit_tests.hpp
	 *
	 * \tparam scheme Discretization scheme
	 * \tparam Grid_dst type of the target grid
	 * \tparam pos target properties
	 *
	 * \param scheme Discretization scheme
	 * \param start point
	 * \param stop point
	 * \param g_dst Destination grid
	 *
	 */
	template<typename scheme, typename Grid_dst ,unsigned int ... pos> void copy(scheme & sc,const long int (& start)[scheme::Sys_eqs_typ::dims], const long int (& stop)[scheme::Sys_eqs_typ::dims], Grid_dst & g_dst)
	{
		if (is_grid_staggered<typename scheme::Sys_eqs_typ>::value())
		{
			if (g_dst.is_staggered() == true)
				copy_staggered_to_staggered<scheme,Grid_dst,pos...>(sc,g_dst);
			else
			{
				// Create a temporal staggered grid and copy the data there
				auto & g_map = sc.getMap();
				staggered_grid_dist<Grid_dst::dims,typename Grid_dst::stype,typename Grid_dst::value_type,typename Grid_dst::decomposition::extended_type, typename Grid_dst::memory_type, typename Grid_dst::device_grid_type> stg(g_dst,g_map.getDecomposition().getGhost(),sc.getPadding());
				stg.setDefaultStagPosition();
				copy_staggered_to_staggered<scheme,decltype(stg),pos...>(sc,stg);

				// sync the ghost and interpolate to the normal grid
				stg.template ghost_get<pos...>();
				stg.template to_normal<Grid_dst,pos...>(g_dst,sc.getPadding(),start,stop);
			}
		}
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
		Vcluster & vcl = *global_v_cluster;

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
