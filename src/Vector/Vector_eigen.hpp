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
#include "Grid/staggered_dist_grid_util.hpp"
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

	/*! \brief Check that the size of the iterators match
	 *
	 * It check the the boxes that the sub iterator defines has same dimensions, for example
	 * if the first sub-iterator, iterate from (1,1) to (5,3) and the second from (2,2) to (6,4)
	 * they match (2,2) to (4,6) they do not match
	 *
	 * \tparam Grid_map type of the map grid
	 * \tparam Grid_dst type of the destination grid
	 *
	 * \param it1 Iterator1
	 * \param it2 Iterator2
	 *
	 * \return true if they match
	 *
	 */
	template<typename Eqs_sys, typename it1_type, typename it2_type> bool checkIterator(const it1_type & it1, const it2_type & it2)
	{
#ifdef SE_CLASS1

		grid_key_dx<Eqs_sys::dims> it1_k = it1.getStop() - it1.getStart();
		grid_key_dx<Eqs_sys::dims> it2_k = it2.getStop() - it2.getStart();

		for (size_t i = 0 ; i < Eqs_sys::dims ; i++)
		{
			if (it1_k.get(i) !=  it2_k.get(i))
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " error src iterator and destination iterator does not match in size\n";
				return false;
			}
		}

		return true;
#else

		return true;

#endif
	}


	/*! \brief Copy in a staggered grid
	 *
	 *
	 */
	template<typename scheme, typename Grid_dst ,unsigned int ... pos> void copy_staggered_to_staggered(scheme & sc,const long int (& start)[scheme::Sys_eqs_typ::dims], const long int (& stop)[scheme::Sys_eqs_typ::dims], Grid_dst & g_dst)
	{
		// get the map
		const auto & g_map = sc.getMap();

		// get the grid padding
		const Padding<scheme::Sys_eqs_typ::dims> & pd = sc.getPadding();

		// shift the start and stop by the padding
		grid_key_dx<scheme::Sys_eqs_typ::dims> start_k = grid_key_dx<scheme::Sys_eqs_typ::dims>(start) + pd.getKP1();
		grid_key_dx<scheme::Sys_eqs_typ::dims> stop_k = grid_key_dx<scheme::Sys_eqs_typ::dims>(stop) + pd.getKP1();

		// sub-grid iterator over the grid map
		auto g_map_it = g_map.getSubDomainIterator(start_k,stop_k);

		// Iterator over the destination grid
		auto g_dst_it = g_dst.getDomainIterator();

		// Check that the 2 iterator has the same size
		checkIterator<typename scheme::Sys_eqs_typ,decltype(g_map_it),decltype(g_dst_it)>(g_map_it,g_dst_it);

		while (g_map_it.isNext() == true)
		{
			typedef typename to_boost_vmpl<pos...>::type vid;
			typedef boost::mpl::size<vid> v_size;

			auto key_src = g_map_it.get();
			size_t lin_id = g_map.template get<0>(key_src);

			// destination point
			auto key_dst = g_dst_it.get();

			// Transform this id into an id for the Eigen vector

			copy_ele<typename scheme::Sys_eqs_typ,Grid_dst,Eigen::Matrix<T, Eigen::Dynamic, 1>> cp(key_dst,g_dst,v,lin_id,g_map.size());

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

			++g_map_it;
			++g_dst_it;
		}
	}

	/*! \brief Copy in a normal grid
	 *
	 *
	 */
	template<typename scheme, typename Grid_dst ,unsigned int ... pos> void copy_staggered_to_normal(scheme & sc,const long int (& start)[scheme::Sys_eqs_typ::dims], const long int (& stop)[scheme::Sys_eqs_typ::dims], Grid_dst & g_dst)
	{
		// get the map
		const auto & g_map = sc.getMap();

		// get the grid padding
		const Padding<scheme::Sys_eqs_typ::dims> & pd = sc.getPadding();

		// set the staggered position of the property
		openfpm::vector<comb<scheme::Sys_eqs_typ::dims>> stag_pos[sizeof...(pos)];

		// interpolation points for each property
		openfpm::vector<std::vector<comb<scheme::Sys_eqs_typ::dims>>> interp_pos[sizeof...(pos)];

		// fill the staggered position vector for each property
		stag_set_position<scheme::Sys_eqs_typ::dims,typename Grid_dst::value_type::type> ssp(stag_pos);

		typedef boost::mpl::vector_c<unsigned int,pos ... > v_pos_type;
		boost::mpl::for_each_ref<v_pos_type>(ssp);

		interp_points<scheme::Sys_eqs_typ::dims,v_pos_type,typename Grid_dst::value_type::type> itp(interp_pos,stag_pos);
		boost::mpl::for_each_ref<v_pos_type>(itp);

		// shift the start and stop by the padding
		grid_key_dx<scheme::Sys_eqs_typ::dims> start_k = grid_key_dx<scheme::Sys_eqs_typ::dims>(start) + pd.getKP1();
		grid_key_dx<scheme::Sys_eqs_typ::dims> stop_k = grid_key_dx<scheme::Sys_eqs_typ::dims>(stop) + pd.getKP1();

		// sub-grid iterator over the grid map
		auto g_map_it = g_map.getSubDomainIterator(start_k,stop_k);

		// Iterator over the destination grid
		auto g_dst_it = g_dst.getDomainIterator();

		// Check that the 2 iterator has the same size
		checkIterator<typename scheme::Sys_eqs_typ,decltype(g_map_it),decltype(g_dst_it)>(g_map_it,g_dst_it);

		while (g_map_it.isNext() == true)
		{
			typedef typename to_boost_vmpl<pos...>::type vid;
			typedef boost::mpl::size<vid> v_size;

			auto key_src = g_map_it.get();

			// destination point
			auto key_dst = g_dst_it.get();

			// Transform this id into an id for the Eigen vector

			interp_ele<scheme,Grid_dst,Eigen::Matrix<T, Eigen::Dynamic, 1>,T,sizeof...(pos)> cp(key_dst,g_dst,v,key_src,g_map,interp_pos);

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
	 *
	 * \return reference to the element vector
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

		if ( it == map.end() )
		{
			// Does not exist

			std::cerr << __FILE__ << ":" << __LINE__ << " value does not exist " << std::endl;

			return invalid;
		}
		else
			return row_val.get(it->second).value();

		return invalid;
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
	 *
	 */
	template<typename scheme, typename Grid_dst ,unsigned int ... pos> void copy(scheme & sc,const long int (& start)[scheme::Sys_eqs_typ::dims], const long int (& stop)[scheme::Sys_eqs_typ::dims], Grid_dst & g_dst)
	{
		if (is_grid_staggered<typename scheme::Sys_eqs_typ>::value())
		{
			if (g_dst.is_staggered() == true)
				copy_staggered_to_staggered<scheme,Grid_dst,pos...>(sc,start,stop,g_dst);
			else
				copy_staggered_to_normal<scheme,Grid_dst,pos...>(sc,start,stop,g_dst);
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
