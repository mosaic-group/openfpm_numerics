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
#include <boost/mpl/vector_c.hpp>


template<typename T>
class Vector<T,Eigen::Matrix<T, Eigen::Dynamic, 1>>
{
	Eigen::Matrix<T, Eigen::Dynamic, 1> v;

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

		// fill the staggered position vector for each property
		stag_set_position<scheme::Sys_eqs_typ::dims,typename Grid_dst::value_type::type> ssp(stag_pos);

		typedef boost::mpl::vector_c<unsigned int,pos ... > v_pos_type;
		boost::mpl::for_each_ref<v_pos_type>(ssp);

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

			interp_ele<scheme,Grid_dst,Eigen::Matrix<T, Eigen::Dynamic, 1>,T,sizeof...(pos)> cp(key_dst,g_dst,v,key_src,g_map,stag_pos);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

			++g_map_it;
			++g_dst_it;
		}
	}

public:

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
	T & get(size_t i)
	{
		return v(i);
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	const Eigen::Matrix<T, Eigen::Dynamic, 1> & getVec() const
	{
		return v;
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	Eigen::Matrix<T, Eigen::Dynamic, 1> & getVec()
	{
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
};


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_ */
