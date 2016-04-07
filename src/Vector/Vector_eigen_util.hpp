/*
 * Vector_util.hpp
 *
 *  Created on: Dec 7, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_UTIL_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_UTIL_HPP_

#include "Grid/grid_dist_key.hpp"
#include "Space/Shape/HyperCube.hpp"

/*!	\brief Copy scalar elements
 *
 * \tparam copy_type Type that should be copied
 * \tparam T property id to copy
 * \tparam Ev Type of source the Vector
 * \tparam sa dimensionality of the array 0 is a scalar
 *
 */
template<typename copy_type, typename T, typename Ev, typename Eqs_sys, int sa>
struct copy_ele_sca_array
{
	template<typename Grid> inline static void copy(Grid & grid_dst, const grid_dist_key_dx<Eqs_sys::dims> & key, const Ev & x,size_t lin_id, size_t base_id, size_t gs_size)
	{
		grid_dst.template get<T::value>(key) = x(lin_id * Eqs_sys::nvar + base_id);
	}
};

/*! \brief Copy 1D array elements
 *
 * spacialization in case of 1D array
 *
 * \tparam copy_type Type that should be copied
 * \tparam T property id to copy
 * \tparam Ev Type of source the Vector
 *
 */
template<typename copy_type, typename T, typename Ev, typename Eqs_sys>
struct copy_ele_sca_array<copy_type,T,Ev,Eqs_sys,1>
{
	template<typename Grid> inline static void copy(Grid & grid_dst, const grid_dist_key_dx<Eqs_sys::dims> & key, const Ev & x,size_t lin_id, size_t base_id, size_t gs_size)
	{
		for (size_t i = 0 ; i < std::extent<copy_type>::value ; i++)
		{
			grid_dst.template get<T::value>(key)[i] = x(lin_id * Eqs_sys::nvar + base_id + i);
		}
	}
};

/*!	\brief Add scalar elements
 *
 * \tparam copy_type Type that should be copied
 * \tparam T property id to copy
 * \tparam Ev Type of source the Vector
 * \tparam sa dimensionality of the array 0 is a scalar
 *
 */
template<typename copy_type, typename T, typename Ev, typename scheme, int sa>
struct interp_ele_sca_array
{
	template<typename Grid> inline static void interp(Grid & grid_dst, const grid_dist_key_dx<scheme::Sys_eqs_typ::dims> & key, const Ev & x, grid_dist_key_dx<scheme::Sys_eqs_typ::dims> key_src, const openfpm::vector<std::vector<comb<scheme::Sys_eqs_typ::dims>>> & interp_pos, const typename scheme::g_map_type & g_map,size_t base_id)
	{
		copy_type division = 0.0;

		for (size_t i = 0 ; i < interp_pos.get(0).size() ; i++)
		{
				auto key_m = key_src.move(interp_pos.get(0)[i]);
				size_t lin_id = g_map.template get<0>(key_m);

				grid_dst.template get<T::value>(key) += x(lin_id * scheme::Sys_eqs_typ::nvar + base_id);

				division += 1.0;
		}
		grid_dst.template get<T::value>(key) /= division;
	}
};

/*! \brief Add 1D array elements
 *
 * spacialization in case of 1D array
 *
 * \tparam copy_type Type that should be copied
 * \tparam T property id to copy
 * \tparam Ev Type of source the Vector
 *
 */
template<typename copy_type, typename T, typename Ev, typename scheme>
struct interp_ele_sca_array<copy_type,T,Ev,scheme,1>
{
	template<typename Grid> inline static void interp(Grid & grid_dst, const grid_dist_key_dx<scheme::Sys_eqs_typ::dims> & key, const Ev & x, grid_dist_key_dx<scheme::Sys_eqs_typ::dims> key_src, const openfpm::vector<std::vector<comb<scheme::Sys_eqs_typ::dims>>> & interp_pos, const typename scheme::g_map_type & g_map,size_t base_id)
	{
		typename std::remove_all_extents<copy_type>::type division;

		for (size_t j = 0 ; j < std::extent<copy_type>::value ; j++)
		{
			division = 0.0;
			for (size_t i = 0 ; i < interp_pos.get(j).size() ; i++)
			{
					auto key_m = key_src.move(interp_pos.get(j)[i]);
					size_t lin_id = g_map.template get<0>(key_m);

					grid_dst.template get<T::value>(key)[j] += x(lin_id * scheme::Sys_eqs_typ::nvar + base_id);

					division += 1.0;
			}
			grid_dst.template get<T::value>(key)[j] /= division;
			base_id++;
		}
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to calculate the interpolation points for each
 * property in a staggered grid
 *
 * \tparam dim Dimensionality
 * \tparam v_prp_id vector of properties id
 * \tparam v_prp_type vector with the properties
 *
 */
template<unsigned int dim, typename v_prp_id, typename v_prp_type>
struct interp_points
{
	// number of properties we are processing
	typedef boost::mpl::size<v_prp_id> v_size;

	// interpolation points for each property
	openfpm::vector<std::vector<comb<dim>>> (& interp_pts)[v_size::value];

	// staggered position for each property
	const openfpm::vector<comb<dim>> (&stag_pos)[v_size::value];

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param inter_pts array that for each property contain the interpolation points for each components
	 * \param staggered position for each property and components
	 *
	 */
	inline interp_points(openfpm::vector<std::vector<comb<dim>>> (& interp_pts)[v_size::value],const openfpm::vector<comb<dim>> (&stag_pos)[v_size::value])
	:interp_pts(interp_pts),stag_pos(stag_pos){};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at_c<v_prp_type,T::value>::type prp_type;

		interp_pts[T::value].resize(stag_pos[T::value].size());

		for (size_t i = 0 ; i < stag_pos[T::value].size() ; i++)
		{
			// Create the interpolation points
			interp_pts[T::value].get(i) = SubHyperCube<dim,dim - std::rank<prp_type>::value>::getCombinations_R(stag_pos[T::value].get(i),0);

			// interp_point are -1,0,1, map the -1 to 0 and 1 to -1
			for (size_t j = 0 ; j < interp_pts[T::value].get(i).size() ; j++)
			{
				for (size_t k = 0 ; k < dim ; k++)
					interp_pts[T::value].get(i)[j].getComb()[k] = - ((interp_pts[T::value].get(i)[j].getComb()[k] == -1)?0:interp_pts[T::value].get(i)[j].getComb()[k]);
			}
		}
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy from the Vector_eigen to the grid target
 * in a generic way for a generic object T
 * with variables number of properties scalar or array of C++ primitives
 *
 * \tparam dim Dimensionality
 * \tparam S type of grid
 * \tparam Ev Eigen vector type
 *
 */
template<typename Eqs_sys, typename S, typename Ev>
struct copy_ele
{
	//! destination element inside the grid
	const grid_dist_key_dx<Eqs_sys::dims> key;

	//! destination grid
	S & grid_dst;

	//! source element inside the Eigen vector
	size_t lin_id;

	//! property id inside the Eigen vector, each property processed increase such id
	size_t prp_id;

	//! number of the elements in the Eigen vector divided the number of variables
	//! It is basically the number of grid points the Eigen vector has inside
	size_t gs_size;

	//! source Eigen vector
	const Ev & x;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_dst grid we are updating
	 * \param v Source Eigen vector
	 *
	 */
	inline copy_ele(const grid_dist_key_dx<Eqs_sys::dims> & key, S & grid_dst, const Ev & x, size_t lin_id, size_t gs_size)
	:key(key),grid_dst(grid_dst),lin_id(lin_id),prp_id(0),gs_size(gs_size),x(x){};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_ele(grid_dist_key_dx<Eqs_sys::dims> & key, S & grid_dst, Ev && x)
	:key(key),grid_dst(grid_dst),x(x)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at_c<typename S::value_type::type,T::value>::type copy_type;

		copy_ele_sca_array<copy_type,T,Ev,Eqs_sys,std::is_array<copy_type>::value>::copy(grid_dst,key,x,lin_id,prp_id,gs_size);
		prp_id += array_extents<copy_type>::mul();
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to interpolate from the staggered Vector_eigen to the grid target
 * in a generic way for a generic object T
 * with variables number of properties scalar or array of C++ primitives
 *
 * \tparam scheme Discretization scheme
 * \tparam S type of grid
 * \tparam Ev Eigen vector type
 *
 */
template<typename scheme, typename S, typename Ev, typename base_type, unsigned int nst_pos>
struct interp_ele
{
	//! destination element inside the grid
	const grid_dist_key_dx<scheme::Sys_eqs_typ::dims> key;

	//! destination grid
	S & grid_dst;

	//! source element inside the Eigen vector
	grid_dist_key_dx<scheme::Sys_eqs_typ::dims> key_src;

	//! g_map
	const typename scheme::g_map_type & g_map;

	//! For each properties [] for each components (openfpm::vector) interpolants points positions (std::vector<comb>)
	openfpm::vector<std::vector<comb<scheme::Sys_eqs_typ::dims>>> (&interp_pos)[nst_pos];

	//! property id inside the Eigen vector, each property processed increase such id
	size_t prp_id;

	//! source Eigen vector
	const Ev & x;

	//! Value type of the vector
	base_type division;

	/*! \brief constructor
	 *
	 * It define the interpolation parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_dst grid we are updating
	 * \param v Source Eigen vector
	 *
	 */
	inline interp_ele(const grid_dist_key_dx<scheme::Sys_eqs_typ::dims> & key, S & grid_dst, const Ev & x, const grid_dist_key_dx<scheme::Sys_eqs_typ::dims> & key_src, const typename scheme::g_map_type & g_map, openfpm::vector<std::vector<comb<scheme::Sys_eqs_typ::dims>>> (&interp_pos)[nst_pos])
	:key(key),grid_dst(grid_dst),key_src(key_src),g_map(g_map),interp_pos(interp_pos),prp_id(0),x(x){};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline interp_ele(grid_dist_key_dx<scheme::Sys_eqs_typ::dims> & key, S & grid_dst, Ev && x)
	:key(key),grid_dst(grid_dst),x(x),interp_pos(interp_pos)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at_c<typename S::value_type::type,T::value>::type copy_type;

		interp_ele_sca_array<copy_type,T,Ev,scheme,std::rank<copy_type>::value>::interp(grid_dst,key,x,key_src,interp_pos[T::value],g_map,prp_id);

		prp_id += array_extents<copy_type>::mul();
	}
};

#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_UTIL_HPP_ */
