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


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy from the Vector to the grid target
 * \note properties can be scalars or arrays of C++ primitives
 *
 * \tparam Eqs_sys System of equation information
 * \tparam S type of destination grid
 * \tparam Ev type of vector
 *
 */
template<typename Eqs_sys, typename S, typename Ev>
struct copy_ele
{
	//! destination grid element
	const grid_dist_key_dx<Eqs_sys::dims> key;

	//! destination grid
	S & grid_dst;

	//! source element inside the Eigen vector
	size_t lin_id;

	//! counter
	size_t prp_id;

	//! It is basically the number of grid points the Eigen vector has inside
	size_t gs_size;

	//! source Eigen vector
	const Ev & x;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key destination position
	 * \param grid_dst grid destination
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


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_UTIL_HPP_ */
