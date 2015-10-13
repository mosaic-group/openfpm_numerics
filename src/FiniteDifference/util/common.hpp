/*
 * common.hpp
 *
 *  Created on: Oct 11, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_UTIL_COMMON_HPP_
#define OPENFPM_NUMERICS_SRC_UTIL_COMMON_HPP_


template<typename T, typename Sfinae = void>
struct has_grid_type: std::false_type {};


/*! \brief has_noPointers check if a type has defined a
 * method called noPointers
 *
 * ### Example
 *
 * \snippet util.hpp Check no pointers
 *
 * return true if T::noPointers() is a valid expression (function pointers)
 * and produce a defined type
 *
 */
template<typename T>
struct has_grid_type<T, typename Void<decltype( T::grid_type )>::type> : std::true_type
{};

/*! \brief is_grid_staggered analyse T if it has a property that define the type of grid
 *
 * This struct define a method value that return id T define a staggered or normal grid
 *
 */
template<typename T, bool has_gt = has_grid_type<T>::type::value >
struct is_grid_staggered
{
	/*! \brief Return true if the grid is staggered
	 *
	 * If T define a staggered grid return true, if define a NORMAL grid or does not
	 * define any grid type return false
	 *
	 * \return true in case the grid type is staggered
	 *
	 */
	static bool value()
	{
		return T::grid_type;
	};
};

/*! \brief is_grid_staggered analyse T if it has a property that define the type of grid
 *
 * This struct define a method value that return id T define a staggered or normal grid
 *
 */
template<typename T>
struct is_grid_staggered<T,false>
{

	/*! \brief Return true if the grid is staggered
	 *
	 * If T define a staggered grid return true, if define a NORMAL grid or does not
	 * define any grid type return false
	 *
	 * \return true in case the grid type is staggered
	 *
	 */
	static bool value()
	{
		return false;
	};
};

#endif /* OPENFPM_NUMERICS_SRC_UTIL_COMMON_HPP_ */
