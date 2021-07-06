/*
 * common.hpp
 *
 *  Created on: Oct 11, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_UTIL_COMMON_HPP_
#define OPENFPM_NUMERICS_SRC_UTIL_COMMON_HPP_

#define STAGGERED_GRID 1
#define NORMAL_GRID 0

#include "util/common.hpp"

template<typename T, typename Sfinae = void>
struct has_grid_type: std::false_type {};


/*! \brief has_grid_type check if T has defined the member grid_type
 *
 * # Define structures with information inside
 * \snippet common_test.hpp Define structures
 * # Usage on the defined structures
 * \snippet common_test.hpp Check has grid_type
 *
 */
template<typename T>
struct has_grid_type<T, typename Void<decltype( T::grid_type )>::type> : std::true_type
{};

/*! \brief is_grid_staggered analyse T if it has a property grid_type defined and indicate that
 *         the grid is staggered
 *
 * This struct define a method value that return true if the grid is staggered, false otherwise
 *
 * # Define structures with information inside
 * \snippet common_test.hpp Define structures
 * # Usage on the defined structures
 * \snippet common_test.hpp Check grid_type staggered
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
		return T::grid_type == STAGGERED_GRID;
	};
};

/*! \brief is_grid_staggered analyse T if it has a property that define the type of grid
 *
 * This struct define a method value that return id T define a staggered or normal grid
 *
 * ### Structures that define a grid type
 * \snippet util_num_unit_tests.hpp grid_type_struct
 * ### Usage on the defined structures
 * \snippet util_num_unit_tests.hpp Usage of is_grid_staggered
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
