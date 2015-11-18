/*
 * util_num.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_UTIL_UTIL_NUM_HPP_
#define OPENFPM_NUMERICS_SRC_UTIL_UTIL_NUM_HPP_

#include "util/common.hpp"
#include "grid_dist_testing.hpp"

template<typename T, typename Sfinae = void>
struct is_const_field: std::false_type {};


/*! \brief is_constant check if a type define a constant field
 *
 * ### Example
 * \snippet util_num.hpp Constant fields struct definition
 * \snippet util_num.hpp Usage of is_const_field
 *
 * return true if T::constant_field exist, it does not matter which type is
 *
 */
template<typename T>
struct is_const_field<T, typename Void<typename T::const_field>::type> : std::true_type
{};


template<typename T, typename Sfinae = void>
struct has_val_pos: std::false_type {};


/*! \brief has_attributes check if a type has defined an
 * internal structure with attributes
 *
 * ### Example
 * \snippet util.hpp Declaration of struct with attributes and without
 * \snippet util.hpp Check has_attributes
 *
 * return true if T::attributes::name[0] is a valid expression
 * and produce a defined type
 *
 */
template<typename T>
struct has_val_pos<T, typename Void<typename T::with_position >::type> : std::true_type
{};



template<typename T, typename Sfinae = void>
struct is_testing: std::false_type {};


/*! \brief has_attributes check if a type has defined an
 * internal structure with attributes
 *
 * ### Example
 * \snippet util_num_unit_tests.hpp Is on test struct definition
 * \snippet util_num_unit_tests.hpp Usage of is_testing
 *
 * return true if T::attributes::name[0] is a valid expression
 * and produce a defined type
 *
 */
template<typename T>
struct is_testing<T, typename Void<typename T::testing >::type> : std::true_type
{};


/*! \brief In case of testing return a stub grid
 *
 * ### Example
 * \snippet util_num_unit_tests.hpp Is on test struct definition
 * \snippet util_num_unit_tests.hpp Usage of is_testing
 *
 * return true if T::attributes::name[0] is a valid expression
 * and produce a defined type
 *
 */
template<typename T, unsigned int dims, typename stype, typename decomposition, bool stub=is_testing<T>::value>
struct stub_or_real
{
	typedef grid_dist_testing<dims> type;
};

template<typename T, unsigned int dims, typename stype, typename decomposition>
struct stub_or_real<T,dims,stype,decomposition,false>
{
	typedef grid_dist_id<dims,stype,scalar<size_t>,decomposition> type;
};

#endif /* OPENFPM_NUMERICS_SRC_UTIL_UTIL_NUM_HPP_ */
