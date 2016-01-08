/*
 * sum.hpp
 *
 *  Created on: Oct 13, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_

#include <boost/mpl/vector.hpp>
#include "config.h"
#include <unordered_map>
#include "util/for_each_ref.hpp"

template<typename v_expr>
struct sum_functor_value
{
	//! Number of elements in the vector v_expr
	typedef boost::mpl::size<v_expr> size;

	//! Last element of sum
	typedef typename boost::mpl::at<v_expr,boost::mpl::int_<size::value-1> >::type last;

	//! sum functor
	std::unordered_map<long int,typename last::stype> & cols;

	//! Grid info
	const grid_sm<last::dims,void> & gs;

	// grid mapping
	const typename stub_or_real<last,last::dims,typename last::stype,typename last::b_grid::decomposition>::type & g_map;

	// grid position
	grid_dist_key_dx<last::dims> & kmap;

	//! position
	grid_key_dx<last::dims> & key;

	//! coefficent
	typename last::stype coeff;

	//! spacing
	typename last::stype (& spacing)[last::dims];


	/*! \brief constructor
	 *
	 */
	sum_functor_value(const typename stub_or_real<last,last::dims,typename last::stype,typename last::b_grid::decomposition>::type & g_map, grid_dist_key_dx<last::dims> & kmap, const grid_sm<last::dims,void> & gs, typename last::stype (& spacing)[last::dims], std::unordered_map<long int,typename last::stype> & cols, typename last::stype coeff)
	:cols(cols),gs(gs),g_map(g_map),kmap(kmap),key(key),coeff(coeff),spacing(spacing)
	{};

	//! It call this function for every expression in the sum
	template<typename T>
	void operator()(T& t) const
	{
		boost::mpl::at<v_expr, boost::mpl::int_<T::value> >::type::value(g_map,kmap,gs,spacing,cols,coeff);
	}

};

/*! \brief It model an expression expr1 + ... exprn
 *
 * \tparam expr.. two or more expression to be summed
 * \tparam Sys_eqs stystem specification
 *
 * ## Example
 *
 * \snippet FDScheme_unit_tests.hpp sum example
 *
 */
template<typename ... expr>
struct sum
{
	// Transform from variadic template to boost mpl vector
	typedef boost::mpl::vector<expr... > v_expr;

	// size of v_expr
	typedef typename boost::mpl::size<v_expr>::type v_sz;

	typedef typename boost::mpl::at<v_expr, boost::mpl::int_<v_sz::type::value - 1> >::type Sys_eqs;

	/*! \brief Calculate which colums of the Matrix are non zero
	 *
	 * \param pos position where the derivative is calculated
	 * \param gs Grid info
	 * \param cols non-zero colums calculated by the function
	 * \param coeff coefficent (constant in front of the derivative)
	 *
	 * ### Example
	 *
	 * \snippet FDScheme_unit_tests.hpp sum example
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap, const grid_sm<Sys_eqs::dims,void> & gs, typename Sys_eqs::stype (& spacing )[Sys_eqs::dims] , std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		// Sum functor
		sum_functor_value<v_expr> sm(g_map,kmap,gs,spacing,cols,coeff);

		// for each element in the expression calculate the non-zero Matrix elements
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_sz::type::value - 1> >(sm);
	}

	/*! \brief Calculate the position in the cell where the sum operator is performed
	 *
	 * it is done for the first element, the others are supposed to be in the same position
	 * it just return the position of the staggered property in the first expression
	 *
	 * \param position where we are calculating the derivative
	 * \param gs Grid info
	 * \param s_pos staggered position of the properties
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
	{
		return boost::mpl::at<v_expr, boost::mpl::int_<0> >::type::position(pos,gs,s_pos);
	}
};


template<typename arg, typename Sys_eqs>
struct minus
{
	/*! \brief Create the row of the Matrix
	 *
	 * \tparam ord
	 *
	 * \snipper FDScheme_unit_tests.hpp minus example
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap, const grid_sm<Sys_eqs::dims,void> & gs, typename Sys_eqs::stype (& spacing )[Sys_eqs::dims] , std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		arg::value(g_map,kmap,gs,spacing,cols,-coeff);
	}

	/*! \brief Calculate the position where the minus is calculated
	 *
	 * it just return the position of the staggered property at first expression
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}
};

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_ */
