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

	const grid_sm<last::dims,void> & gs;

	//! position
	grid_key_dx<last::dims> & key;

	//! coefficent
	typename last::stype coeff;

	/*! \brief constructor
	 *
	 */
	sum_functor_value(grid_key_dx<last::dims> & key, const grid_sm<last::dims,void> & gs, typename last::stype coeff)
	:cols(cols),gs(gs),key(key),coeff(coeff)
	{};



	//! It call this function for every expression in the sum
	template<typename T>
	void operator()(T& t) const
	{
		boost::mpl::at<v_expr, boost::mpl::int_<T::value> >::type::value(key,gs,cols,coeff);
	}
};

/*! \brief It model an expression expr1 + ... exprn
 *
 * \tparam expr.. two or more expression to be summed
 * \tparam Sys_eqs stystem specification
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

	/*! \brief Create the row of the Matrix
	 *
	 * \tparam ord
	 *
	 */
	inline static void value(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		// Sum functor
		sum_functor_value<v_expr> sm(pos,gs,coeff);

		// for each element in the expression calculate the non-zero Matrix elements
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_sz::type::value - 1> >(sm);
	}

	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case on non staggered case this function just return pos, in case of staggered,
	 *  it calculate where the operator is calculated on a staggered grid
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_ */
