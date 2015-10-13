/*
 * sum.hpp
 *
 *  Created on: Oct 13, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_

#include <boost/mpl/vector.hpp>

/*! \brief It model an expression expr1 + ... exprn
 *
 * \tparam expr.. two or more expression to be summed
 * \tparam Sys_eqs stystem specification
 *
 */
template<typename ... expr>
class sum
{
	// Transform from variadic template to boost mpl vector

	typedef boost::mpl::vector<expr... > v_expr;


};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_ */
