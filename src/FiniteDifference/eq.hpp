/*
 * eq.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_HPP_

#define EQS_FIELD 0
#define EQS_POS 1

#define STAGGERED_GRID 1
#define NORMAL_GRID 0

#define PERIODIC true
#define NON_PERIODIC false

#include "data_type/scalar.hpp"
#include "util/util_num.hpp"

/*! \brief Equation
 *
 * It model an equation like expr1 = expr2
 *
 * \tparam expr1
 * \tparam expr2
 *
 */
template<typename expr1,typename expr2,typename Sys_eqs>
class Eq
{
	/*! \brief Create the row of the Matrix
	 *
	 * \tparam ord
	 *
	 */
	template<unsigned int ord=EQS_FIELD> static void value(const grid_key_dx<Sys_eqs::dims> & pos)
	{
		if (EQS_FIELD)
			value_f(pos);
		else
			value_s(pos);
	}

	/*! \brief fill the row
	 *
	 *
	 */
	static openfpm::vector<cval<typename Sys_eqs::stype>> value_s(grid_key_dx<Sys_eqs::dims> & it)
	{
		return expr1::value_s(it) - expr2::value_s(it);
	}

	/*! \brief fill the row
	 *
	 *
	 */
	static void value_f(grid_key_dx<Sys_eqs::dims> & it)
	{
		return expr1::value_s(it) - expr2::value_s(it);
	}
};


// spatial position + value

template<unsigned int dim,typename T>
struct pos_val
{
	/*! \brief Initialize to zero the value
	 *
	 */
	pos_val()
	{
		value = 0.0;
	}

	grid_key_dx<dim> pos;
	T value;
};

template<unsigned int f, typename Sys_eqs>
class Field
{
	typedef typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type map_grid;

public:

	/*! \brief fill the row
	 *
	 *
	 */
	static void value(const map_grid & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap, const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		if (Sys_eqs::ord == EQS_FIELD)
			cols[g_map.template get<0>(kmap)*Sys_eqs::nvar + f] += coeff;
		else
			cols[g_map.template get<0>(kmap) + f * gs.size()] += coeff;
	}
};

class ConstField
{

};

inline size_t mat_factor(size_t nvar, size_t sz, const size_t ord)
{
	if (ord == EQS_FIELD)
		return nvar;
	else
		return sz;
}

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_HPP_ */
