/*
 * vector_dist_operators_extensions.hpp
 *
 *  Created on: Jul 18, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_EXTENSIONS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_EXTENSIONS_HPP_


/*! \Create an expression from a Point
 *
 * \tpatam prp property
 * \param v
 *
 */
template <unsigned int dim, typename T> inline vector_dist_expression<16384,Point<dim,T> > getVExpr(Point<dim,T> & v)
{
	vector_dist_expression<(unsigned int)16384,Point<dim,T>> exp_v(v);

	return exp_v;
}


/*! \brief Main class that encapsulate a vector properties
 *
 * \tparam prp property involved
 * \tparam vector involved
 *
 */
template<typename point>
class vector_dist_expression<16384,point>
{
	point p;

public:

	vector_dist_expression(point p)
	:p(p)
	{}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 */
	inline point value(const vect_dist_key_dx & k) const
	{
		return p;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_EXTENSIONS_HPP_ */
