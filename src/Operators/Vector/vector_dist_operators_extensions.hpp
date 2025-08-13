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
template <unsigned int dim, typename T>
inline vector_dist_expression<16384,Point<dim,T> > getVExpr(Point<dim,T> & v)
{
	vector_dist_expression<(unsigned int)16384,Point<dim,T>> exp_v(v);

	return exp_v;
}


/*! \brief This class represent a constant parameter in a vector expression
 *
 * \tparam point type of point
 *
 */
template<typename point>
class vector_dist_expression<16384,point>
{
	//! constant point stored
	point p;

public:

	typedef void vtype;

	//! NN_type
	typedef void NN_type;

	//! vector expression from a constant point
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
	 * \param k where to evaluate the expression (ignored)
	 *
	 * \return the point stored
	 *
	 */
	__device__ __host__ inline point value(const vect_dist_key_dx & k) const
	{
		return p;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_EXTENSIONS_HPP_ */
