/*
 * vector_dist_operators_functions.hpp
 *
 *  Created on: Jul 17, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_FUNCTIONS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_FUNCTIONS_HPP_

/*! A macro to define single value function specialization that apply the function component-wise
 *
 * \param fun function name
 * \param ID function ID
 *
 */
#define CREATE_VDIST_ARG_FUNC(fun_base,fun_name,OP_ID) \
\
\
template <typename exp1>\
class vector_dist_expression_op<exp1,void,OP_ID>\
{\
	const exp1 o1;\
\
public:\
\
	vector_dist_expression_op(const exp1 & o1)\
	:o1(o1)\
	{}\
\
	inline void init() const\
	{\
		o1.init();\
	}\
\
	template<typename r_type=typename std::remove_reference<decltype(fun_base(o1.value(vect_dist_key_dx(0))))>::type > inline r_type value(const vect_dist_key_dx & key) const\
	{\
		return fun_base(o1.value(key));\
	}\
};\
\
\
template<typename exp1, typename exp2_, unsigned int op1>\
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,void,OP_ID>\
fun_name(const vector_dist_expression_op<exp1,exp2_,op1> & va)\
{\
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,void,OP_ID> exp_sum(va);\
\
	return exp_sum;\
}\
\
\
template<unsigned int prp1, typename v1>\
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,void,OP_ID>\
fun_name(const vector_dist_expression<prp1,v1> & va)\
{\
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,void,OP_ID> exp_sum(va);\
\
	return exp_sum;\
}


CREATE_VDIST_ARG_FUNC(norm,norm,VECT_NORM)
CREATE_VDIST_ARG_FUNC(norm2,norm2,VECT_NORM2)
CREATE_VDIST_ARG_FUNC(abs,abs,POINT_ABS)
CREATE_VDIST_ARG_FUNC(exp,exp,POINT_EXP)
CREATE_VDIST_ARG_FUNC(exp2,exp2,POINT_EXP2)
CREATE_VDIST_ARG_FUNC(expm1,expm1,POINT_EXPM1)
CREATE_VDIST_ARG_FUNC(log,log,POINT_LOG)
CREATE_VDIST_ARG_FUNC(log10,log10,POINT_LOG10)
CREATE_VDIST_ARG_FUNC(log2,log2,POINT_LOG2)
CREATE_VDIST_ARG_FUNC(log1p,log1p,POINT_LOG1P)
CREATE_VDIST_ARG_FUNC(sqrt,sqrt,POINT_SQRT)
CREATE_VDIST_ARG_FUNC(cbrt,cbrt,POINT_CBRT)
CREATE_VDIST_ARG_FUNC(sin,sin,POINT_SIN)
CREATE_VDIST_ARG_FUNC(cos,cos,POINT_COS)
CREATE_VDIST_ARG_FUNC(tan,tan,POINT_TAN)
CREATE_VDIST_ARG_FUNC(asin,asin,POINT_ASIN)
CREATE_VDIST_ARG_FUNC(acos,acos,POINT_ACOS)
CREATE_VDIST_ARG_FUNC(atan,atan,POINT_ATAN)
CREATE_VDIST_ARG_FUNC(sinh,sinh,POINT_SINH)
CREATE_VDIST_ARG_FUNC(cosh,cosh,POINT_COSH)
CREATE_VDIST_ARG_FUNC(tanh,tanh,POINT_TANH)
CREATE_VDIST_ARG_FUNC(asinh,asinh,POINT_ASINH)
CREATE_VDIST_ARG_FUNC(acosh,acosh,POINT_ACOSH)
CREATE_VDIST_ARG_FUNC(atanh,atanh,POINT_ATANH)
CREATE_VDIST_ARG_FUNC(erf,erf,POINT_ERF)
CREATE_VDIST_ARG_FUNC(erfc,erfc,POINT_ERFC)
CREATE_VDIST_ARG_FUNC(tgamma,tgamma,POINT_TGAMMA)
CREATE_VDIST_ARG_FUNC(lgamma,lgamma,POINT_LGAMMA)
CREATE_VDIST_ARG_FUNC(ceil,ceil,POINT_CEIL)
CREATE_VDIST_ARG_FUNC(floor,floor,POINT_FLOOR)
CREATE_VDIST_ARG_FUNC(trunc,trunc,POINT_TRUNC)
CREATE_VDIST_ARG_FUNC(round,round,POINT_ROUND)
CREATE_VDIST_ARG_FUNC(nearbyint,nearbyint,POINT_NEARBYINT)
CREATE_VDIST_ARG_FUNC(rint,rint,POINT_RINT)


/*! A macro to define single value function specialization that apply the function component-wise
 *
 * \param fun function name
 * \param ID function ID
 *
 */
#define CREATE_VDIST_ARG2_FUNC(fun_base,fun_name,OP_ID) \
\
\
template <typename exp1,typename exp2>\
class vector_dist_expression_op<exp1,exp2,OP_ID>\
{\
	const exp1 o1;\
	const exp2 o2;\
\
public:\
\
	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)\
	:o1(o1),o2(o2)\
	{}\
\
	inline void init() const\
	{\
		o1.init();\
		o2.init();\
	}\
\
	template<typename r_type=typename std::remove_reference<decltype(fun_base(o1.value(vect_dist_key_dx(0)),o2.value(vect_dist_key_dx(0)) ))>::type > inline r_type value(const vect_dist_key_dx & key) const\
	{\
		return fun_base(o1.value(key),o2.value(key));\
	}\
};\
\
\
template<unsigned int p1, unsigned int p2, typename v1, typename v2>\
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,OP_ID>\
fun_name(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)\
{\
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename exp1 , typename exp2, unsigned int op1, unsigned int prp1, typename v1>\
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,OP_ID>\
fun_name(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<prp1,v1> & vb)\
{\
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename exp1 , typename exp2, unsigned int op1, unsigned int prp1, typename v1>\
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,OP_ID>\
fun_name(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)\
{\
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename exp1 , typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>\
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,OP_ID>\
fun_name(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)\
{\
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<unsigned int prp1 , typename v1>\
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,OP_ID>\
fun_name(const vector_dist_expression<prp1,v1> & va, double d)\
{\
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,OP_ID> exp_sum(va,vector_dist_expression<0,double>(d));\
\
	return exp_sum;\
}\
\
template<unsigned int prp1 , typename v1>\
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,OP_ID>\
fun_name(double d, const vector_dist_expression<prp1,v1> & vb)\
{\
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,OP_ID> exp_sum(vector_dist_expression<0,double>(d),vb);\
\
	return exp_sum;\
}\
\
template<typename exp1 , typename exp2, unsigned int op1>\
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,OP_ID>\
fun_name(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)\
{\
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,OP_ID> exp_sum(va,vector_dist_expression<0,double>(d));\
\
	return exp_sum;\
}\
\
template<typename exp1 , typename exp2, unsigned int op1>\
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression_op<exp1,exp2,op1>,OP_ID>\
fun_name(double d, const vector_dist_expression_op<exp1,exp2,op1> & va)\
{\
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,OP_ID> exp_sum(vector_dist_expression<0,double>(d),va);\
\
	return exp_sum;\
}

CREATE_VDIST_ARG2_FUNC(pmul,pmul,VECT_PMUL)

////////// Special function reduce /////////////////////////


/*! \brief expression that encapsulate a vector reduction expression
 *
 * \tparam exp1 expression 1
 * \tparam vector_type type of vector on which the expression is acting
 *
 */
template <typename exp1, typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_SUM_REDUCE>
{

	//! expression on which apply the reduction
	const exp1 o1;

	//! return type of this expression
	typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx()))>::rtype rtype;

	//! return type of the calculated value (without reference)
	mutable typename std::remove_reference<rtype>::type val;

	//! vector on which we apply the reduction expression
	const vector_type & vd;

public:

	//! constructor from an epxression exp1 and a vector vd
	vector_dist_expression_op(const exp1 & o1, const vector_type & vd)
	:o1(o1),val(0),vd(vd)
	{}

	//! sum reduction require initialization where we calculate the reduction
	// this produce a cache for the calculated value
	inline void init() const
	{
		o1.init();

		val = 0.0;

		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			val += o1.value(key);

			++it;
		}
	}

	//! it return the result of the expression
	inline typename std::remove_reference<rtype>::type get()
	{
		init();
		return value(vect_dist_key_dx());
	}

	//! it return the result of the expression (precalculated before)
	template<typename r_type= typename std::remove_reference<rtype>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return val;
	}
};

//! Reduce function (it generate an expression)
template<typename exp1, typename exp2_, unsigned int op1, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,vector_type,VECT_SUM_REDUCE>
rsum(const vector_dist_expression_op<exp1,exp2_,op1> & va, const vector_type & vd)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,vector_type,VECT_SUM_REDUCE> exp_sum(va,vd);

	return exp_sum;
}

//! Reduce function (It generate an expression)
template<unsigned int prp1, typename v1, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_type,VECT_SUM_REDUCE>
rsum(const vector_dist_expression<prp1,v1> & va, const vector_type & vd)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_type,VECT_SUM_REDUCE> exp_sum(va,vd);

	return exp_sum;
}

namespace openfpm
{
	/*! \brief General distance formula
	 *
	 *
	 */
	template <typename T, typename P> auto distance(T exp1, P exp2) -> decltype(norm(exp1 - exp2))
	{
		return norm(exp1 - exp2);
	}
}


#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_FUNCTIONS_HPP_ */
