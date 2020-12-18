/*
 * vector_dist_operators_functions.hpp
 *
 *  Created on: Jul 17, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_FUNCTIONS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_FUNCTIONS_HPP_

#ifdef __NVCC__
#include "cuda/vector_dist_operators_cuda.cuh"
#endif

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
	typedef typename exp1::is_ker is_ker;\
\
	typedef typename vector_result<typename exp1::vtype,void>::type vtype;\
\
	typedef typename vector_is_sort_result<exp1::is_sort::value,false>::type is_sort;\
\
	typedef typename nn_type_result<typename exp1::NN_type,void>::type NN_type;\
\
	vector_dist_expression_op(const exp1 & o1)\
	:o1(o1)\
	{}\
\
	inline NN_type * getNN()\
	{\
		return nn_type_result<typename exp1::NN_type,void>::getNN(o1);\
	}\
\
	const vtype & getVector()\
	{\
		return vector_result<typename exp1::vtype,void>::getVector(o1);\
	}\
\
	inline const exp1 & getExpr() const\
	{\
		return o1;\
	}\
	\
	inline void init() const\
	{\
		o1.init();\
	}\
\
	template<typename r_type=typename std::remove_reference<decltype(fun_base(o1.value(vect_dist_key_dx(0))))>::type > \
	inline r_type value(const vect_dist_key_dx & key) const\
	{\
		return fun_base(o1.value(key));\
	}\
\
	template<typename r_type=typename std::remove_reference<decltype(fun_base(o1.value(vect_dist_key_dx(0))))>::type > \
	__device__ __host__ inline r_type value(const unsigned int & key) const\
	{\
		return fun_base(o1.value(key));\
	}\
\
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
	typedef std::integral_constant<bool,exp1::is_ker::value || exp1::is_ker::value> is_ker;\
\
	typedef typename vector_result<typename exp1::vtype,typename exp2::vtype>::type vtype;\
\
	typedef typename vector_is_sort_result<exp1::is_sort::value,exp2::is_sort::value>::type is_sort;\
\
	typedef typename nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::type NN_type;\
\
	inline NN_type * getNN()\
	{\
		return nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::getNN(o1,o2);\
	}\
\
	vtype & getVector()\
	{\
		return vector_result<typename exp1::vtype,typename exp2::vtype>::getVector(o1,o2);\
	}\
\
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

template<typename val_type, bool is_sort, bool is_scalar = is_Point<val_type>::type::value>
struct point_scalar_process
{
	typedef aggregate<val_type> type;

	template<typename vector_type, typename expression>
	static void process(val_type & val, vector_type & ve, expression & o1)
	{
#ifdef __NVCC__

		auto vek = ve.toKernel();
		vector_dist_op_compute_op<0,is_sort,comp_dev>::compute_expr_v(vek,o1);

		exp_tmp2[0].resize(sizeof(val_type));

		auto & v_cl = create_vcluster<CudaMemory>();

		mgpu::reduce((val_type *)ve.template getDeviceBuffer<0>(), ve.size(), (val_type *)(exp_tmp2[0].getDevicePointer()), mgpu::plus_t<val_type>(), v_cl.getmgpuContext());

		exp_tmp2[0].deviceToHost();

		val = *(val_type *)(exp_tmp2[0].getPointer());
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error: to make expression work on GPU the file must be compiled on GPU" << std::endl;
#endif
	}
};

template<typename val_type, bool is_sort>
struct point_scalar_process<val_type,is_sort,true>
{
	typedef val_type type;

	template<typename vector_type, typename expression>
	static void process(val_type & val, vector_type & ve, expression & o1)
	{
#ifdef __NVCC__

//		auto ite = ve.getGPUIterator(256);

//		compute_expr_ker_vv<0,val_type::dims><<<ite.wthr,ite.thr>>>(ve.toKernel(),o1);

		auto vek = ve.toKernel();
		vector_dist_op_compute_op<0,is_sort,comp_dev>::template compute_expr_vv<val_type::dims>(vek,o1);

		exp_tmp2[0].resize(sizeof(val_type));

		size_t offset = 0;

		auto & v_cl = create_vcluster<CudaMemory>();

		for (size_t i = 0 ; i < val_type::dims ; i++)
		{
			mgpu::reduce(&((typename val_type::coord_type *)ve.template getDeviceBuffer<0>())[offset],
						 ve.size(),
						 (typename val_type::coord_type *)(exp_tmp2[0].getDevicePointer()),
						 mgpu::plus_t<typename val_type::coord_type>(),
						 v_cl.getmgpuContext());

			exp_tmp2[0].deviceToHost();

			val.get(i) = *(typename val_type::coord_type *)exp_tmp2[0].getPointer();

			offset += ve.capacity();
		}

#else
		std::cout << __FILE__ << ":" << __LINE__ << " error: to make expression work on GPU the file must be compiled on GPU" << std::endl;
#endif
	}
};

/*! \brief expression that encapsulate a vector reduction expression
 *
 * \tparam exp1 expression 1
 * \tparam vector_type type of vector on which the expression is acting
 *
 */
template <typename exp1>
class vector_dist_expression_op<exp1,void,VECT_SUM_REDUCE>
{

	//! expression on which apply the reduction
	const exp1 o1;

	//! return type of this expression
	typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx()))>::rtype rtype;

	//! return type of the calculated value (without reference)
	mutable typename std::remove_reference<rtype>::type val;

public:

	//! Indicate if it is an in kernel expression
	typedef typename exp1::is_ker is_ker;

	//! return the vector type on which this expression operate
	typedef typename vector_result<typename exp1::vtype,void>::type vtype;

	//! result for is sort
	typedef typename vector_is_sort_result<exp1::is_sort::value,false>::type is_sort;

	//! NN_type
	typedef typename nn_type_result<typename exp1::NN_type,void>::type NN_type;

	//! constructor from an epxression exp1 and a vector vd
	vector_dist_expression_op(const exp1 & o1)
	:o1(o1),val(0)
	{}

	//! sum reduction require initialization where we calculate the reduction
	// this produce a cache for the calculated value
	inline void init() const
	{
		if (exp1::is_ker::value == true)
		{

#ifdef __NVCC__
			typedef decltype(val) val_type;

			// we have to do it on GPU

			openfpm::vector<typename point_scalar_process<val_type,is_sort::value>::type,
							CudaMemory,
							typename memory_traits_inte<typename point_scalar_process<val_type,is_sort::value>::type>::type,
							memory_traits_inte,
							openfpm::grow_policy_identity> ve;

			auto & orig_v = o1.getVector();

			if (exp_tmp.ref() == 0)
			{exp_tmp.incRef();}

			ve.setMemory(exp_tmp);
			ve.resize(orig_v.size_local());

			point_scalar_process<val_type,is_sort::value>::process(val,ve,o1);
#else
			std::cout << __FILE__ << ":" << __LINE__ << " error, to use expression on GPU you must compile with nvcc compiler " << std::endl;
#endif
		}
		else
		{
			const auto & orig_v = o1.getVector();

			o1.init();

			val = 0.0;

			auto it = orig_v.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				val += o1.value(key);

				++it;
			}
		}
	}

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_type * getNN() const
	{
		return nn_type_result<typename exp1::NN_type,void>::getNN(o1);
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
template<typename exp1, typename exp2_, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,void,VECT_SUM_REDUCE>
rsum(const vector_dist_expression_op<exp1,exp2_,op1> & va)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,void,VECT_SUM_REDUCE> exp_sum(va);

	return exp_sum;
}

//! Reduce function (It generate an expression)
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,void,VECT_SUM_REDUCE>
rsum(const vector_dist_expression<prp1,v1> & va)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,void,VECT_SUM_REDUCE> exp_sum(va);

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
