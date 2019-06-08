/*
 * vector_dist_operators_cuda.cuh
 *
 *  Created on: May 31, 2019
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_OPERATORS_CUDA_CUH_
#define VECTOR_DIST_OPERATORS_CUDA_CUH_

constexpr unsigned int PROP_POS =(unsigned int)-1;

/*! \brief selector for position or properties left side expression
 *
 * \tparam vector type of the original vector
 *
 * \tparam prp property id
 *
 */
template <typename vector, unsigned int prp>
struct pos_or_propL
{
	//! return the value (position or property) of the particle k in the vector v
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}
};

/*! \brief selector for position or properties left side expression
 *
 * \tparam vector type of the original vector
 *
 * \tparam prp property id
 *
 */
template <typename vector, unsigned int prp>
struct pos_or_propL_ker
{
	//! return the value (position or property) of the particle k in the vector v
	__device__ static inline auto value(vector & v, const unsigned int & k) -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}
};


/*! \brief selector for position or properties left side
 *
 * \tparam vector type of the original vector
 *
 * \tparam prp property id
 *
 */
template <typename vector>
struct pos_or_propL<vector,PROP_POS>
{
#ifdef SE_CLASS3

	//! return the value (position or property) of the particle k in the vector v
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(getExprL(v.getPos(k).getReference()))
	{
		return getExprL(v.getPos(k).getReference());
	}

#else

	//! return the value (position or property) of the particle k in the vector v
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(ger<vector::dims,typename vector::stype>::getExprL(v.getPos(k)))
	{
		return ger<vector::dims,typename vector::stype>::getExprL(v.getPos(k));
	}

#endif
};

/*! \brief selector for position or properties left side
 *
 * \tparam vector type of the original vector
 *
 * \tparam prp property id
 *
 */
template <typename vector>
struct pos_or_propL_ker<vector,PROP_POS>
{
#ifdef SE_CLASS3

	//! return the value (position or property) of the particle k in the vector v
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(getExprL(v.getPos(k).getReference()))
	{
		return getExprL(v.getPos(k).getReference());
	}

#else

	//! return the value (position or property) of the particle k in the vector v
	__device__ static inline auto value(vector & v, const unsigned int & k) -> decltype(ger<vector::dims,typename vector::stype>::getExprL(v.getPos(k)))
	{
		return ger<vector::dims,typename vector::stype>::getExprL(v.getPos(k));
	}

#endif
};

/*! \brief selector for position or properties right side position
 *
 * \tparam vector type of the original vector
 *
 * \tparam prp property id
 *
 */
template <typename vector, unsigned int prp>
struct pos_or_propR
{
	//! return the value (position or property) of the particle k in the vector v
	__device__ __host__ static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}

	//! return the value (position or property) of the particle k in the vector v
	__device__ __host__ static inline auto value(vector & v, const unsigned int & k) -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}
};


/*! \brief selector for position or properties right side
 *
 * \tparam vector type of the original vector
 *
 * \tparam prp property id
 *
 */
template <typename vector>
struct pos_or_propR<vector,PROP_POS>
{
	//! return the value (position or property) of the particle k in the vector v
	__device__ __host__ static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(ger<vector::dims,typename vector::stype>::getExprR(v.getPos(k)))
	{
		return ger<vector::dims,typename vector::stype>::getExprR(v.getPos(k));
	}

	//! return the value (position or property) of the particle k in the vector v
	__device__ __host__ static inline auto value(vector & v, const unsigned int & k) -> decltype(ger<vector::dims,typename vector::stype>::getExprR(v.getPos(k)))
	{
		return ger<vector::dims,typename vector::stype>::getExprR(v.getPos(k));
	}
};

template<unsigned int prp ,int impl>
struct vector_dist_op_compute_op
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{}

	template<typename vector>
	static void compute_const(vector & v,double d)
	{}
};

template<unsigned int prp>
struct vector_dist_op_compute_op<prp,comp_host>
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			pos_or_propL<vector,prp>::value(v,key) = v_exp.value(key);

			++it;
		}
	}

	template<typename vector>
	static void compute_const(vector & v,double d)
	{
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			pos_or_propL<vector,prp>::value(v,key) = d;

			++it;
		}
	}
};

#ifdef __NVCC__

template<unsigned int prp, unsigned int dim ,typename vector, typename expr>
__global__ void compute_expr_ker_vv(vector vd, expr v_exp)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size())	{return;}

	for (unsigned int i = 0 ; i < dim ; i++)
	{
		vd.template get<prp>(p)[i] = v_exp.value(p).get(i);
	}
}

template<unsigned int prp, typename vector, typename expr>
__global__ void compute_expr_ker_v(vector vd, expr v_exp)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size())	{return;}

	vd.template get<prp>(p) = v_exp.value(p);
}

template<unsigned int prp, typename vector, typename expr>
__global__ void compute_expr_ker(vector vd, expr v_exp)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size_local())	{return;}

	pos_or_propL_ker<vector,prp>::value(vd,p) = v_exp.value(p);
}

template<unsigned int prp, typename vector>
__global__ void compute_double_ker(vector vd, double d)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size_local())	{return;}

	pos_or_propL_ker<vector,prp>::value(vd,p) = d;
}

template<unsigned int prp>
struct vector_dist_op_compute_op<prp,comp_dev>
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getDomainIteratorGPU(256);

		compute_expr_ker<prp><<<ite.wthr,ite.thr>>>(v,v_exp);
	}

	template<typename vector>
	static void compute_const(vector & v,double d)
	{
		auto ite  = v.getDomainIteratorGPU(256);

		compute_double_ker<prp><<<ite.wthr,ite.thr>>>(v,d);
	}
};

#endif


#endif /* VECTOR_DIST_OPERATORS_CUDA_CUH_ */
