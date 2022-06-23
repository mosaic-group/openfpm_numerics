/*
 * vector_dist_operators_cuda.cuh
 *
 *  Created on: May 31, 2019
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_OPERATORS_CUDA_CUH_
#define VECTOR_DIST_OPERATORS_CUDA_CUH_

#include "Space/Shape/Point.hpp"
#include "util/cuda_launch.hpp"
#include <utility>

#ifdef SE_CLASS1
template<bool is_subset>
struct SubsetSelector_impl{
    template<typename particle_type,typename subset_type>
    static void check(particle_type &particles,subset_type &particle_subset)
    {
    }
};

template<>
struct SubsetSelector_impl<true>
{
    template<typename particle_type,typename subset_type>
    static void check(particle_type &particles,subset_type &particle_subset){
        //This getMapCtr needs to be created or fixed for cuda!
       /* if(particles.getMapCtr()!=particle_subset.getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a subset update after map."<<std::endl;
        }*/
    }
};
#endif

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
	typedef typename boost::mpl::at<typename vector::value_type::type, boost::mpl::int_<prp>>::type property_act;

	//! return the value (position or property) of the particle k in the vector v
	__device__ __host__ static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}

	//! return the value (position or property) of the particle k in the vector v
	static inline auto value_type(vector && v, const vect_dist_key_dx & k) -> decltype(v.template getProp<prp>(k))
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
	typedef typename Point<vector::dims,typename vector::stype>::type_native property_act;

#ifdef SE_CLASS3

	//! return the value (position or property) of the particle k in the vector v
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(getExprL(v.getPos(k).getReference()))
	{
		return getExprL(v.getPos(k).getReference());
	}

#else

	//! return the value (position or property) of the particle k in the vector v
	__device__ __host__ static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(ger<vector::dims,typename vector::stype>::getExprL(v.getPos(k)))
	{
		return ger<vector::dims,typename vector::stype>::getExprL(v.getPos(k));
	}

	//! return the value (position or property) of the particle k in the vector v
	static inline auto value_type(vector && v, const vect_dist_key_dx & k) -> decltype(v.getPos(k))
	{
		return v.getPos(k);
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

template<unsigned int prp ,bool is_sort, int impl>
struct vector_dist_op_compute_op
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{}

	template<unsigned int n, typename vector, typename expr>
	static void compute_expr_slice(vector & v,expr & v_exp, int (& comp)[n])
	{}

	template<typename vector>
	static void compute_const(vector & v,double d)
	{}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<unsigned int, bool is_valid>
struct get_vector_dist_expression_op
{
	template<typename exp_type>
	__device__ __host__ inline static auto get(exp_type & o1, const vect_dist_key_dx & key) -> decltype(o1.value(vect_dist_key_dx(0)))
	{
		return o1.value(key);
	}

	template<unsigned int prop, typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig)
	{
		pos_or_propL<vector_type,exp_type::prop>::value(v,key) = o1.value(key_orig);
	}

	template<unsigned int prop, typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const unsigned int & key_orig)
	{
		pos_or_propL<vector_type,exp_type::prop>::value(v,key) = o1.value(key_orig);
	}

	template<unsigned int prop, typename vector_type>
	inline static void assign_double(double d, vector_type & v, const vect_dist_key_dx & key)
	{
		pos_or_propL<vector_type,prop>::value(v,key) = d;
	}
};

template<>
struct get_vector_dist_expression_op<1,false>
{
	template<typename exp_type>
	__device__ __host__ static int get(exp_type & o1, const vect_dist_key_dx & key, const int (& comp)[1])
	{
		printf("ERROR: Slicer, the expression is incorrect, please check it\n");
		return 0;
	}

	template<unsigned int prop, typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig, const int (& comp)[1])
	{
		printf("ERROR: Slicer, the expression is incorrect, please check it\n");
	}

	template<unsigned int prop, typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const vect_dist_key_dx & key_orig, const int (& comp)[1])
	{
		printf("ERROR: Slicer, the expression is incorrect, please check it\n");
	}

        template<unsigned int prop,typename exp_type, typename vector_type>
        __device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig, const Point<1,int> & comp)
        {
                printf("ERROR: Slicer, the expression is incorrect, please check it\n");
        }

        template<unsigned int prop,typename exp_type, typename vector_type>
        __device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const unsigned int & key_orig, const Point<1,int> & comp)
        {
                printf("ERROR: Slicer, the expression is incorrect, please check it\n");
        }

	template<unsigned int prop, typename vector_type>
	inline static void assign_double(double d, vector_type & v, const vect_dist_key_dx & key, const int (& comp)[1])
	{
		printf("ERROR: Slicer, the expression is incorrect, please check it\n");
	}
};

template<>
struct get_vector_dist_expression_op<1,true>
{
	template<typename exp_type>
	__device__ __host__ static auto get(exp_type & o1, const vect_dist_key_dx & key, const int (& comp)[1]) -> decltype(o1.value(vect_dist_key_dx(0))[0])
	{
		return o1.value(key)[comp[0]];
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig, const int (& comp)[1])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]] = o1.value(key_orig);
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const unsigned int & key_orig, const int (& comp)[1])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]] = o1.value(key_orig);
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig, const Point<1,int> & comp)
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]] = o1.value(key_orig);
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const unsigned int & key_orig, const Point<1,int> & comp)
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]] = o1.value(key_orig);
	}

	template<unsigned int prop, typename vector_type>
	inline static void assign_double(double d, vector_type & v, const vect_dist_key_dx & key, const int (& comp)[1])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]] = d;
	}
};

template<>
struct get_vector_dist_expression_op<2,true>
{
	template<typename exp_type>
	__device__ __host__ static auto get(exp_type & o1, const vect_dist_key_dx & key, const int (& comp)[2]) -> decltype(o1.value(vect_dist_key_dx(0))[0][0])
	{
		return o1.value(key)[comp[0]][comp[1]];
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig, const int (& comp)[2])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = o1.value(key_orig);
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const vect_dist_key_dx & key_orig, const Point<2,int> & comp)
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = o1.value(key_orig);
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const unsigned int & key_orig, const int (& comp)[2])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = o1.value(key_orig);
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	__device__ __host__ inline static void assign(exp_type & o1, vector_type & v, const unsigned int & key, const unsigned int & key_orig, const Point<2,int> & comp)
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = o1.value(key_orig);
	}

	template<unsigned int prop, typename vector_type>
	inline static void assign_double(double d, vector_type & v, const vect_dist_key_dx & key, const int (& comp)[2])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = d;
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<unsigned int prp>
struct vector_dist_op_compute_op<prp,false,comp_host>
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{
		v_exp.init();

        auto it = v.getDomainIterator();

        while (it.isNext())
        {
                auto key = it.get();
                auto key_orig = v.getOriginKey(key);

                pos_or_propL<vector,prp>::value(v,key) = v_exp.value(key_orig);

                ++it;
		}
	}

	template<unsigned int n, typename vector, typename expr>
	static void compute_expr_slice(vector & v,expr & v_exp, int (& comp)[n])
	{
		typedef typename std::remove_const<typename std::remove_reference<decltype(pos_or_propL<vector,prp>::value_type(std::declval<vector>(),vect_dist_key_dx(0)))>::type>::type property_act;

		v_exp.init();

#ifdef SE_CLASS1
        auto &v2=v_exp.getVector();

        SubsetSelector_impl<std::remove_reference<decltype(v)>::type::is_it_a_subset::value>::check(v2,v);
#endif

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();
			auto key_orig = v.getOriginKey(key);

			get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<prp>(v_exp,v,key,key_orig,comp);

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
                auto key_orig = v.getOriginKey(key);

                pos_or_propL<vector,prp>::value(v,key) = d;

                ++it;

		}
	}
};


#ifdef __NVCC__

template<unsigned int prp, unsigned int dim ,typename vector, typename expr>
__global__ void compute_expr_ker_vv(vector vd, expr v_exp)
{
	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size())	{return;}

	for (unsigned int i = 0 ; i < dim ; i++)
	{
		vd.template get<prp>(p)[i] = v_exp.value(p).get(i);
	}
}

template<unsigned int prp, typename vector, typename expr>
__global__ void compute_expr_ker_v(vector vd, expr v_exp)
{
	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size())	{return;}

	vd.template get<prp>(p) = v_exp.value(p);
}

template<unsigned int prp, typename vector, typename expr>
__global__ void compute_expr_ker(vector vd, expr v_exp)
{
	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size_local())	{return;}

	pos_or_propL_ker<vector,prp>::value(vd,p) = v_exp.value(p);
}

namespace openfpm
{

  template<typename _Tp, typename _Up = _Tp&&>
    __device__ __host__ _Up
    __declval(int);

  template<typename _Tp>
    __device__ __host__ _Tp
    __declval(long);

  template<typename _Tp>
    __device__ __host__ auto declval() noexcept -> decltype(__declval<_Tp>(0))
    {
      return __declval<_Tp>(0);
    }
}

template<unsigned int prp, unsigned int n, typename vector, typename expr>
__global__ void compute_expr_ker_slice(vector vd, expr v_exp, Point<n,int> comp)
{
	typedef typename std::remove_const<typename std::remove_reference<decltype(pos_or_propL<vector,prp>::value_type(openfpm::declval<vector>(),vect_dist_key_dx(0)))>::type>::type property_act;

	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size_local())	{return;}

	get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<prp>(v_exp,vd,p,p,comp);
}

template<unsigned int prp, typename vector>
__global__ void compute_double_ker(vector vd, double d)
{
	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size_local())	{return;}

	pos_or_propL_ker<vector,prp>::value(vd,p) = d;
}

/////////// SORTED VERSION //

template<unsigned int prp, unsigned int dim ,typename vector, typename NN_type, typename expr>
__global__ void compute_expr_ker_sort_vv(vector vd, NN_type NN, expr v_exp)
{
	unsigned int p;

    GET_PARTICLE_SORT(p,NN);

	for (unsigned int i = 0 ; i < dim ; i++)
	{
		vd.template get<prp>(p)[i] = v_exp.value(p).get(i);
	}
}

template<unsigned int prp, typename vector, typename NN_type, typename expr>
__global__ void compute_expr_ker_sort_v(vector vd, NN_type NN, expr v_exp)
{
	unsigned int p;

    GET_PARTICLE_SORT(p,NN);

	vd.template get<prp>(p) = v_exp.value(p);
}

template<unsigned int prp, typename vector, typename expr, typename NN_type>
__global__ void compute_expr_ker_sort(vector vd, NN_type NN, expr v_exp)
{
	unsigned int p;

    GET_PARTICLE_SORT(p,NN);

	pos_or_propL_ker<vector,prp>::value(vd,p) = v_exp.value(p);
}


/////////////////////////////

template<unsigned int prp>
struct vector_dist_op_compute_op<prp,false,comp_dev>
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getDomainIteratorGPU(256);

		CUDA_LAUNCH((compute_expr_ker<prp>),ite,v,v_exp);
	}

	template<unsigned int n, typename vector, typename expr>
	static void compute_expr_slice(vector & v,expr & v_exp, int (& comp)[n])
	{
		v_exp.init();

		auto ite  = v.getDomainIteratorGPU(256);

		Point<n,int> comp_;
		for (int i = 0 ; i < n ; i++)	{comp_[i] = comp[i];}

		CUDA_LAUNCH((compute_expr_ker_slice<prp,n>),ite,v,v_exp,comp);
	}

	template<typename vector, typename expr>
	static void compute_expr_v(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getGPUIterator(256);

		CUDA_LAUNCH((compute_expr_ker_v<prp>),ite,v,v_exp);
	}

	template<unsigned int dim, typename vector, typename expr>
	static void compute_expr_vv(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getGPUIterator(256);

		CUDA_LAUNCH((compute_expr_ker_vv<prp,dim>),ite,v,v_exp);
	}

	template<typename vector>
	static void compute_const(vector & v,double d)
	{
		auto ite  = v.getDomainIteratorGPU(256);

		CUDA_LAUNCH((compute_double_ker<prp>),ite,v,d);
	}
};

template<unsigned int prp>
struct vector_dist_op_compute_op<prp,true,comp_dev>
{
	template<typename vector, typename expr>
	static void compute_expr(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getDomainIteratorGPU(256);

		auto NN = v_exp.getNN();

		CUDA_LAUNCH((compute_expr_ker_sort<prp>),ite,v,*NN,v_exp);
	}

	template<typename vector, typename expr>
	static void compute_expr_v(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getGPUIterator(256);

		auto NN = v_exp.getNN();

		CUDA_LAUNCH((compute_expr_ker_sort_v<prp>),ite,v,*NN,v_exp);
	}

	template<unsigned int dim, typename vector, typename expr>
	static void compute_expr_vv(vector & v,expr & v_exp)
	{
		v_exp.init();

		auto ite  = v.getGPUIterator(256);

		auto NN = v_exp.getNN();

		CUDA_LAUNCH((compute_expr_ker_sort_vv<prp,dim>),ite,v,*NN,v_exp);
	}
};

#endif


#endif /* VECTOR_DIST_OPERATORS_CUDA_CUH_ */
