/*
 * vector_dist_operators.hpp
 *
 *  Created on: Jun 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_

#include "Vector/vector_dist.hpp"
#include "lib/pdata.hpp"
#include "cuda/vector_dist_operators_cuda.cuh"

#define PROP_CUSTOM (unsigned int)-2

#define VECT_SUM 1
#define VECT_SUB 2
#define VECT_MUL 3
#define VECT_DIV 4

#define VECT_APPLYKER_IN 7
#define VECT_APPLYKER_OUT 8
#define VECT_APPLYKER_REDUCE 9
#define VECT_APPLYKER_IN_GEN 10
#define VECT_APPLYKER_OUT_GEN 11
#define VECT_APPLYKER_REDUCE_GEN 12
#define VECT_APPLYKER_IN_SIM 13
#define VECT_APPLYKER_OUT_SIM 14
#define VECT_APPLYKER_REDUCE_SIM 15

#define VECT_APPLYKER_IN_GEN_SORT 16
#define VECT_APPLYKER_IN_SORT 17
#define VECT_APPLYKER_IN_SIM_SORT 18

#define VECT_NORM 56
#define VECT_NORM2 57
#define VECT_ABS 58
#define VECT_EXP 59
#define VECT_EXP2 60
#define VECT_EXPM1 61
#define VECT_LOG 62
#define VECT_LOG10 63
#define VECT_LOG2 64
#define VECT_LOG1P 65
#define VECT_SQRT 67
#define VECT_CBRT 68
#define VECT_SIN 69
#define VECT_COS 70
#define VECT_TAN 71
#define VECT_ASIN 72
#define VECT_ACOS 73
#define VECT_ATAN 74
#define VECT_SINH 75
#define VECT_COSH 76
#define VECT_TANH 77
#define VECT_ASINH 78
#define VECT_ACOSH 79
#define VECT_ATANH 80
#define VECT_ERF 81
#define VECT_ERFC 82
#define VECT_TGAMMA 83
#define VECT_LGAMMA 84
#define VECT_CEIL 85
#define VECT_FLOOR 86
#define VECT_TRUNC 87
#define VECT_ROUND 88
#define VECT_NEARBYINT 89
#define VECT_RINT 90
#define VECT_PMUL 91
#define VECT_SUB_UNI 92
#define VECT_SUM_REDUCE 93
#define VECT_COMP 94
#define VECT_NORM_INF 95


#define VECT_DCPSE 100
#define VECT_DCPSE_V 101
#define VECT_DCPSE_V_SUM 102
#define VECT_DCPSE_V_DOT 103
#define VECT_DCPSE_V_DIV 104
#define VECT_DCPSE_V_CURL2D 105
#define VECT_COPY_1_TO_N 300
#define VECT_COPY_N_TO_N 301
#define VECT_COPY_N_TO_1 302
#define VECT_PMUL 91
#define VECT_SUB_UNI 92




template<bool cond, typename exp1, typename exp2>
struct first_or_second
{
    typedef typename exp2::vtype vtype;

    static auto getVector(const exp1 & o1, const exp2 & o2) -> decltype(o2.getVector())
    {
        return o2.getVector();
    }
};

template<typename exp1, typename exp2>
struct first_or_second<true,exp1,exp2>
{
    typedef typename exp1::vtype vtype;

    static auto getVector(const exp1 & o1, const exp2 & o2) -> decltype(o1.getVector())
    {
        return o1.getVector();
    }
};

template<typename T, typename Sfinae = void>
struct has_vtype: std::false_type {};

/*! \brief has_data check if a type has defined a member data
 *
 * ### Example
 *
 * \snippet util_test.hpp Check has_data
 *
 * return true if T::type is a valid type
 *
 */
template<typename T>
struct has_vtype<T, typename Void<typename T::vtype>::type> : std::true_type
{};




/*! \brief has_init check if a type has defined a
 * method called init
 *
 *
 * return true if T::init() is a valid expression (function pointers)
 * and produce a defined type
 *
 */

template<typename ObjType, typename Sfinae = void>
struct has_init: std::false_type {};

template<typename ObjType>
struct has_init<ObjType, typename Void<typename ObjType::has_init>::type> : std::true_type
{};

/*! \brief Call the init function if a type T has the function init
 *
 * \tparam type T
 *
 */
template <typename T, bool has_init = has_init<T>::value >
struct call_init_if_needed
{
	//! it call the function init for r_exp if T has the function init
	static inline void call(T & r_exp)
	{
		r_exp.init();
	}
};

/*! \brief Call the init function if a type T has the function init
 *
 * \tparam type T
 *
 */
template <typename T>
struct call_init_if_needed<T,false>
{
	//! it call the function init for r_exp if T has the function init
	static inline void call(T & r_exp)
	{
	}
};



/*! \brief Unknown operation specialization
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1, typename exp2, unsigned int op>
class vector_dist_expression_op
{

};

template<typename v1_type, typename v2_type>
struct vector_result
{
	typedef v1_type type;

	template<typename exp1, typename exp2>
	static const type & getVector(const exp1 & o1, const exp2 & o2)
	{
		return o1.getVector();
	}

	template<typename exp1>
	static const type & getVector(const exp1 & o1)
	{
		return o1.getVector();
	}
};


template<typename v2_type>
struct vector_result<void,v2_type>
{
	typedef v2_type type;

	template<typename exp1, typename exp2>
	static const type & getVector(const exp1 & o1, const exp2 & o2)
	{
		return o2.getVector();
	}

	template<typename exp2>
	static const type & getVector(exp2 & o2)
	{
		return o2.getVector();
	}
};

template<typename NN1_type, typename NN2_type>
struct nn_type_result
{
	typedef NN1_type type;

	template<typename exp1, typename exp2>
	static type * getNN(exp1 & o1, exp2 & o2)
	{
		return o1.getNN();
	}

	template<typename exp1>
	static type * getNN(exp1 & o1)
	{
		return o1.getNN();
	}
};


template<typename NN2_type>
struct nn_type_result<void,NN2_type>
{
	typedef NN2_type type;

	template<typename exp1, typename exp2>
	static type * getNN(exp1 & o1, exp2 & o2)
	{
		return o2.getNN();
	}

	template<typename exp2>
	static type * getNN(exp2 & o2)
	{
		return o2.getNN();
	}
};

template<bool s1, bool s2>
struct vector_is_sort_result
{
	typedef boost::mpl::bool_<s1 | s2> type;
};

/*! \brief Sum operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,VECT_SUM>
{
	//! expression 1
	const exp1 o1;

	//! expression 2
	const exp2 o2;

public:

	//! indicate if this vector is kernel type
	typedef typename exp1::is_ker is_ker;

	//! return the vector type on which this expression operate
	typedef typename first_or_second<has_vtype<exp1>::value,exp1,exp2>::vtype vtype;

	//! result for is sort
	typedef typename vector_is_sort_result<exp1::is_sort::value,exp2::is_sort::value>::type is_sort;

	//! NN_type
	typedef typename nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::type NN_type;

	//! constructor of the expression to sum two expression
	inline vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_type * getNN() const
	{
		return nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::getNN(o1,o2);
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) + o2.value(vect_dist_key_dx()))>::type >
	__device__ __host__ inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) + o2.value(key);
	}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
	{
		o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff, comp);
		o2.template value_nz<Sys_eqs>(p_map,key,cols,coeff, comp);
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype & getVector()
	{
		return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype & getVector() const
	{
        	return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0) + o2.value(0))>::type >
	__device__ __host__ inline r_type value(const unsigned int & key) const
	{
		return o1.value(key) + o2.value(key);
	}

};

/*! \brief Subtraction operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,VECT_SUB>
{
	//! expression 1
	const exp1 o1;

	//! expression 2
	const exp2 o2;

public:

	//! The type of the internal vector
	typedef typename first_or_second<has_vtype<exp1>::value,exp1,exp2>::vtype vtype;

	typedef typename exp1::is_ker is_ker;

	//! result for is sort
	typedef typename vector_is_sort_result<exp1::is_sort::value,exp2::is_sort::value>::type is_sort;

	//! NN_type
	typedef typename nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::type NN_type;

	//! Costruct a subtraction expression out of two expressions
	inline vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_type * getNN() const
	{
		return nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::getNN(o1,o2);
	}


	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) - o2.value(vect_dist_key_dx()))>::type >
	inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) - o2.value(key);
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) - o2.value(vect_dist_key_dx()))>::type >
	__device__ __host__ inline r_type value(const unsigned int & key) const
	{
		return o1.value(key) - o2.value(key);
	}

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    const vtype & getVector()
    {
        return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
    }

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    const vtype & getVector() const
    {
        return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
    }

    template<typename Sys_eqs,typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
        o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp);
        coeff_type tmp = -coeff;
        o2.template value_nz<Sys_eqs>(p_map,key,cols,tmp,comp);
    }
};

/*! \brief Multiplication operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,VECT_MUL>
{
	//! expression 1
	const exp1 o1;

	//! expression 2
	const exp2 o2;

public:

	//! The type of the internal vector
	typedef typename first_or_second<has_vtype<exp1>::value,exp1,exp2>::vtype vtype;

	typedef typename exp1::is_ker is_ker;

	//! result for is sort
	typedef typename vector_is_sort_result<exp1::is_sort::value,exp2::is_sort::value>::type is_sort;

	//! NN_type
	typedef typename nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::type NN_type;

	//! constructor from two expressions
	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) * o2.value(vect_dist_key_dx()))>::type >
	__device__ __host__ inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) * o2.value(key);
	}

	template<typename Sys_eqs,typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
	{
        	//o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp);
    		auto coeff_tmp = o1.value(key) * coeff;
        	o2.template value_nz<Sys_eqs>(p_map,key,cols,coeff_tmp,comp);
    	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	const vtype & getVector()
	{
        	return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	*/
	const vtype & getVector() const
	{
        	return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
 	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) * o2.value(vect_dist_key_dx()))>::type >
	__device__ __host__ inline r_type value(const unsigned int & key) const
	{
		return o1.value(key) * o2.value(key);
	}
};

/*! \brief Division operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1, typename exp2>
class vector_dist_expression_op<exp1,exp2,VECT_DIV>
{
	//! expression 1
	const exp1 o1;

	//! expression 2
	const exp2 o2;

public:

	//! The type of the internal vector
	typedef typename first_or_second<has_vtype<exp1>::value,exp1,exp2>::vtype vtype;

	typedef typename exp1::is_ker is_ker;

	//! result for is sort
	typedef typename vector_is_sort_result<exp1::is_sort::value,exp2::is_sort::value>::type is_sort;

	//! NN_type
	typedef typename nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::type NN_type;

	//! constructor from two expressions
	vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_type * getNN() const
	{
		return nn_type_result<typename exp1::NN_type,typename exp2::NN_type>::getNN(o1,o2);
	}


	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) / o2.value(vect_dist_key_dx()))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) / o2.value(key);
	}

    template<typename Sys_eqs,typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
    	std::cout << __FILE__ << ":" << __LINE__ << " You are trying to divide by an operator,  this is not possible " << std::endl;
    }

    /*! \brief Return the vector on which is acting
 *
 * It return the vector used in getVExpr, to get this object
 *
 * \return the vector
 *
 */
    const vtype & getVector()
    {
        return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
    }

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    const vtype & getVector() const
    {
        return first_or_second<has_vtype<exp1>::value,exp1,exp2>::getVector(o1,o2);
    }

        /*! \brief Evaluate the expression
         *
         * \param key where to evaluate the expression
         *
         * \return the result of the expression
         *
         */
        template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) / o2.value(vect_dist_key_dx()))>::type >
        __device__ __host__ inline r_type value(const unsigned int & key) const
        {
                return o1.value(key) / o2.value(key);
        }
};


/*! \brief it take an expression and create the negatove of this expression
 *
 *
 */
template <typename exp1>
class vector_dist_expression_op<exp1,void,VECT_SUB_UNI>
{
	//! expression 1
	const exp1 o1;

public:


	typedef typename exp1::is_ker is_ker;

	//! return the vector type on which this expression operate
	typedef typename vector_result<typename exp1::vtype,void>::type vtype;

	//! result for is sort
	typedef typename vector_is_sort_result<exp1::is_sort::value,false>::type is_sort;

	//! NN_type
	typedef typename nn_type_result<typename exp1::NN_type,void>::type NN_type;

	//! constructor from an expresssion
	vector_dist_expression_op(const exp1 & o1)
	:o1(o1)
	{}

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_type * getNN() const
	{
		return nn_type_result<typename exp1::NN_type,void>::getNN(o1);
	}

	/*! \brief Return the underlying vector
	 *
	 * \return the vector
	 *
	 */
	const vtype & getVector()
	{
		return vector_result<typename exp1::vtype,void>::getVector(o1);
	}

	//! initialize the expression tree
	inline void init() const
	{
		o1.init();
	}

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    const vtype & getVector() const
    {
        return o1.getVector();
    }


	//! return the result of the expression
	template<typename r_type=typename std::remove_reference<decltype(-(o1.value(vect_dist_key_dx(0))))>::type >
	__device__ __host__  inline r_type value(const vect_dist_key_dx & key) const
	{
		return -(o1.value(key));
	}

    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
	    coeff_type coeff_tmp = -coeff;
        o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff_tmp, comp);
    }
	//! return the result of the expression
	template<typename r_type=typename std::remove_reference<decltype(-(o1.value(vect_dist_key_dx(0))))>::type >
	__device__ __host__ inline r_type value(const unsigned int & key) const
	{
		return -(o1.value(key));
	}
};

/*! \brief Expression implementation computation selector
 *
 *
 */
template<int impl, bool vect_ker>
struct vector_dist_expression_comp_sel
{
	typedef boost::mpl::int_<impl> type;
};

template<>
struct vector_dist_expression_comp_sel<comp_host,true>
{
	typedef boost::mpl::int_<-1> type;
};

template<>
struct vector_dist_expression_comp_sel<comp_dev,false>
{
	typedef boost::mpl::int_<-1> type;
};

/*! \brief Expression implementation computation selector
 *
 */
template<bool cond>
struct vector_dist_expression_comp_proxy_sel
{
    template<bool cond_, typename v_type, typename exp_type>
    static void compute(v_type &v,exp_type &v_exp)
    { vector_dist_op_compute_op<0,false,vector_dist_expression_comp_sel<comp_dev,cond_>::type::value>
        ::compute_expr(v,v_exp);}
};
template<>
struct vector_dist_expression_comp_proxy_sel<false>
{
    template<bool cond, typename v_type, typename exp_type>
    static void compute(v_type &v, exp_type &v_exp)
    {   auto v_ker=v.toKernel();
        vector_dist_op_compute_op<0,false,vector_dist_expression_comp_sel<comp_dev,cond>::type::value>
        ::compute_expr(v_ker,v_exp);}
};

template<typename vector, bool is_ker = has_vector_kernel<vector>::type::value>
struct vector_expression_transform
{
	typedef vector& type;
};

template<typename vector>
struct vector_expression_transform<vector,true>
{
	typedef vector type;
};

template<typename vector>
struct v_mem_mutable
{
	vector v;

	v_mem_mutable(vector & v)
	:v(v)
	{}
};



/*! \brief Main class that encapsulate a vector properties operand to be used for expressions construction
 *
 * \tparam prp property involved
 * \tparam vector involved
 *
 */
template<unsigned int prp, typename vector>
class vector_dist_expression
{
	//! The vector
	mutable v_mem_mutable<typename vector_expression_transform<vector>::type> v;

	vector_dist_ker_list<vector> * vdl;

public:

	typedef typename has_vector_kernel<vector>::type is_ker;

	//! The type of the internal vector
	typedef vector vtype;

	//! result for is sort
	typedef boost::mpl::bool_<false> is_sort;

	//! NN_type
	typedef void NN_type;

	//! Property id of the point
	static const unsigned int prop = prp;

	int var_id = 0;

	void setVarId(int var_id)
	{
		this->var_id = var_id;
	}

	//! constructor for an external vector
	vector_dist_expression(vector & v)
	:v(v),vdl(NULL)
	{}

	//! constructor for an external vector
	~vector_dist_expression()
	{
		if (vdl != NULL)
		{vdl->remove(v.v);}
	}

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline void * getNN() const
	{
		return NULL;
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	__device__ __host__ const vector & getVector() const
	{
		return v.v;
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	__device__ __host__ vector & getVector()
	{
		return v.v;
	}

	/*! \brief set vector_dist_ker_list
	 *
	 * \param vdkl vector_dist_ker_list
	 *
	 */
	void set_vector_dist_ker_list(vector_dist_ker_list<vector> & vdkl, bool is_sort)
	{
		vdkl.add(v.v,is_sort);
		vdl = &vdkl;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param k where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__host__ __device__ inline auto value(const vect_dist_key_dx & k) const -> decltype(pos_or_propR<vector,prp>::value(v.v,k))
	{
		return pos_or_propR<vector,prp>::value(v.v,k);
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__host__ __device__ inline auto value(const vect_dist_key_dx & k) -> decltype(pos_or_propR<vector,prp>::value(v.v,k))
	{
		return pos_or_propR<vector,prp>::value(v.v,k);
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__device__ __host__ inline auto value(const unsigned int & k) const -> decltype(pos_or_propR<vector,prp>::value(v.v,k))
	{
		return pos_or_propR<vector,prp>::value(v.v,k);
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__device__ __host__ inline auto value(const unsigned int & k) -> decltype(pos_or_propR<vector,prp>::value(v.v,k))
	{
		return pos_or_propR<vector,prp>::value(v.v,k);
	}

    inline auto get(const vect_dist_key_dx & key) const -> decltype(value(key))
    {
        return this->value(key);
    }

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<unsigned int prp2> vector & operator=(const vector_dist_expression<prp2,vector> & v_exp)
	{
        if (v_exp.getVector().isSubset() == true)
        {
                        std::cout << __FILE__ << ":" << __LINE__ << " error on the right hand side of the expression you have to use non-subset properties" << std::endl;
                        return v.v;
        }

		if (has_vector_kernel<vector>::type::value == false)
		{
			vector_dist_op_compute_op<prp,false,vector_dist_expression_comp_sel<comp_host,
																	   	  has_vector_kernel<vector>::type::value>::type::value>
			::compute_expr(v.v,v_exp);
		}
		else
		{
			vector_dist_op_compute_op<prp,false,vector_dist_expression_comp_sel<comp_dev,
		   	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  has_vector_kernel<vector>::type::value>::type::value>
			::compute_expr(v.v,v_exp);
		}

		return v.v;
	}

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<typename T,typename memory,template <typename> class layout_base > vector & operator=(const vector_dist_expression<0,openfpm::vector<aggregate<T>, memory, layout_base>> & v_exp)
	{
		//vector_dist_op_compute_op<prp,false,vector_dist_expression_comp_sel<comp_host,has_vector_kernel<vector>::type::value>::type::value>
		//::compute_expr(v.v,v_exp);
            vector_dist_op_compute_op<prp,false,vector_dist_expression_comp_sel<comp_host,
                    has_vector_kernel<vector>::type::value>::type::value>
            ::compute_expr(v.v,v_exp);

        
		return v.v;
	}

    /*! \brief Fill the vector property with the evaluated expression
 *
 * \param v_exp expression to evaluate
 *
 * \return itself
 *
 */
    template<typename T> vector & operator=(const vector_dist_expression<0,openfpm::vector_gpu<aggregate<T>>> & v_exp)
    {
            vector_dist_op_compute_op<prp,false,vector_dist_expression_comp_sel<comp_dev,
                    has_vector_kernel<vector>::type::value>::type::value>
            ::compute_expr(v.v,v_exp);

        return v.v;
    }

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<typename exp1, typename exp2, unsigned int op>
	vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
	{
        if (v_exp.getVector().isSubset() == true)
        {
        	std::cout << __FILE__ << ":" << __LINE__ << " error on the right hand side of the expression you have to use non-subset properties" << std::endl;
            return v.v;
        }

		if (has_vector_kernel<vector>::type::value == false)
		{
			vector_dist_op_compute_op<prp,
									  vector_dist_expression_op<exp1,exp2,op>::is_sort::value,
									  vector_dist_expression_comp_sel<comp_host,
																	  has_vector_kernel<vector>::type::value>::type::value>
			::compute_expr(v.v,v_exp);
		}
		else
		{
			vector_dist_op_compute_op<prp,
									  vector_dist_expression_op<exp1,exp2,op>::is_sort::value,
									  vector_dist_expression_comp_sel<comp_dev,
		   	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  has_vector_kernel<vector>::type::value>::type::value>
			::compute_expr(v.v,v_exp);
		}

		return v.v;
	}

	/*! \brief Fill the vector property with the double
	 *
	 * \param d value to fill
	 *
	 * \return the internal vector
	 *
	 */
	vector & operator=(double d)
	{
		if (has_vector_kernel<vector>::type::value == false)
		{
			vector_dist_op_compute_op<prp,
									  false,
									  vector_dist_expression_comp_sel<comp_host,
																	  has_vector_kernel<vector>::type::value>::type::value>
			::compute_const(v.v,d);
		}
		else
		{
			vector_dist_op_compute_op<prp,
									  false,
									  vector_dist_expression_comp_sel<comp_dev,
		   	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  has_vector_kernel<vector>::type::value>::type::value>
			::compute_const(v.v,d);
		}

		return v.v;
	}


    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
            cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + var_id + comp] += coeff;
    }

    inline vector_dist_expression_op<vector_dist_expression<prp,vector>,boost::mpl::int_<1>,VECT_COMP> operator[](int comp)
    {
        int comp_n[1];

        comp_n[0] = comp;

        vector_dist_expression_op<vector_dist_expression<prp,vector>,boost::mpl::int_<1>,VECT_COMP> v_exp(*this,comp_n,var_id);

        return v_exp;
    }
};

/*! \brief Main class that encapsulate a vector properties operand to be used for expressions construction
 *  Temporal Expressions
 * \tparam prp property involved
 * \tparam vector involved
 *
 */
template<typename vector_type>
class vector_dist_expression_impl
{
    //! Internal vector
    typedef vector_type vector;
    typedef typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<0>>::type T;
    //! The temporal vector
    mutable vector v;

public:

    typedef T * iterator;
    typedef const  T * const_iterator;

    typedef typename has_vector_kernel<vector>::type is_ker;

    //! The type of the internal vector
    typedef vector vtype;

    //! The type of the internal value
    typedef T value_type;

    //! result for is sort
    typedef boost::mpl::bool_<false> is_sort;

    //! NN_type
    typedef void NN_type;

    //! Property id of the point
    static const unsigned int prop = 0;

    int var_id = 0;

    void setVarId(int var_id)
    {
        this->var_id = var_id;
    }

    ///////// BOOST ODEINT interface
    iterator begin()
    { return &v.template get<0>(0); }

    const_iterator begin() const
    { return &v.template get<0>(0); }

    iterator end()
    { return &v.template get<0>(v.size()-1)+1; }

    const_iterator end() const
    { return &v.template get<0>(v.size()-1)+1; }

    size_t size() const
    { return v.size(); }

    void resize(size_t n)
    {
        // Here

        v.resize(n);
    }

/*	T * begin() {
        return &v.template get<0>(0);
    }

    T * end() {
	    return &v.template get<0>(v.size()-1);
	}*/

    // ... [ implement container interface ]
//]
    //const double & operator[]( const size_t n ) const
    //{ return m_v[n]; }

    //double & operator[]( const size_t n )
    //{ return m_v[n]; }


    ////////////////////////////////////

    vector_dist_expression_impl()
    {}

    template<typename exp1, typename exp2, unsigned int op>
    vector_dist_expression_impl(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
    {
        this->operator=(v_exp);
    }

    /*! \brief get the NN object
     *
     * \return the NN object
     *
     */
    inline void * getNN() const
    {
        return NULL;
    }

    /*! \brief Return the vector on which is acting
     *
     * It return the vector used in getVExpr, to get this object
     *
     * \return the vector
     *
     */
    __device__ __host__ const vector & getVector() const
    {
        return v;
    }

    /*! \brief Return the vector on which is acting
     *
     * It return the vector used in getVExpr, to get this object
     *
     * \return the vector
     *
     */
    __device__ __host__ vector & getVector()
    {
        return v;
    }

    /*! \brief This function must be called before value
     *
     * it initialize the expression if needed
     *
     */
    inline void init() const
    {}

    /*! \brief Evaluate the expression
     *
     * \param k where to evaluate the expression
     *
     * \return the result of the expression
     *
     */
    __host__ __device__ inline auto value(const vect_dist_key_dx & k) const -> decltype(v.template get<0>(k.getKey()))
    {
        return v.template get<0>(k.getKey());
    }


    /*! \brief Fill the vector property with the evaluated expression
     *
     * \param v_exp expression to evaluate
     *
     * \return itself
     *
     */
    template<unsigned int prp2, typename vector2> vector & operator=(const vector_dist_expression<prp2,vector2> & v_exp)
    {
        if (v_exp.getVector().isSubset() == true)
        {
            std::cout << __FILE__ << ":" << __LINE__ << " error on the right hand side of the expression you have to use non-subset properties" << std::endl;
            return v;
        }

        v.resize(v_exp.getVector().size_local());
        constexpr bool cond=has_vector_kernel<vector>::type::value || std::is_same<vector,openfpm::vector<aggregate<T>,CudaMemory,memory_traits_inte>>::value;
        //std::cout<<cond<<std::endl;
        //std::cout<< (vector_dist_expression_comp_sel<comp_host,has_vector_kernel<vector>::type::value>::type::value || std::is_same<vector,openfpm::vector<aggregate<T>,CudaMemory,memory_traits_inte>>::value)<<std::endl;
        //std::cout<<(vector_dist_expression_comp_sel<2,
           //     has_vector_kernel<vector>::type::value>::type::value || std::is_same<vector,openfpm::vector<aggregate<T>,CudaMemory,memory_traits_inte>>::value)<<std::endl;
        //std::cout<<has_vector_kernel<vector>::type::value<<std::endl;
        //std::cout<<vector_dist_expression_comp_sel<2,false>::type::value<<std::endl;
        //std::cout<<!std::is_same<vector,openfpm::vector<aggregate<T>,CudaMemory,memory_traits_inte>>::value<<std::endl;
        if (has_vector_kernel<vector>::type::value == false && !std::is_same<vector,openfpm::vector<aggregate<T>,CudaMemory,memory_traits_inte>>::value)
        {
            vector_dist_op_compute_op<0,false,vector_dist_expression_comp_sel<comp_host,cond>::type::value>
            ::compute_expr(v,v_exp);
        }
        else
        {
            vector_dist_expression_comp_proxy_sel<!std::is_same<vector,openfpm::vector<aggregate<T>,CudaMemory,memory_traits_inte>>::value>::template compute<cond>(v,v_exp);
        }

        return v;
    }


    /*! \brief Fill the vector property with the evaluated expression
     *
     * \param v_exp expression to evaluate
     *
     * \return itself
     *
     */
    template<typename exp1, typename exp2, unsigned int op>
    vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
    {
        if (v_exp.getVector().isSubset() == true)
        {
            std::cout << __FILE__ << ":" << __LINE__ << " error on the right hand side of the expression you have to use non-subset properties" << std::endl;
            return v;
        }

        v.resize(v_exp.getVector().size_local());

        if (has_vector_kernel<vector>::type::value == false)
        {
            vector_dist_op_compute_op<0,
                    vector_dist_expression_op<exp1,exp2,op>::is_sort::value,
                    vector_dist_expression_comp_sel<comp_host,
                            has_vector_kernel<vector>::type::value>::type::value>
            ::compute_expr(v,v_exp);
        }
        else
        {
            vector_dist_op_compute_op<0,
                    vector_dist_expression_op<exp1,exp2,op>::is_sort::value,
                    vector_dist_expression_comp_sel<comp_dev,
                            has_vector_kernel<vector>::type::value>::type::value>
            ::compute_expr(v,v_exp);
        }

        return v;
    }

    /*! \brief Fill the vector property with the double
     *
     * \param d value to fill
     *
     * \return the internal vector
     *
     */
    vector & operator=(double d)
    {
        std::cout << __FILE__ << ":" << __LINE__ << " Error: temporal with constants is unsupported" << std::endl;
    }


    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
        std::cout << __FILE__ << ":" << __LINE__ << " Error: use of temporal is not supported to construct equations";
    }

    inline vector_dist_expression_op<vector_dist_expression<0,vector>,boost::mpl::int_<1>,VECT_COMP> operator[](int comp)
    {
        int comp_n[1];

        comp_n[0] = comp;

        vector_dist_expression_op<vector_dist_expression<0,vector>,boost::mpl::int_<1>,VECT_COMP> v_exp(*this,comp_n,var_id);

        return v_exp;
    }
};

/*! \brief Sub class that encapsulate a vector properties operand to be used for expressions construction
 *  Temporal Expressions
 * \tparam prp property involved
 * \tparam vector involved
 *
 */
template<typename T, typename memory,template <typename> class layout_base >
class vector_dist_expression<0,openfpm::vector<aggregate<T>,memory, layout_base> > : public vector_dist_expression_impl<openfpm::vector<aggregate<T>,memory, layout_base>>
{
    typedef openfpm::vector<aggregate<T>,memory, layout_base> vector;
    typedef vector_dist_expression_impl<vector> base;
public:
    template<unsigned int prp2, typename vector2> vector & operator=(const vector_dist_expression<prp2,vector2> & v_exp)
    {
        return base::operator=(v_exp);
    }
    template<typename exp1, typename exp2, unsigned int op>
    vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
    {
        return base::operator=(v_exp);
    }
};

/*! \brief Sub class that encapsulate a GPU vector properties operand to be used for expressions construction
 *  Temporal Expressions
 * \tparam prp property involved
 * \tparam vector involved
 *
 */
template<typename T>
class vector_dist_expression<0,openfpm::vector_gpu<aggregate<T>>> : public vector_dist_expression_impl<openfpm::vector_gpu<aggregate<T>>>
{
    typedef openfpm::vector_gpu<aggregate<T>> vector;
    typedef vector_dist_expression_impl<vector> base;
public:
    template<unsigned int prp2, typename vector2> vector & operator=(const vector_dist_expression<prp2,vector2> & v_exp)
    {
        return base::operator=(v_exp);
    }
    template<typename exp1, typename exp2, unsigned int op>
    vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
    {
        return base::operator=(v_exp);
    }

};

template<typename T> using texp_v = vector_dist_expression<0,openfpm::vector<aggregate<T>>>;
template<typename T> using texp_v_gpu = vector_dist_expression<0,openfpm::vector_gpu<aggregate<T>>>;


template<typename vector, unsigned int impl>
struct switcher_get_v
{
	typedef vector type;

	static vector & get(vector & v)	{return v;};

	template<typename exp_type, typename vector_klist>
	static void register_vector(exp_type & exp_v, vector_klist & v,bool is_sort)
	{
	}
};

template<typename vector>
struct switcher_get_v<vector,comp_dev>
{
	typedef decltype(std::declval<vector>().toKernel()) type;

	static type get(vector & v)	{return v.toKernel();};

	template<typename exp_type, typename vector_klist>
	static void register_vector(exp_type & exp_v, vector_klist & v,bool is_sort)
	{
		exp_v.set_vector_dist_ker_list(v.private_get_vector_dist_ker_list(),is_sort);
	}
};

/*! \brief it take an expression and create the negatove of this expression
 *
 *
 */
template <typename exp1,int n>
class vector_dist_expression_op<exp1,boost::mpl::int_<n>,VECT_COMP>
{
	//! expression 1
	exp1 o1;

	//! component
	int comp[n];

	int var_id = 0;
    void setVarId(int var_id)
    {
        this->var_id = var_id;
    }

	typedef vector_dist_expression_op<exp1,boost::mpl::int_<n>,VECT_COMP> myself;

public:

        typedef std::false_type is_ker;

        //! result for is sort
        typedef boost::mpl::bool_<false> is_sort;

        //! result for is sort
        typedef boost::mpl::bool_<false> NN_type;

	typedef typename exp1::vtype vtype;

	//! constructor from an expresssion

	vector_dist_expression_op(const exp1 & o1, int (& comp)[n], int var_id)
	:o1(o1),var_id(var_id)
	{
		for (int i = 0 ; i < n ; i++)
		{this->comp[i] = comp[i];}
	}

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    __device__ __host__ const vtype & getVector() const
    {
        return o1.getVector();
    }

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    __device__ __host__ vtype & getVector()
    {
        return o1.getVector();
    }

	//! initialize the expression tree
	inline void init() const
	{
		o1.init();
	}

	//! property on which this view is acting
	//typedef typename boost::mpl::at<typename vtype::value_type::type,boost::mpl::int_<exp1::prop>>::type property_act;
	typedef typename pos_or_propL<vtype,exp1::prop>::property_act property_act;

	/*! \brief Return the result of the expression
	 *
	 * \note this function must be deactivated on transitional objects. Suppose we are slicing a tensor of rank 2
	 *            an object of rank 1 is implicitly created for such object we have to deactivate this function
	 *            because ill-formed
	 *
	 * \param key point where to evaluate
	 *
	 *
	 */
	__host__ __device__ inline auto value(const vect_dist_key_dx & key) const -> decltype(get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,vect_dist_key_dx(0),comp))
	{
		return get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,key,comp);
	}

	/*! \brief Return the result of the expression
	 *
	 * \note this function must be deactivated on transitional objects. Suppose we are slicing a tensor of rank 2
	 *            an object of rank 1 is implicitly created for such object we have to deactivate this function
	 *            because ill-formed
	 *
	 * \param key point where to evaluate
	 *
	 *
	 */
	__host__ __device__ inline auto value(const vect_dist_key_dx & key) -> decltype(get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,vect_dist_key_dx(0),comp))
	{
		return get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,key,comp);
	}

	/*! \brief Return the result of the expression
	 *
	 * \note this function must be deactivated on transitional objects. Suppose we are slicing a tensor of rank 2
	 *            an object of rank 1 is implicitly created for such object we have to deactivate this function
	 *            because ill-formed
	 *
	 * \param key point where to evaluate
	 *
	 *
	 */
	inline auto get(const vect_dist_key_dx & key) const -> decltype(value(key))
	{
		return this->value(key);
	}

    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp_) const
    {
#ifdef SE_CLASS1

    	if (n != 1)
    	{
    		std::cout << __FILE__ << ":" << __LINE__ << " Error it only work for tensore of rank 1 ... like vectors " << std::endl;
    	}

#endif

        o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp_ + var_id + comp[0]);
    }

    inline vector_dist_expression_op<exp1,boost::mpl::int_<2>,VECT_COMP> operator[](int comp_)
    {
    	int comp_n[n+1];

    	for (int i = 0 ; i < n ; i++)
    	{comp_n[i] = comp[i];}
    	comp_n[n] = comp_;

    	vector_dist_expression_op<exp1,boost::mpl::int_<2>,VECT_COMP> v_exp(o1,comp_n,var_id);

    	return v_exp;
    }

    /*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
    template<typename T, typename memory> vtype & operator=(const vector_dist_expression<0,openfpm::vector<aggregate<T>,memory>> & v_exp)
    {
        v_exp.init();

        auto & v = getVector();
/*#ifdef SE_CLASS1
        auto &v2=v_exp.getVector();

        SubsetSelector_impl<std::remove_reference<decltype(v)>::type::is_it_a_subset::value>::check(v2,v);
#endif*/

        auto it = v.getDomainIterator();

        while (it.isNext())
        {
            auto key = it.get();
            auto key_orig = v.getOriginKey(key);

            get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<exp1::prop>(v_exp,v,key,key_orig,comp);

            ++it;
        }

        return v;
    }

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<unsigned int prp2> vtype & operator=(const vector_dist_expression<prp2,vtype> & v_exp)
	{
		v_exp.init();

		auto & v = getVector();
#ifdef SE_CLASS1
		auto &v2=v_exp.getVector();

        SubsetSelector_impl<std::remove_reference<decltype(v)>::type::is_it_a_subset::value>::check(v2,v);
#endif
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();
			auto key_orig = v.getOriginKey(key);

			get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<exp1::prop>(v_exp,v,key,key_orig,comp);

			++it;
		}

		return v;
	}

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<typename exp1_, typename exp2_, unsigned int op> vtype & operator=(const vector_dist_expression_op<exp1_,exp2_,op> & v_exp)
	{
        if (v_exp.getVector().isSubset() == true)
        {
            std::cout << __FILE__ << ":" << __LINE__ << " error on the right hand side of the expression you have to use non-subset properties" << std::endl;
            return this->getVector();
        }

		if (has_vector_kernel<vtype>::type::value == false)
		{
			vector_dist_op_compute_op<exp1::prop,false,vector_dist_expression_comp_sel<comp_host,
																	   	  has_vector_kernel<vtype>::type::value>::type::value>
			::compute_expr_slice(o1.getVector(),v_exp,comp);
		}
		else
		{
			vector_dist_op_compute_op<exp1::prop,false,vector_dist_expression_comp_sel<comp_dev,
		   	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  has_vector_kernel<vtype>::type::value>::type::value>
			::compute_expr_slice(o1.getVector(),v_exp,comp);
		}

		return this->getVector();
	}

	/*! \brief Fill the vector property with the double
	 *
	 * \param d value to fill
	 *
	 * \return the internal vector
	 *
	 */
	vtype & operator=(double d)
	{
		auto & v = getVector();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			//pos_or_propL<vtype,exp1::prp>::value(v,key) = d;
			get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign_double<exp1::prop>(d,v,key,comp);


			++it;
		}

		return v;
	}
};

/*! \Create an expression from a vector property
 *
 * \tpatam prp property
 * \param v
 *
 */
template <unsigned int prp,unsigned int impl = comp_host, typename vector>
inline vector_dist_expression<prp,typename switcher_get_v<vector,impl>::type > getV(vector & v)
{
	decltype(switcher_get_v<vector,impl>::get(v)) vk = switcher_get_v<vector,impl>::get(v);
	vector_dist_expression<prp,typename switcher_get_v<vector,impl>::type > exp_v(vk);

	switcher_get_v<vector,impl>::register_vector(exp_v,v,false);

	return exp_v;
}


/*! \Create an expression from a vector property
 *
 * \tpatam prp property
 * \param v
 *
 */
template <unsigned int prp,typename vector>
inline vector_dist_expression<prp, typename switcher_get_v<vector,comp_dev>::type > getV_sort(vector & v)
{
	auto vk = v.toKernel_sorted();
	vector_dist_expression<prp,typename switcher_get_v<vector, comp_dev>::type > exp_v(vk);

	exp_v.set_vector_dist_ker_list(v.private_get_vector_dist_ker_list(),true);

	return exp_v;
}

/*! \brief Main class that encapsulate a double constant
 *
 * \param prp no meaning
 *
 */
template<unsigned int prp>
class vector_dist_expression<prp,double>
{
	//! constant parameter
	double d;

public:

	typedef std::false_type is_ker;

	//! result for is sort
	typedef boost::mpl::bool_<false> is_sort;

	typedef void NN_type;

	//! constructor from a constant expression
	inline vector_dist_expression(const double & d)
	:d(d)
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
	 * \param k ignored position in the vector
	 *
	 * It just return the value set in the constructor
	 *
	 * \return the constant value
	 *
	 */
	__device__ __host__ inline double value(const vect_dist_key_dx & k) const
	{
		return d;
	}


    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
        cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += coeff;
    }
	/*! \brief Evaluate the expression
	 *
	 * \param k ignored position in the vector
	 *
	 * It just return the value set in the constructor
	 *
	 * \return the constant value
	 *
	 */
	__host__ __device__ inline double value(const unsigned int & k) const
	{
		return d;
	}
};


/*! \brief Main class that encapsulate a float constant
 *
 * \param prp no meaning
 *
 */
template<unsigned int prp>
class vector_dist_expression<prp,float>
{
	//! constant value
	float d;

public:

	typedef std::false_type is_ker;


	//! result for is sort
	typedef boost::mpl::bool_<false> is_sort;

	typedef void NN_type;

	//! constrictor from constant value
	inline vector_dist_expression(const float & d)
	:d(d)
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
	 * \param k ignored position in the vector
	 *
	 * It just return the value set in the constructor
	 *
	 * \return the constant value set in the constructor
	 *
	 */
	inline float value(const vect_dist_key_dx & k) const
	{
		return d;
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k ignored position in the vector
	 *
	 * It just return the value set in the constructor
	 *
	 * \return the constant value set in the constructor
	 *
	 */
	__device__ __host__ inline float value(const unsigned int & k) const
	{
		return d;
	}
};

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, unsigned int p2, typename v1, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,VECT_SUM>
operator+(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,VECT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1, unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,VECT_SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,VECT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1, unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,VECT_SUM>
operator+(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,VECT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,VECT_SUM>
operator+(const vector_dist_expression<prp1,v1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,VECT_SUM> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, unsigned int prp1, typename v1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,float>,VECT_SUM>
operator+(const vector_dist_expression<prp1,v1> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,float>,VECT_SUM> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_SUM>
operator+(double d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_SUM> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<prp1,v1>,VECT_SUM>
operator+(float d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<prp1,v1>,VECT_SUM> exp_sum(vector_dist_expression<0,float>(d),vb);

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_SUM> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief sum two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, typename exp1 , typename exp2, unsigned int op1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_SUM> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, unsigned int p2, typename v1, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,VECT_SUB>
operator-(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,VECT_SUB> exp_sum(va,vb);

	return exp_sum;
}


/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p2,v2>,VECT_SUB>
operator-(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p2,v2>,VECT_SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief minus of a distributed vector expression operator
 *
 * \param va vector expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2_, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,void,VECT_SUB_UNI>
operator-(const vector_dist_expression_op<exp1,exp2_,op1> & va)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2_,op1>,void,VECT_SUB_UNI> exp_sum(va);

	return exp_sum;
}

/* \brief minus of a distributed vector expression
 *
 * \param va vector expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,void,VECT_SUB_UNI>
operator-(const vector_dist_expression<p1,v1> & va)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,void,VECT_SUB_UNI> exp_sum(va);

	return exp_sum;
}


/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression_op<exp1,exp2,op1>,VECT_SUB>
operator-(const vector_dist_expression<p2,v2> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p2,v2>, vector_dist_expression_op<exp1,exp2,op1>,VECT_SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename exp3, typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_SUB>
operator-(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,VECT_SUB>
operator-(const vector_dist_expression<prp1,v1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,VECT_SUB> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
//template<unsigned int prp1, typename v1>
template<typename T, unsigned int prp1,typename v1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,float>,VECT_SUB>
operator-(const vector_dist_expression<prp1,v1> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,float>,VECT_SUB> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_SUB>
operator-(double d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_SUB> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief subtract two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, unsigned int prp1,typename v1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<prp1,v1>,VECT_SUB>
operator-(T d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<prp1,v1>,VECT_SUB> exp_sum(vector_dist_expression<0,float>(d),vb);

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<p2,v2>,VECT_MUL>
operator*(double d, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<p2,v2>,VECT_MUL> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, unsigned int p2,typename v2, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<p2,v2>,VECT_MUL>
operator*(T d, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<p2,v2>,VECT_MUL> exp_sum(vector_dist_expression<0,float>(d),vb);

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression<0,double>,VECT_MUL>
operator*(const vector_dist_expression<p2,v2> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression<0,double>,VECT_MUL> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, unsigned int p2,typename v2, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression<0,float>,VECT_MUL>
operator*(const vector_dist_expression<p2,v2> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression<p2,v2>,vector_dist_expression<0,float>,VECT_MUL> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1,unsigned int p2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,VECT_MUL>
operator*(const vector_dist_expression<p1,v1> & va, const vector_dist_expression<p2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression<p2,v2>,VECT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1, typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression_op<exp1,exp2,op1>,VECT_MUL>
operator*(const vector_dist_expression<p1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<p1,v1>,vector_dist_expression_op<exp1,exp2,op1>,VECT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1, typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p1,v1>,VECT_MUL>
operator*(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<p1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<p1,v1>,VECT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_MUL>
operator*(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply a distributed vector expression by a number
 *
 * \param va vector expression
 * \param d number
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_MUL>
operator*(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_MUL> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief Multiply a distributed vector expression by a number
 *
 * \param va vector expression
 * \param d number
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, typename exp1 , typename exp2, unsigned int op1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_MUL>
operator*(const vector_dist_expression_op<exp1,exp2,op1> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_MUL> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief Multiply a distributed vector expression by a number
 *
 * \param d number
 * \param vb vector expression
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression_op<exp1,exp2,op1>,VECT_MUL>
operator*(double d, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression_op<exp1,exp2,op1>,VECT_MUL> exp_sum(vector_dist_expression<0,double>(d),vb);

	return exp_sum;
}

/* \brief Multiply a distributed vector expression by a number
 *
 * \param d number
 * \param vb vector expression
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, typename exp1 , typename exp2, unsigned int op1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression_op<exp1,exp2,op1>,VECT_MUL>
operator*(T d, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression_op<exp1,exp2,op1>,VECT_MUL> exp_sum(vector_dist_expression<0,float>(d),vb);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_DIV> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, typename exp1 , typename exp2, unsigned int op1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_DIV> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_DIV>
operator/(double d, const vector_dist_expression_op<exp1,exp2,op1> & va)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_DIV> exp_sum(vector_dist_expression<0,double>(d),va);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, typename exp1 , typename exp2, unsigned int op1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_DIV>
operator/(T d, const vector_dist_expression_op<exp1,exp2,op1> & va)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,float>,VECT_DIV> exp_sum(vector_dist_expression<0,float>(d),va);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,VECT_DIV>
operator/(const vector_dist_expression<prp1,v1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,double>,VECT_DIV> exp_sum(va,vector_dist_expression<0,double>(d));

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, unsigned int prp1,typename v1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,float>,VECT_DIV>
operator/(const vector_dist_expression<prp1,v1> & va, T d)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<0,float>,VECT_DIV> exp_sum(va,vector_dist_expression<0,float>(d));

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_DIV>
operator/(double d, const vector_dist_expression<prp1,v1> & va)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_DIV> exp_sum(vector_dist_expression<0,double>(d),va);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, unsigned int prp1,typename v1, typename sfinae = typename std::enable_if<std::is_same<T,float>::value>::type >
inline vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<prp1,v1>,VECT_DIV>
operator/(T d, const vector_dist_expression<prp1,v1> & va)
{
	vector_dist_expression_op<vector_dist_expression<0,float>,vector_dist_expression<prp1,v1>,VECT_DIV> exp_sum(vector_dist_expression<0,float>(d),va);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1, unsigned int prp2, typename v2>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<prp2,v2>,VECT_DIV>
operator/(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression<prp2,v2> & vb)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression<prp2,v2>,VECT_DIV> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1, typename exp1,typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,VECT_DIV>
operator/(const vector_dist_expression<prp1,v1> & va, const vector_dist_expression_op<exp1,exp2,op1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,vector_dist_expression_op<exp1,exp2,op1>,VECT_DIV> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1, typename v1, typename exp1,typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,VECT_DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<prp1,v1>,VECT_DIV> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1,typename exp2, unsigned int op1, typename exp3, typename exp4, unsigned int op2>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_DIV>
operator/(const vector_dist_expression_op<exp1,exp2,op1> & va, const vector_dist_expression_op<exp3,exp4,op2> & vb)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression_op<exp3,exp4,op2>,VECT_DIV> exp_sum(va,vb);

	return exp_sum;
}

#include "vector_dist_operators_apply_kernel.hpp"
#include "vector_dist_operators_functions.hpp"
#include "vector_dist_operators_extensions.hpp"
#include "Operators/Vector/vector_dist_operator_assign.hpp"


#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_ */
