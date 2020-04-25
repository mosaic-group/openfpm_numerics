/*
* vector_dist_operators.hpp
*
*  Created on: Jun 11, 2016
*      Author: i-bird
*/

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_

#include "Vector/vector_dist.hpp"

#define PROP_POS (unsigned int)-1
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

    //! The type of the internal vector
    typedef typename first_or_second<has_vtype<exp1>::value,exp1,exp2>::vtype vtype;

	//! constructor of the expression to sum two expression
	inline vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
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
	 * \return return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) + o2.value(vect_dist_key_dx()))>::type >
	inline r_type value(const vect_dist_key_dx & key) const
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
    vtype & getVector()
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

	//! Costruct a subtraction expression out of two expressions
	inline vector_dist_expression_op(const exp1 & o1, const exp2 & o2)
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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) - o2.value(vect_dist_key_dx()))>::type > inline r_type value(const vect_dist_key_dx & key) const
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
    vtype & getVector()
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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) * o2.value(vect_dist_key_dx()))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) * o2.value(key);
	}

    template<typename Sys_eqs,typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
        o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp);
        o2.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp);
    }

    /*! \brief Return the vector on which is acting
 *
 * It return the vector used in getVExpr, to get this object
 *
 * \return the vector
 *
 */
    vtype & getVector()
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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()) / o2.value(vect_dist_key_dx()))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) / o2.value(key);
	}

    template<typename Sys_eqs,typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
        o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp);
        o2.template value_nz<Sys_eqs>(p_map,key,cols,coeff,comp);
    }

    /*! \brief Return the vector on which is acting
 *
 * It return the vector used in getVExpr, to get this object
 *
 * \return the vector
 *
 */
    vtype & getVector()
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
};

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
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(v.template getProp<prp>(k))
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
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(getExprL(v.getPos(k)))
	{
		return getExprL(v.getPos(k));
	}

#endif
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
	static inline auto value(vector & v, const vect_dist_key_dx & k) -> decltype(getExprR(v.getPos(k)))
	{
		return getExprR(v.getPos(k));
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

	typedef typename exp1::vtype vtype;

	//! constructor from an expresssion
	vector_dist_expression_op(const exp1 & o1)
	:o1(o1)
	{}

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
	template<typename r_type=typename std::remove_reference<decltype(-(o1.value(vect_dist_key_dx(0))))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return -(o1.value(key));
	}

    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
	    coeff_type coeff_tmp = -coeff;
        o1.template value_nz<Sys_eqs>(p_map,key,cols,coeff_tmp, comp);
    }
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
	vector & v;

public:

	//! The type of the internal vector
	typedef vector vtype;

	//! Property id of the point
	static const unsigned int prop = prp;

	int var_id = 0;

	void setVarId(int var_id)
	{
		this->var_id = var_id;
	}

	//! constructor for an external vector
	vector_dist_expression(vector & v)
	:v(v)
	{}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	vector & getVector()
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
    const vector & getVector() const
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
	inline auto value(const vect_dist_key_dx & k) const -> decltype(pos_or_propR<vector,prp>::value(v,k))
	{
		return pos_or_propR<vector,prp>::value(v,k);
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
		v_exp.init();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			pos_or_propL<vector,prp>::value(v,key) = v_exp.value(key);

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
	template<typename exp1, typename exp2, unsigned int op> vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
	{
		v_exp.init();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			pos_or_propL<vector,prp>::value(v,key) = v_exp.value(key);

			++it;
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
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			pos_or_propL<vector,prp>::value(v,key) = d;

			++it;
		}

		return v;
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

template<unsigned int, bool is_valid>
struct get_vector_dist_expression_op
{
	template<typename exp_type>
	inline static auto get(exp_type & o1, const vect_dist_key_dx & key) -> decltype(o1.value(vect_dist_key_dx(0)))
	{
		return o1.value(key);
	}

	template<unsigned int prop, typename exp_type, typename vector_type>
	inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key)
	{
		pos_or_propL<vector_type,exp_type::prop>::value(v,key) = o1.value(key);
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
	static int get(exp_type & o1, const vect_dist_key_dx & key, const int (& comp)[1])
	{
		return 0;
	}

	template<unsigned int prop, typename exp_type, typename vector_type>
	inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key)
	{
	}

	template<unsigned int prop, typename vector_type>
	inline static void assign_double(double d, vector_type & v, const vect_dist_key_dx & key)
	{
	}
};

template<>
struct get_vector_dist_expression_op<1,true>
{
	template<typename exp_type>
	static auto get(exp_type & o1, const vect_dist_key_dx & key, const int (& comp)[1]) -> decltype(o1.value(vect_dist_key_dx(0))[0])
	{
		return o1.value(key)[comp[0]];
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const int (& comp)[1])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]] = o1.value(key);
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
	static auto get(exp_type & o1, const vect_dist_key_dx & key, const int (& comp)[2]) -> decltype(o1.value(vect_dist_key_dx(0))[0][0])
	{
		return o1.value(key)[comp[0]][comp[1]];
	}

	template<unsigned int prop,typename exp_type, typename vector_type>
	inline static void assign(exp_type & o1, vector_type & v, const vect_dist_key_dx & key, const int (& comp)[2])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = o1.value(key);
	}

	template<unsigned int prop, typename vector_type>
	inline static void assign_double(double d, vector_type & v, const vect_dist_key_dx & key, const int (& comp)[2])
	{
		pos_or_propL<vector_type,prop>::value(v,key)[comp[0]][comp[1]] = d;
	}
};

/*! \brief like std::rank but it also work for openfpm structures like Point where it return 1
 *
 * \tparam T structure to check
 *
 */
template<typename T, bool is_point = is_Point<T>::value>
struct rank_gen
{
	typedef boost::mpl::int_<std::rank<T>::value> type;
};

template<typename T>
struct rank_gen<T,true>
{
	typedef boost::mpl::int_<1> type;
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
    const vtype & getVector() const
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
    vtype & getVector()
    {
        return o1.getVector();
    }

	//! initialize the expression tree
	inline void init() const
	{
		o1.init();
	}

	//! property on which this view is acting
	typedef typename boost::mpl::at<typename vtype::value_type::type,boost::mpl::int_<exp1::prop>>::type property_act;

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
	inline auto value(const vect_dist_key_dx & key) const -> decltype(get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,vect_dist_key_dx(0),comp))
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
	template<unsigned int prp2> vtype & operator=(const vector_dist_expression<prp2,vtype> & v_exp)
	{
		v_exp.init();

		auto & v = getVector();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<exp1::prop>(v_exp,v,key,comp);

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
		v_exp.init();

		auto & v = getVector();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			get_vector_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<exp1::prop>(v_exp,v,key,comp);

			++it;
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
template <unsigned int prp,typename vector> inline vector_dist_expression<prp,vector > getV(vector & v)
{
	vector_dist_expression<prp,vector > exp_v(v);

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
	inline double value(const vect_dist_key_dx & k) const
	{
		return d;
	}


    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
    {
        cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += coeff;
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

	//! type of object the structure return then evaluated
	typedef float vtype;

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
template<typename exp1 , typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_SUM>
operator+(const vector_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,vector_dist_expression<0,double>,VECT_SUM> exp_sum(va,vector_dist_expression<0,double>(d));

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
template<unsigned int prp1, typename v1>
inline vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_SUB>
operator-(double d, const vector_dist_expression<prp1,v1> & vb)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,vector_dist_expression<prp1,v1>,VECT_SUB> exp_sum(vector_dist_expression<0,double>(d),vb);

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
