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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)) + o2.value(vect_dist_key_dx(0)))>::type >
	inline r_type value(const vect_dist_key_dx & key) const
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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)) - o2.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) - o2.value(key);
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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)) * o2.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
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
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)) / o2.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) / o2.value(key);
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

	//! constructor from an expresssion
	vector_dist_expression_op(const exp1 & o1)
	:o1(o1)
	{}

	//! initialize the expression tree
	inline void init() const
	{
		o1.init();
	}

	//! return the result of the expression
	template<typename r_type=typename std::remove_reference<decltype(-(o1.value(vect_dist_key_dx(0))))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return -(o1.value(key));
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
