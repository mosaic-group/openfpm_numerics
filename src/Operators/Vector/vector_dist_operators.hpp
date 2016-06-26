/*
 * vector_dist_operators.hpp
 *
 *  Created on: Jun 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_

#include "Vector/vector_dist.hpp"

#define VECT_SUM 1
#define VECT_SUB 2
#define VECT_MUL 3
#define VECT_DIV 4
#define VECT_NORM 5
#define VECT_NORM2 6
#define VECT_APPLYKER_IN 7
#define VECT_APPLYKER_OUT 8
#define VECT_APPLYKER_REDUCE 9
#define VECT_APPLYKER_IN_GEN 10
#define VECT_APPLYKER_OUT_GEN 11
#define VECT_APPLYKER_REDUCE_GEN 12
#define VECT_APPLYKER_MULTI_IN 13
#define VECT_APPLYKER_MULTI_OUT 14
#define VECT_APPLYKER_MULTI_REDUCE 15
#define VECT_APPLYKER_MULTI_IN_GEN 16
#define VECT_APPLYKER_MULTI_OUT_GEN 17
#define VECT_APPLYKER_MULTI_REDUCE_GEN 18


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

/*! \brief Call the init function if needed
 *
 * \param r_exp expression
 *
 */
template <typename T, bool has_init = has_init<T>::value >
struct call_init_if_needed
{
	static inline void call(T & r_exp)
	{
		r_exp.init();
	}
};

template <typename T>
struct call_init_if_needed<T,false>
{
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
	const exp1 o1;
	const exp2 o2;

public:

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
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)) + o2.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
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
	const exp1 o1;
	const exp2 o2;

public:

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
	const exp1 o1;
	const exp2 o2;

public:

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
	const exp1 o1;
	const exp2 o2;

public:

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
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)) / o2.value(vect_dist_key_dx(0)))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return o1.value(key) / o2.value(key);
	}
};

/*! \brief norm operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1>
class vector_dist_expression_op<exp1,void,VECT_NORM>
{
	const exp1 o1;

public:

	vector_dist_expression_op(const exp1 & o1)
	:o1(o1)
	{}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(norm(o1.value(vect_dist_key_dx(0))))>::type > inline r_type value(const vect_dist_key_dx & key) const
	{
		return norm(o1.value(key));
	}
};


/*! \brief Main class that encapsulate a vector properties
 *
 * \tparam prp property involved
 * \tparam vector involved
 *
 */
template<unsigned int prp, typename vector>
class vector_dist_expression
{
	vector & v;

public:

	vector_dist_expression(vector & v)
	:v(v)
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
	inline auto value(const vect_dist_key_dx & k) const -> decltype(v.template getProp<prp>(k))
	{
		return v.template getProp<prp>(k);
	}


	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param v_exp expression to evaluate
	 *
	 */
	template<typename exp1, typename exp2, unsigned int op> vector & operator=(const vector_dist_expression_op<exp1,exp2,op> & v_exp)
	{
		v_exp.init();

		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			v.template getProp<prp>(key) = v_exp.value(key);

			++it;
		}

		return v;
	}

	/*! \brief Fill the vector property with the double
	 *
	 * \param d value to fill
	 *
	 */
	vector & operator=(double d)
	{
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			v.template getProp<prp>(key) = d;

			++it;
		}

		return v;
	}
};

/*! \Create an expression from a vector
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
	double d;

public:

	inline vector_dist_expression(double & d)
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
	 * It just return the velue set in the constructor
	 *
	 */
	inline double value(const vect_dist_key_dx & k) const
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
/////////// NORM operator ///////////////////////

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,void,VECT_NORM>
norm(const vector_dist_expression_op<exp1,exp2,op1> & va)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,void,VECT_NORM> exp_sum(va);

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
inline vector_dist_expression_op<vector_dist_expression<0,double>,void,VECT_NORM>
norm(double d)
{
	vector_dist_expression_op<vector_dist_expression<0,double>,void,VECT_NORM> exp_sum(vector_dist_expression<0,double>(d));

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
inline vector_dist_expression_op<vector_dist_expression<prp1,v1>,void,VECT_NORM>
norm(const vector_dist_expression<prp1,v1> & va)
{
	vector_dist_expression_op<vector_dist_expression<prp1,v1>,void,VECT_NORM> exp_sum(va);

	return exp_sum;
}

#include "vector_dist_operators_apply_kernel.hpp"

#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_ */
