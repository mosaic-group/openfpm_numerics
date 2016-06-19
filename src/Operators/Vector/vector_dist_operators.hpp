/*
 * vector_dist_operators.hpp
 *
 *  Created on: Jun 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_


#define VECT_SUM 1
#define VECT_SUB 2
#define VECT_MUL 3
#define VECT_DIV 4
#define VECT_NORM 5
#define VECT_NORM2 6
#define VECT_APPLYKER 7

/*! \brief is_expression check if a type is simple a type or is just an encapsulation of an expression
 *
 * return true if T::is_expression is a valid expression
 *
 */

template<typename ObjType, typename Sfinae = void>
struct is_expression: std::false_type {};

template<typename ObjType>
struct is_expression<ObjType, typename Void<typename ObjType::is_expression>::type> : std::true_type
{};

template<typename exp, bool is_exp = is_expression<exp>::value>
struct apply_kernel_rtype
{
	typedef typename exp::orig_type rtype;
};

template<typename exp>
struct apply_kernel_rtype<exp,false>
{
	typedef exp rtype;
};


/*! \brief Apply the kernel to particle differently that is a number or is an expression
 *
 *
 */
template<typename T, typename vector, typename exp,typename NN_type, typename Kernel, typename rtype, bool is_exp=is_expression<T>::value>
struct apply_kernel_is_number_or_expression
{
	inline static rtype apply(const vector & vd, NN_type & cl, const exp & v_exp, const vect_dist_key_dx & key, Kernel & lker)
	{
	    // accumulator
	    rtype pse = static_cast<typename rtype::coord_type>(0.0);

	    // position of particle p
	    Point<vector::dims,typename vector::stype> p = vd.getPos(key);

	    // property of the particle x
	    rtype prp_p = v_exp.value(key);

	    // Get the neighborhood of the particle
	    auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
	    while(NN.isNext())
	    {
	    	auto nnp = NN.get();

	    	// Calculate contribution given by the kernel value at position p,
	    	// given by the Near particle, exclude itself
	    	if (nnp != key.getKey())
	    	{
	    	    // property of the particle x
	    		rtype prp_q = v_exp.value(nnp);

	    	    // position of the particle q
	    		Point<vector::dims,typename vector::stype> q = vd.getPos(nnp);

	    	    for (size_t i = 0 ; i < T::nvals ; i++)
	    	    {
	    	    	float val = lker.value(p,q,prp_p.value(i),prp_q.value(i));

	    	    	// W(x-y)
	    	    	pse.value(i) += val;
	    	    }
	    	}

	    	// Next particle
	    	++NN;
	    }

	    return pse;
	}
};

template<typename T, typename vector, typename exp,typename NN_type, typename Kernel, typename rtype>
struct apply_kernel_is_number_or_expression<T,vector,exp,NN_type,Kernel,rtype,false>
{

	/*! \brief Apply the kernel
	 *
	 * \param cl Neighborhood of particles
	 * \param v_exp vector expression
	 * \param key particle id
	 * \param lker kernel function
	 *
	 */
	inline static rtype apply(const vector & vd, NN_type & cl, const exp & v_exp, const vect_dist_key_dx & key, const Kernel & lker)
	{
	    // accumulator
	    rtype pse = 0;

	    // Get f(x) at the position of the particle
	    rtype & prp_p = v_exp.value(key);

	    // position of particle p
	    auto & p = vd.getPos(key);

	    // Get the neighborhood of the particle
	    auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
	    while(NN.isNext())
	    {
	    	auto nnp = NN.get();

	    	// Calculate contribution given by the kernel value at position p,
	    	// given by the Near particle, exclude itself
	    	if (nnp != key.getKey())
	    	{
	    	    // Get f(x) at the position of the particle
	    	    rtype & prp_q = vd.getPos(nnp);

	    	    // position of particle q
	    	    auto & q = v_exp.value_pos(key);

	    		// W(x-y)
	    		auto ker = lker.value(p,q,prp_p,prp_q);

	    		pse += ker;
	    	}

	    	// Next particle
	    	++NN;
	    }

	    return pse;
	}
};

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

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER>
{
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	const exp1 o1;
	NN & cl;
	Kernel & ker;
	const vector_orig & vd;

public:

	vector_dist_expression_op(const exp1 & o1, NN & cl, Kernel & ker, const vector_orig & vd)
	:o1(o1),cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 */
	template<typename r_type=typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx(0)))>::rtype > inline r_type value(const vect_dist_key_dx & key) const
	{
		typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx(0)))>::rtype rtype;

		return apply_kernel_is_number_or_expression<decltype(o1.value(key)),vector_orig,exp1,NN,Kernel,rtype>::apply(vd,cl,o1,key,ker);
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
		auto it = v.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			auto exp = v_exp.value(key);
			call_init_if_needed<decltype(exp)>::call(exp);
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


///////////////////////////////////////////// Apply kernel operator ////
////////////////////////////////////////////////////////////////////////

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,boost::mpl::vector<NN,Kernel,vector_type>,VECT_APPLYKER>
applyKernel(const vector_dist_expression_op<exp1,exp2,op1> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,boost::mpl::vector<NN,Kernel,vector_type>,VECT_APPLYKER> exp_sum(va,cl,ker,vd);

	return exp_sum;
}



#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_HPP_ */
