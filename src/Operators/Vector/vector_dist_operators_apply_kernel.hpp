/*
 * vector_dist_operators_apply_kernel.hpp
 *
 *  Created on: Jun 19, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_APPLY_KERNEL_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_APPLY_KERNEL_HPP_

//////// SET of small macro to make to define integral easy

#define DEFINE_INTERACTION_3D(name) struct name \
{\
\
	Point<3,double> value(const Point<3,double> & xp, const Point<3,double> xq)\
	{

#define END_INTERACTION }\
						};

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

//! It give the return type of the expression if applicable
template<typename exp, bool is_exp = is_expression<exp>::value>
struct apply_kernel_rtype
{
	//! indicate the return type of the expression exp
	typedef typename exp::return_type rtype;
};

//! It give the return type of the expression if applicable
template<typename exp>
struct apply_kernel_rtype<exp,false>
{
	//! indicate the return type of the expression exp
	typedef exp rtype;
};


/*! \brief Meta-function to return a compatible zero-element
 *
 *
 */
template<typename rtype>
struct set_zero
{
	//! return 0
	__device__ __host__ static rtype create()
	{
		return 0.0;
	}
};

//! Create a point with all compunent set to zero
template<unsigned int dim, typename T>
struct set_zero<Point<dim,T>>
{
	//! return a point with all the coordinates set to zero
	__device__ __host__ static Point<dim,T> create()
	{
		Point<dim,T> ret;

		ret.zero();
		return ret;
	}
};

constexpr int NN_index_sort = 1;
constexpr int NN_index_unsort = 0;

template<unsigned int impl>
struct NN_index
{
	template<typename NN_type>
	__device__ __host__ static auto get(NN_type & NN) -> decltype(NN.get())
	{
		return NN.get();
	}
};

template<>
struct NN_index<NN_index_sort>
{
	template<typename NN_type>
	__device__ __host__ static auto get(NN_type & NN) -> decltype(NN.get_sort())
	{
		return NN.get_sort();
	}
};

/*! \brief Apply the kernel to particle differently that is a number or is an expression
 *
 *
 */
template<unsigned int impl, typename T, typename vector, typename exp,typename NN_type, typename Kernel, typename rtype, bool is_exp=is_expression<T>::value>
struct apply_kernel_is_number_or_expression
{
	/*! \brief Apply the kernel expression to a particle
	 *
	 * \param vd vector with the particles positions
	 * \param cl Cell-list
	 * \param v_exp expression to execute for each particle
	 * \param key to which particle to apply the expression
	 * \param lker kernel function
	 *
	 * \return the result of apply the kernel on the particle key
	 *
	 */
	inline __host__ __device__ static typename std::remove_reference<rtype>::type apply(const vector & vd, NN_type & cl, const exp & v_exp, const vect_dist_key_dx & key, Kernel & lker)
	{
	    // accumulator
		typename std::remove_reference<rtype>::type pse = set_zero<typename std::remove_reference<rtype>::type>::create();

	    // position of particle p
	    Point<vector::dims,typename vector::stype> p = vd.getPos(key);

	    // property of the particle x
	    rtype prp_p = v_exp.value(key);

	    // Get the neighborhood of the particle
	    auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
	    while(NN.isNext())
	    {
	    	auto nnp = NN_index<impl>::get(NN);

	    	// Calculate contribution given by the kernel value at position p,
	    	// given by the Near particle, exclude itself
	    	if (nnp != key.getKey())
	    	{
	    		vect_dist_key_dx nnp_k;
	    		nnp_k.setKey(nnp);

	    	    // property of the particle x
	    		rtype prp_q = v_exp.value(nnp_k);

	    	    // position of the particle q
	    		Point<vector::dims,typename vector::stype> q = vd.getPos(nnp);

	    	    pse += lker.value(p,q,prp_p,prp_q);
	    	}

	    	// Next particle
	    	++NN;
	    }

	    return pse;
	}
};


/*! \brief Apply the kernel to particle differently that is a number or is an expression
 *
 *
 */
template<unsigned int impl, typename vector, typename exp,typename NN_type, typename Kernel, typename rtype>
struct apply_kernel_is_number_or_expression_sim
{
	/*! \brief Apply the kernel expression to a particle
	 *
	 * \param vd vector with the particles positions
	 * \param cl Cell-list (symmetric version)
	 * \param key to which particle to apply the expression
	 * \param lker kernel function
	 *
	 * \return the result of apply the kernel on the particle key
	 *
	 */
	__device__ __host__ inline static typename std::remove_reference<rtype>::type apply(const vector & vd, NN_type & cl, const vect_dist_key_dx & key, Kernel & lker)
	{
	    // accumulator
		typename std::remove_reference<rtype>::type pse = set_zero<typename std::remove_reference<rtype>::type>::create();

	    // position of particle p
	    Point<vector::dims,typename vector::stype> p = vd.getPos(key);

	    // Get the neighborhood of the particle
	    auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
	    while(NN.isNext())
	    {
	    	auto nnp = NN_index<impl>::get(NN);

	    	// Calculate contribution given by the kernel value at position p,
	    	// given by the Near particle, exclude itself
	    	if (nnp != key.getKey())
	    	{
	    	    // position of the particle q
	    		Point<vector::dims,typename vector::stype> q = vd.getPos(nnp);

	    	    pse += lker.value(p,q);
	    	}

	    	// Next particle
	    	++NN;
	    }

	    return pse;
	}
};



/*! \brief Apply the kernel to particle differently that is a number or is an expression
 *
 *
 */
template<unsigned int impl, typename T, typename vector, typename exp,typename NN_type, typename Kernel, typename rtype, bool is_exp=is_expression<T>::value>
struct apply_kernel_is_number_or_expression_gen
{
	/*! \brief Apply the kernel expression to a particle
	 *
	 * \param vd vector with the particles positions
	 * \param cl Cell-list
	 * \param v_exp expression to execute for each particle
	 * \param key to which particle to apply the expression
	 * \param lker kernel function
	 *
	 * \return the result of apply the kernel on the particle key
	 *
	 */
	__device__ __host__ inline static typename std::remove_reference<rtype>::type apply(const vector & vd, NN_type & cl, const exp & v_exp, const vect_dist_key_dx & key, Kernel & lker)
	{
	    // accumulator
	    typename std::remove_reference<rtype>::type pse = set_zero<typename std::remove_reference<rtype>::type>::create();

	    // property of the particle x
	    rtype prp_p = v_exp.value(key);

	    // position of particle p
	    Point<vector::dims,typename vector::stype> p = vd.getPos(key);

	    // Get the neighborhood of the particle
	    auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
	    while(NN.isNext())
	    {
	    	auto nnp = NN_index<impl>::get(NN);

	    	// Calculate contribution given by the kernel value at position p,
	    	// given by the Near particle, exclude itself
	    	if (nnp != key.getKey())
	    	{
	    		vect_dist_key_dx nnp_k;
	    		nnp_k.setKey(nnp);

	    	    // property of the particle x
	    		rtype prp_q = v_exp.value(nnp_k);

	    	    pse += lker.value(key.getKey(),nnp,prp_p,prp_q,vd);
	    	}

	    	// Next particle
	    	++NN;
	    }

	    return pse;
	}
};

template<typename T, bool mm>
struct mutable_or_not
{
	T obj;

	mutable_or_not(T & obj)
	:obj(obj)
	{}
};

template<typename T>
struct mutable_or_not<T,true>
{
	mutable T obj;

	mutable_or_not(T & obj)
	:obj(obj)
	{}
};

template<typename T, bool is_reference = std::is_reference<T>::value>
struct add_const_reference
{
	typedef const T type;
};

template<typename T>
struct add_const_reference<T,true>
{
	typedef typename std::remove_reference<T>::type const & type;
};

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER_IN>
{
	//! Get the type of the Cell-list
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	//! Get the type of the kernel
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;
	//! Get the type that contain the particle positions
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	//! Return the type of the Cell-list
	typedef typename std::remove_reference<NN>::type NN_nr;
	//! Return the type of the kernel
	typedef typename std::remove_reference<Kernel>::type Kernel_nr;
	//! Return the vector containing the position of the particles
	typedef typename std::remove_reference<vector_orig>::type vector_orig_nr;

	//! expression
	const exp1 o1;

	//! Cell-list
	mutable_or_not<NN,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> cl;

	//! kernel
	mutable_or_not<Kernel,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> ker;

	//! The vector that contain the particles
	typename add_const_reference<vector_orig>::type vd;

	//! Get the return type of applying the kernel to a particle
	typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx()))>::rtype rtype;

public:

	//! indicate if this expression operate on kernel expression
	typedef typename has_vector_kernel<vector_orig_nr>::type is_ker;

	//! vector type on which this expression work
	typedef vector_orig_nr vtype;

	//! indicate that is apply_kernel is not a sorted version
	typedef boost::mpl::bool_<false> is_sort;

	//! type of NN structure
	typedef NN_nr NN_type;

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_nr * getNN() const
	{
		return &cl.obj;
	}

	/*! \brief get the vector
	 *
	 * \return the vector
	 *
	 */
	inline const vtype & getVector() const
	{
		return vd;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
	}

	/*! \brief Constructor
	 *
	 * \param o1 expression
	 * \param cl Cell-list
	 * \param ker kernel function
	 * \param vd vector containing the particle positions
	 *
	 */
	vector_dist_expression_op(const exp1 & o1, NN_nr & cl, Kernel & ker, const vector_orig_nr & vd)
	:o1(o1),cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__device__ __host__ inline typename std::remove_reference<rtype>::type value(const vect_dist_key_dx & key) const
	{
		return apply_kernel_is_number_or_expression<NN_index_unsort,
													decltype(o1.value(key)),
													vector_orig_nr,
													exp1,
													NN_nr,
													Kernel_nr,
													rtype>::apply(vd,cl.obj,o1,key,ker.obj);
	}
};

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER_IN_SORT>
{
	//! Get the type of the Cell-list
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	//! Get the type of the kernel
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;
	//! Get the type that contain the particle positions
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	//! Return the type of the Cell-list
	typedef typename std::remove_reference<NN>::type NN_nr;
	//! Return the type of the kernel
	typedef typename std::remove_reference<Kernel>::type Kernel_nr;
	//! Return the vector containing the position of the particles
	typedef typename std::remove_reference<vector_orig>::type vector_orig_nr;

	//! expression
	const exp1 o1;

	//! Cell-list
	mutable_or_not<NN,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> cl;

	//! kernel
	mutable_or_not<Kernel,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> ker;

	//! The vector that contain the particles
	typename add_const_reference<vector_orig>::type vd;

	//! Get the return type of applying the kernel to a particle
	typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx()))>::rtype rtype;

public:

	//! indicate if this expression operate on kernel expression
	typedef typename has_vector_kernel<vector_orig_nr>::type is_ker;

	//! vector type on which this expression work
	typedef vector_orig_nr vtype;

	//! indicate that is apply_kernel is not a sorted version
	typedef boost::mpl::bool_<true> is_sort;

	//! type of NN structure
	typedef NN_nr NN_type;

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_nr * getNN() const
	{
		return &cl.obj;
	}

	/*! \brief get the vector
	 *
	 * \return the vector
	 *
	 */
	inline const vtype & getVector() const
	{
		return vd;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
	}

	/*! \brief Constructor
	 *
	 * \param o1 expression
	 * \param cl Cell-list
	 * \param ker kernel function
	 * \param vd vector containing the particle positions
	 *
	 */
	vector_dist_expression_op(const exp1 & o1, NN_nr & cl, Kernel & ker, const vector_orig_nr & vd)
	:o1(o1),cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__device__ __host__ inline typename std::remove_reference<rtype>::type value(const vect_dist_key_dx & key) const
	{
		return apply_kernel_is_number_or_expression<NN_index_sort,
													decltype(o1.value(key)),
													vector_orig_nr,
													exp1,
													NN_nr,
													Kernel_nr,
													rtype>::apply(vd,cl.obj,o1,key,ker.obj);
	}
};

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER_IN_SIM>
{
	//! Return the type of the Cell-list
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	//! Return the type of the kernel
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;
	//! Return the vector containing the position of the particles
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	//! Return the type of the Cell-list
	typedef typename std::remove_reference<NN>::type NN_nr;
	//! Return the type of the kernel
	typedef typename std::remove_reference<Kernel>::type Kernel_nr;
	//! Return the vector containing the position of the particles
	typedef typename std::remove_reference<vector_orig>::type vector_orig_nr;

	//! Cell-list
	mutable_or_not<NN,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> cl;

	//! kernel
	mutable_or_not<Kernel,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> ker;

	//! vector with the particle positions
	const vector_orig vd;

	//! Get the return type of the expression
	typedef typename apply_kernel_rtype<decltype(std::declval<Kernel_nr>().value(Point<vector_orig_nr::dims,typename vector_orig_nr::stype>(0.0), Point<vector_orig_nr::dims,typename vector_orig_nr::stype>(0.0) ) )>::rtype rtype;


public:

	//! indicate if this expression operate on kernel expression
	typedef typename has_vector_kernel<vector_orig_nr>::type is_ker;

	//! vector type on which this expression work
	typedef vector_orig_nr vtype;

	//! indicate that is apply_kernel is not a sorted version
	typedef boost::mpl::bool_<false> is_sort;

	//! type of NN structure
	typedef NN_nr NN_type;

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_nr * getNN() const
	{
		return &cl.obj;
	}

	/*! \brief get the vector
	 *
	 * \return the vector
	 *
	 */
	inline const vtype & getVector() const
	{
		return vd;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
	}

	/*! \brief Constructor
	 *
	 * \param cl cell-list
	 * \param ker kernel
	 * \param vd Vector containing the particle positions
	 *
	 */
	vector_dist_expression_op(NN & cl, Kernel & ker, const vector_orig & vd)
	:cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result produced by the expression
	 *
	 */
	__device__ __host__ inline typename std::remove_reference<rtype>::type value(const vect_dist_key_dx & key) const
	{
		return apply_kernel_is_number_or_expression_sim<NN_index_unsort,
														vector_orig_nr,
														exp1,
														NN_nr,
														Kernel_nr,
														rtype>::apply(vd,cl.obj,key,ker.obj);
	}
};

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER_IN_SIM_SORT>
{
	//! Return the type of the Cell-list
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	//! Return the type of the kernel
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;
	//! Return the vector containing the position of the particles
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	//! Return the type of the Cell-list
	typedef typename std::remove_reference<NN>::type NN_nr;
	//! Return the type of the kernel
	typedef typename std::remove_reference<Kernel>::type Kernel_nr;
	//! Return the vector containing the position of the particles
	typedef typename std::remove_reference<vector_orig>::type vector_orig_nr;

	//! Cell-list
	mutable_or_not<NN,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> cl;

	//! kernel
	mutable_or_not<Kernel,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> ker;

	//! vector with the particle positions
	const vector_orig vd;

	//! Get the return type of the expression
	typedef typename apply_kernel_rtype<decltype(std::declval<Kernel_nr>().value(Point<vector_orig_nr::dims,typename vector_orig_nr::stype>(0.0), Point<vector_orig_nr::dims,typename vector_orig_nr::stype>(0.0) ) )>::rtype rtype;


public:

	//! indicate if this expression operate on kernel expression
	typedef typename has_vector_kernel<vector_orig_nr>::type is_ker;

	//! vector type on which this expression work
	typedef vector_orig_nr vtype;

	//! indicate that is apply_kernel is not a sorted version
	typedef boost::mpl::bool_<true> is_sort;

	//! type of NN structure
	typedef NN_nr NN_type;

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_nr * getNN() const
	{
		return &cl.obj;
	}

	/*! \brief get the vector
	 *
	 * \return the vector
	 *
	 */
	inline const vtype & getVector() const
	{
		return vd;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
	}

	/*! \brief Constructor
	 *
	 * \param cl cell-list
	 * \param ker kernel
	 * \param vd Vector containing the particle positions
	 *
	 */
	vector_dist_expression_op(NN & cl, Kernel & ker, const vector_orig & vd)
	:cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result produced by the expression
	 *
	 */
	__device__ __host__ inline typename std::remove_reference<rtype>::type value(const vect_dist_key_dx & key) const
	{
		return apply_kernel_is_number_or_expression_sim<NN_index_sort,
														vector_orig_nr,
														exp1,
														NN_nr,
														Kernel_nr,
														rtype>::apply(vd,cl.obj,key,ker.obj);
	}
};

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER_IN_GEN>
{
	//! Get the type of the Cell-list
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	//! Get the type of the kernel
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;

	//! Get the type of the vector containing the set of particles
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	//! Return the type of the Cell-list
	typedef typename std::remove_reference<NN>::type NN_nr;
	//! Return the type of the kernel
	typedef typename std::remove_reference<Kernel>::type Kernel_nr;
	//! Return the vector containing the position of the particles
	typedef typename std::remove_reference<vector_orig>::type vector_orig_nr;

	//! Expression
	const exp1 o1;

	//! Cell-list
	mutable_or_not<NN,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> cl;

	//! kernel
	mutable_or_not<Kernel,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> ker;

	//! Vector containing the particles
	const vector_orig vd;

	//! Return type of the expression
	typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx()))>::rtype rtype;

public:

	//! indicate if this expression operate on kernel expression
	typedef typename has_vector_kernel<vector_orig_nr>::type is_ker;

	//! vector type on which this expression work
	typedef vector_orig_nr vtype;

	//! indicate that is apply_kernel is not a sorted version
	typedef boost::mpl::bool_<false> is_sort;

	//! type of NN structure
	typedef NN_nr NN_type;

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_nr * getNN() const
	{
		return &cl.obj;
	}

	/*! \brief get the vector
	 *
	 * \return the vector
	 *
	 */
	inline const vtype & getVector() const
	{
		return vd;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
	}

	/*! \brief Constructor
	 *
	 * \param o1 expression
	 * \param cl Cell-list
	 * \param ker Kernel
	 * \param vd vector containing the set of particles
	 *
	 */
	vector_dist_expression_op(const exp1 & o1, NN & cl, Kernel & ker, const vector_orig & vd)
	:o1(o1),cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__device__ __host__ inline typename std::remove_reference<rtype>::type value(const vect_dist_key_dx & key) const
	{
		return apply_kernel_is_number_or_expression_gen<NN_index_unsort,
														decltype(o1.value(key)),
														vector_orig_nr,
														exp1,
														NN_nr,
														Kernel_nr,
														rtype>::apply(vd,cl.obj,o1,key,ker.obj);
	}
};

/*! \brief Apply kernel operation
 *
 * \tparam exp1 expression1
 * \tparam NN list
 *
 */
template <typename exp1,typename vector_type>
class vector_dist_expression_op<exp1,vector_type,VECT_APPLYKER_IN_GEN_SORT>
{
	//! Get the type of the Cell-list
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<0>>::type NN;
	//! Get the type of the kernel
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<1>>::type Kernel;

	//! Get the type of the vector containing the set of particles
	typedef typename boost::mpl::at<vector_type,boost::mpl::int_<2>>::type vector_orig;

	//! Return the type of the Cell-list
	typedef typename std::remove_reference<NN>::type NN_nr;
	//! Return the type of the kernel
	typedef typename std::remove_reference<Kernel>::type Kernel_nr;
	//! Return the vector containing the position of the particles
	typedef typename std::remove_reference<vector_orig>::type vector_orig_nr;

	//! Expression
	const exp1 o1;

	//! Cell-list
	mutable_or_not<NN,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> cl;

	//! kernel
	mutable_or_not<Kernel,is_gpu_celllist<NN>::type::value || is_gpu_ker_celllist<NN>::type::value> ker;

	//! Vector containing the particles
	const vector_orig vd;

	//! Return type of the expression
	typedef typename apply_kernel_rtype<decltype(o1.value(vect_dist_key_dx()))>::rtype rtype;

public:

	//! indicate if this expression operate on kernel expression
	typedef typename has_vector_kernel<vector_orig_nr>::type is_ker;

	//! vector type on which this expression work
	typedef vector_orig_nr vtype;

	//! indicate that is apply_kernel is not a sorted version
	typedef boost::mpl::bool_<true> is_sort;

	//! type of NN structure
	typedef NN_nr NN_type;

	/*! \brief get the NN object
	 *
	 * \return the NN object
	 *
	 */
	inline NN_nr * getNN() const
	{
		return &cl.obj;
	}

	/*! \brief get the vector
	 *
	 * \return the vector
	 *
	 */
	inline const vtype & getVector() const
	{
		return vd;
	}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const
	{
		o1.init();
	}

	/*! \brief Constructor
	 *
	 * \param o1 expression
	 * \param cl Cell-list
	 * \param ker Kernel
	 * \param vd vector containing the set of particles
	 *
	 */
	vector_dist_expression_op(const exp1 & o1, NN & cl, Kernel & ker, const vector_orig & vd)
	:o1(o1),cl(cl),ker(ker),vd(vd)
	{}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	__device__ __host__ inline typename std::remove_reference<rtype>::type value(const vect_dist_key_dx & key) const
	{
		return apply_kernel_is_number_or_expression_gen<NN_index_sort,
														decltype(o1.value(key)),
														vector_orig_nr,
														exp1,
														NN_nr,
														Kernel_nr,
														rtype>::apply(vd,cl.obj,o1,key,ker.obj);
	}
};

template<typename cl_type, int impl=2*is_gpu_celllist<cl_type>::type::value + is_gpu_ker_celllist<cl_type>::type::value >
struct cl_selector_impl
{
	typedef cl_type& ctype;

	static ctype & get(cl_type & cl)
	{
		return cl;
	}
};

template<typename cl_type>
struct cl_selector_impl<cl_type,2>
{
	typedef decltype(std::declval<cl_type>().toKernel()) ctype;

	static ctype get(cl_type & cl)
	{
		return cl.toKernel();
	}
};

template<typename cl_type>
struct cl_selector_impl<cl_type,1>
{
	typedef cl_type ctype;

	static ctype & get(cl_type & cl)
	{
		return cl;
	}
};

template<typename kl_type, typename cl_type, int impl=2*is_gpu_celllist<cl_type>::type::value + is_gpu_ker_celllist<cl_type>::type::value >
struct kl_selector_impl
{
	typedef kl_type& ktype;
};

template<typename kl_type, typename cl_type>
struct kl_selector_impl<kl_type,cl_type,2>
{
	typedef kl_type ktype;
};

template<typename kl_type, typename cl_type>
struct kl_selector_impl<kl_type,cl_type,1>
{
	typedef kl_type ktype;
};

template<bool is_sorted, bool uwk, typename T>
struct unroll_with_to_kernel
{
	typedef T type;

	static T & get(T & v)
	{
		return v;
	}
};

template<typename T>
struct unroll_with_to_kernel<false,true,T>
{
	typedef decltype(std::declval<T>().toKernel()) type;

	static type get(T & v)
	{
		return v.toKernel();
	}
};

template<typename T>
struct unroll_with_to_kernel<true,true,T>
{
	typedef decltype(std::declval<T>().toKernel_sorted()) type;

	static type get(T & v)
	{
		return v.toKernel_sorted();
	}
};

template<bool is_sorted, typename vl_type, typename cl_type, int impl=2*is_gpu_celllist<cl_type>::type::value + is_gpu_ker_celllist<cl_type>::type::value>
struct vl_selector_impl
{
	typedef vl_type& vtype;

	static vtype & get(vl_type & v)
	{
		return v;
	}
};

template<typename vl_type, typename cl_type>
struct vl_selector_impl<false,vl_type,cl_type,2>
{
	typedef decltype(std::declval<vl_type>().toKernel()) vtype;

	static vtype & get(vl_type & v)
	{
		return v.toKernel();
	}
};

template<typename vl_type, typename cl_type>
struct vl_selector_impl<true,vl_type,cl_type,2>
{
	typedef decltype(std::declval<vl_type>().toKernel()) vtype;

	static vtype & get(vl_type & v)
	{
		return v.toKernel_sorted();
	}
};

template<bool is_sorted, typename vl_type, typename cl_type>
struct vl_selector_impl<is_sorted, vl_type,cl_type,1>
{
	typedef typename unroll_with_to_kernel<is_sorted,has_toKernel<vl_type>::type::value,vl_type>::type vtype;

	static vtype get(vl_type & v)
	{
		return unroll_with_to_kernel<is_sorted,has_toKernel<vl_type>::type::value,vl_type>::get(v);
	}
};


///////////////////////////////////////////// Apply kernel operator ////
////////////////////////////////////////////////////////////////////////

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<false,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN>
applyKernel_in(const vector_dist_expression_op<exp1,exp2,op1> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<false,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN> exp_sum(va,
									  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	ker,
									  	  	  	  	  	vl_selector_impl<false,vector_type,NN>::get(vd));

	return exp_sum;
}


///////////////////////////////////////////// Apply kernel operator ////
////////////////////////////////////////////////////////////////////////

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<true,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_SORT>
applyKernel_in_sort(const vector_dist_expression_op<exp1,exp2,op1> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<true,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_SORT> exp_sum(va,
									  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	ker,
									  	  	  	  	  	vl_selector_impl<true,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<false,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_GEN>
applyKernel_in_gen(const vector_dist_expression_op<exp1,exp2,op1> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<false,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_GEN> exp_sum(va,
									  	  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	  	ker,
									  	  	  	  	  	  	vl_selector_impl<false,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, unsigned int op1, typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<true,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_GEN_SORT>
applyKernel_in_gen_sort(const vector_dist_expression_op<exp1,exp2,op1> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression_op<exp1,exp2,op1>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<true,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_GEN_SORT> exp_sum(va,
									  	  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	  	ker,
									  	  	  	  	  	  	vl_selector_impl<true,vector_type,NN>::get(vd));

	return exp_sum;
}

//////////////////////////////////////// For vector expression ///////////////////////

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp, typename NN, typename Kernel, typename vector_type, typename vector_type_ker>
inline vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<false,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN>
applyKernel_in(const vector_dist_expression<prp,vector_type_ker> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<false,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN> exp_sum(va,
									  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	ker,
									  	  	  	  	  	vl_selector_impl<false,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp, typename NN, typename Kernel, typename vector_type, typename vector_type_ker>
inline vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<true,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_SORT>
applyKernel_in_sort(const vector_dist_expression<prp,vector_type_ker> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<true,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_SORT> exp_sum(va,
									  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	ker,
									  	  	  	  	  	vl_selector_impl<true,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp, typename NN, typename Kernel, typename vector_type, typename vector_type_ker>
inline vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<false,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_GEN>
applyKernel_in_gen(const vector_dist_expression<prp,vector_type_ker> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<false,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_GEN> exp_sum(va,
									  	  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	  	ker,
									  	  	  	  	  	  	vl_selector_impl<false,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp, typename NN, typename Kernel, typename vector_type, typename vector_type_ker>
inline vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<true,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_GEN_SORT>
applyKernel_in_gen_sort(const vector_dist_expression<prp,vector_type_ker> & va, vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<vector_dist_expression<prp,vector_type_ker>,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<true,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_GEN_SORT> exp_sum(va,
									  	  	  	  	  	  	cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	  	ker,
									  	  	  	  	  	  	vl_selector_impl<true,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<void,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<false,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_SIM>
applyKernel_in_sim(vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<void,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<false,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_SIM> exp_sum(cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	  	ker,
									  	  	  	  	  	  	vl_selector_impl<false,vector_type,NN>::get(vd));

	return exp_sum;
}

/* \brief Apply kernel expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename NN, typename Kernel, typename vector_type>
inline vector_dist_expression_op<void,
								 boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
								 	 	 	 	 	typename kl_selector_impl<Kernel,NN>::ktype,
								 	 	 	 	 	typename vl_selector_impl<true,vector_type,NN>::vtype>,
								 VECT_APPLYKER_IN_SIM_SORT>
applyKernel_in_sim_sort(vector_type & vd, NN & cl, Kernel & ker)
{
	vector_dist_expression_op<void,
							  boost::mpl::vector<typename cl_selector_impl<NN>::ctype,
							  	  	  	  	  	 typename kl_selector_impl<Kernel,NN>::ktype,
							  	  	  	  	  	 typename vl_selector_impl<true,vector_type,NN>::vtype>,
							  VECT_APPLYKER_IN_SIM_SORT> exp_sum(cl_selector_impl<NN>::get(cl),
									  	  	  	  	  	  	ker,
									  	  	  	  	  	  	vl_selector_impl<true,vector_type,NN>::get(vd));

	return exp_sum;
}

#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_APPLY_KERNEL_HPP_ */
