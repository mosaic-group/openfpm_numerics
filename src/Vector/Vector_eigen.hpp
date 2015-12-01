/*
 * Vector_eigen.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_

/*!	\brief Copy scalar element
 *
 *
 *
 */
template<typename copy_type, typename T, typename Eqs_sys, bool sa>
struct copy_ele_sca_array
{
	template<typename Grid> inline void copy(Grid & grid_dst,size_t lin_id, size_t base_id, size_t gs_size)
	{
		if (Eqs_sys::ord == EQS_FIELD)
			grid_dst.template get<T::value>(key) = x(lin_id * Eqs_sys::nvar + base_id);
		else
			grid_dst.template get<T::value>(key) = x(base_id * gs_size + lin_id);
	}
};


/*! \brief Copy array element
 *
 *
 *
 */
template<typename copy_type, typename T, typename Eqs_sys>
struct copy_ele_sca_array<copy_type,T,Eqs_sys,true>
{
	template<typename Grid> inline void copy(Grid & grid_dst,size_t lin_id, size_t base_id, size_t gs_size)
	{
		for (size_t i = 0 ; i < std::extents<copy_type>::value ; i++)
		{
			if ()
				grid_dst.template get<T::value>(key)[i] = x(lin_id * Eqs_sys::nvar + base_id + i);
			else
				grid_dst.template get<T::value>(key)[i] = x(base_id * gs_size + lin_id);
		}
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy from the Vector_eigen to the grid target
 * in a generic way for a generic object T with variable number of properties
 *
 * \tparam dim Dimensionality
 * \tparam S type of grid
 * \tparam Ev Eigen vector type
 *
 */

template<unsigned int dim, typename S, typename Ev>
struct copy_ele
{
	//! where to copy the object in the grid
	const grid_dist_key_dx<dim> key;

	//! destination grid
	S & grid_dst;

	//! source Eigen vector
	const Ev & x;

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param key which element we are modifying
	 * \param grid_dst grid we are updating
	 * \param v Source Eigen vector
	 *
	 */
	inline copy_ele(const grid_dist_key_dx<dim> & key, S & grid_dst, const Ev & x)
	:key(key),grid_dst(grid_dst),x(x){};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline copy_ele(grid_dist_key_dx<dim> & key, S & grid_dst, Ev && x)
	:key(key),grid_dst(grid_dst),x(x)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at_c<S::value_type::type,T::value>::type copy_type;

/*		for (size_t i = 0 ; i < std::extents<copy_type>::value ; i++)
		{
			if ()
				grid_dst.template get<T::value>(key)[i] = x(lin_id * eqs_sys::nvar + base_id + i);
			else
				grid_dst.template get<T::value>(key)[i] = x(lin_id * eqs_sys::nvar + base_id + i);
		}*/
	}
};



template<typename T>
class Vector<T,Eigen::Matrix<T, Eigen::Dynamic, 1>>
{
	Eigen::Matrix<T, Eigen::Dynamic, 1> v;

public:

	/*! \brief Resize the Vector
	 *
	 * \param row numbers of row
	 *
	 */
	void resize(size_t row)
	{
		v.resize(row);
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	T & get(size_t i)
	{
		return v(i);
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	const Eigen::Matrix<T, Eigen::Dynamic, 1> & getVec() const
	{
		return v;
	}

	/*! \brief Get the Eigen Vector object
	 *
	 * \return the Eigen Vector
	 *
	 */
	Eigen::Matrix<T, Eigen::Dynamic, 1> & getVec()
	{
		return v;
	}

	/*! \brief Copy the vector into the grid
	 *
	 *
	 */
	template<typename Eqs_sys,typename Grid ,unsigned int ... pos> void copy(Grid & g, const grid_dist_iterator_sub<Eqs_sys::dims,typename Grid::d_grid> & g_i)
	{
		auto g_it = g_i;

		while (g_it.isNext() == true)
		{
			typedef typename to_boost_vmpl<pos...>::type vid;
			typedef boost::mpl::size<vid> v_size;

			copy_ele<Eqs_sys::dims,Grid,Eigen::Matrix<T, Eigen::Dynamic, 1>> cp(g_it.get(),g,v);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

			++g_it;
		}
	}
};


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_ */
