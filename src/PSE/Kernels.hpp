/*
 * Kernels.hpp
 *
 *  Created on: Jan 12, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_PSE_KERNELS_HPP_
#define OPENFPM_NUMERICS_SRC_PSE_KERNELS_HPP_

#include <boost/math/constants/constants.hpp>

// Gaussian kernel
#define KER_GAUSSIAN 1

/*! \brief Implementation of the Laplacian kernels for PSE
 *
 * \tparam dim Dimension
 * \tparam T type
 * \ord order pf approximation (default 2)
 * \impl TYPE of kernel
 *
 */
template<unsigned int dim, typename T, unsigned int ord=2, unsigned int impl=KER_GAUSSIAN>
struct Lap
{
	T epsilon;

	Lap(T epsilon)
	:epsilon(epsilon)
	{}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[dim], T (&y)[dim])
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " The laplacian for order:" << ord << " in dimension " << dim << " has not been implemented";
		return 0.0;
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[dim], const Point<dim,T> & y)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " The laplacian for order:" << ord << " in dimension " << dim << " has not been implemented";
		return 0.0;
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	double value(const Point<dim,T> & x, T (&y)[dim])
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " The laplacian for order:" << ord << " in dimension " << dim << " has not been implemented";
		return 0.0;
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	double value(const Point<dim,T> & x, const Point<dim,T> & y)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " The laplacian for order:" << ord << " in dimension " << dim << " has not been implemented";
		return 0.0;
	}
};

template<typename T>
struct Lap<1,T,2,KER_GAUSSIAN>
{
	T epsilon;

	inline Lap(T epsilon)
	:epsilon(epsilon)
	{}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y[i]) * (x[i] - y[i]);
		d = sqrt(d) / epsilon;

		return 4.0 / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y.get(i)) * (x[i] - y.get(i));
		d = sqrt(d) / epsilon;

		return 4.0 / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y[i]) * (x.get(i) - y[i]);
		d = sqrt(d) / epsilon;

		return 4.0 / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y.get(i)) * (x.get(i) - y.get(i));
		d = sqrt(d) / epsilon;

		return 4.0 / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d);
	}
};

template<typename T>
struct Lap<1,T,4,KER_GAUSSIAN>
{
	T epsilon;

	inline Lap(T epsilon)
	:epsilon(epsilon)
	{}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y[i]) * (x[i] - y[i]);
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (10.0-4.0*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y.get(i)) * (x[i] - y.get(i));
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (10.0-4.0*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y[i]) * (x.get(i) - y[i]);
		d = sqrt(d) / epsilon;

		return d / epsilon/ boost::math::constants::root_pi<T>() * exp(-d*d) * (10.0-4.0*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y.get(i)) * (x.get(i) - y.get(i));
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (10.0-4.0*d*d);
	}
};

template<typename T>
struct Lap<1,T,6,KER_GAUSSIAN>
{
	T epsilon;

	inline Lap(T epsilon)
	:epsilon(epsilon)
	{}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y[i]) * (x[i] - y[i]);
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (35.0/2.0-14.0*d*d+2.0*d*d*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y.get(i)) * (x[i] - y.get(i));
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (35.0/2.0-14.0*d*d+2.0*d*d*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y[i]) * (x.get(i) - y[i]);
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (35.0/2.0-14.0*d*d+2.0*d*d*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y.get(i)) * (x.get(i) - y.get(i));
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (35.0/2.0-14.0*d*d+2.0*d*d*d*d);
	}
};

template<typename T>
struct Lap<1,T,8,KER_GAUSSIAN>
{
	T epsilon;

	inline Lap(T epsilon)
	:epsilon(epsilon)
	{}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y[i]) * (x[i] - y[i]);
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (105.0/4.0-63.0/2.0*d*d+9.0*d*d*d*d-2.0/3.0*d*d*d*d*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(T (&x)[1], const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x[i] - y.get(i)) * (x[i] - y.get(i));
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (105.0/4.0-63.0/2.0*d*d+9.0*d*d*d*d-2.0/3.0*d*d*d*d*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, T (&y)[1])
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y[i]) * (x.get(i) - y[i]);
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (105.0/4.0-63.0/2.0*d*d+9.0*d*d*d*d-2.0/3.0*d*d*d*d*d*d);
	}

	/*! \brief From a kernel centered in x, it give the value of the kernel in y
	 *
	 * \param x center of the kernel
	 * \param y where we calculate the kernel
	 *
	 */
	inline double value(const Point<1,T> & x, const Point<1,T> & y)
	{
		double d = 0.0;
		for (size_t i = 0 ; i < 1 ; i++)
			d += (x.get(i) - y.get(i)) * (x.get(i) - y.get(i));
		d = sqrt(d) / epsilon;

		return d / epsilon / boost::math::constants::root_pi<T>() * exp(-d*d) * (105.0/4.0-63.0/2.0*d*d+9.0*d*d*d*d-2.0/3.0*d*d*d*d*d*d);
	}
};

#endif /* OPENFPM_NUMERICS_SRC_PSE_KERNELS_HPP_ */
