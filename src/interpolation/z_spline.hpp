/*
 * z_spline.hpp
 *
 *  Created on: May 8, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_INTERPOLATION_Z_SPLINE_HPP_
#define OPENFPM_NUMERICS_SRC_INTERPOLATION_Z_SPLINE_HPP_

#include "mp4_kernel.hpp"

template<typename st, unsigned int ord>
class z_kernel
{
public:

	static const int np = 4;

	static inline st value(st x, size_t i)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << ": Error this order has not been implemented" << std::endl;
		return 0.0;
	}
};

template<typename st>
class z_kernel<st,1>
{
public:

	static const int np = 2;

	static inline st value(st x, size_t i)
	{
		if (i == 0)
			return x + 1;
		else if (i == 1)
			return 1-x;
		return 0.0;
	}
};

template<typename st>
class z_kernel<st,2>: public mp4_kernel<st>
{

};

template<typename st>
class z_kernel<st,3>
{
public:

	static const int np = 6;

	static inline st value(st x, size_t i)
	{
		if (i == 0)
			return 18.0 + ( 38.25 + (31.875 + ( 313.0/24.0 + (2.625 + 5.0/24.0*x)*x)*x)*x)*x;
		else if (i == 1)
			return -4.0 + (-18.75 + (-30.625 + (-545.0/24.0 + ( -7.875 - 25.0/24.0*x)*x)*x)*x)*x;
		else if (i == 2)
			return 1.0 + (-1.25 +(35.0/12.0 +(5.25 + 25.0/12.0*x)*x)*x)*x*x;
		else if (i == 3)
			return 1.0 + (-1.25 +(-35.0/12.0 +(5.25 - 25.0/12.0*x)*x)*x)*x*x;
		else if (i == 4)
			return -4.0 + (18.75 + (-30.625 + (545.0/24.0 + ( -7.875 + 25.0/24.0*x)*x)*x)*x)*x;
		else if (i == 5)
			return 18.0 + (-38.25 + (31.875 + (-313.0/24.0 + (2.625 - 5.0/24.0*x)*x)*x)*x)*x;

		return 0.0;
	}
};

template<typename st>
class z_kernel<st,4>
{
public:

	static const int np = 8;

	static inline st value(st x, size_t i)
	{
		if (i == 0)
			return 726.4 + ( 1491.2 + (58786.0/45.0 + ( 633.0 + (26383.0/144.0 + (22807.0/720.0 + (727.0/240.0 + 89.0/720.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 1)
			return -440 +( -1297.45 + ( -117131/72.0 + ( -1123.5 + ( -66437.0/144.0 + ( -81109.0/720.0 + ( -727.0/48.0 - 623.0/720.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 2)
			return 27.6 + (8617.0/60.0 +(321.825 +(395.5 +( 284.8125 + (119.7875 +(27.2625 + 623.0/240.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 3)
			return 1.0 + (( -49.0/36.0 + (( -959.0/144.0 + ( -2569.0/144.0 + ( -727.0/48.0 - 623.0/144.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 4)
			return 1.0 + (( -49.0/36.0 + (( -959.0/144.0 + ( 2569.0/144.0 + ( -727.0/48.0 + 623.0/144.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 5)
			return 27.6 + (-8617.0/60.0 +(321.825 +(-395.5 +( 284.8125 + (-119.7875 +(27.2625 - 623.0/240.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 6)
			return -440 +( 1297.45 + ( -117131/72.0 + ( 1123.5 + ( -66437.0/144.0 + ( 81109.0/720.0 + ( -727.0/48.0 + 623.0/720.0*x)*x)*x)*x)*x)*x)*x;
		else if (i == 7)
			return 726.4 + ( -1491.2 + (58786.0/45.0 + ( -633.0 + (26383.0/144.0 + (-22807.0/720.0 + (727.0/240.0 - 89.0/720.0*x)*x)*x)*x)*x)*x)*x;

		return 0.0;
	}
};

#endif /* OPENFPM_NUMERICS_SRC_INTERPOLATION_Z_SPLINE_HPP_ */
