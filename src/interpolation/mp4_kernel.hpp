/*
 * mp4_kernel.hpp
 *
 *  Created on: May 5, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_INTERPOLATION_MP4_KERNEL_HPP_
#define OPENFPM_NUMERICS_SRC_INTERPOLATION_MP4_KERNEL_HPP_

#include <iostream>

template<typename st>
class mp4_kernel
{
public:

	static const int np = 4;

	static inline st value(st x, size_t i)
	{
		if (i == 0)
			return st(2.0) + (st(-4.0)+(st(2.5)-st(0.5)*-x)*-x)*-x;
		else if (i == 1)
			return st(1.0) + (st(-2.5)+st(1.5)*-x)*x*x;
		else  if (i == 2)
			return st(1.0) + (st(-2.5)+st(1.5)*x)*x*x;
		else if (i == 3)
			return st(2.0) + (st(-4.0)+(st(2.5)-st(0.5)*x)*x)*x;
		return 0.0;
	}
};

#endif /* OPENFPM_NUMERICS_SRC_INTERPOLATION_MP4_KERNEL_HPP_ */
