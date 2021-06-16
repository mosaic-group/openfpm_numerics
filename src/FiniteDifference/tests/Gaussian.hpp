//
// Created by jstark on 14.06.21.
//

#ifndef OPENFPM_NUMERICS_GAUSSIAN_HPP
#define OPENFPM_NUMERICS_GAUSSIAN_HPP

#include "cmath"


double hermite_polynomial(double x, double sigma, int order)
{
	double h;
	switch(order)
	{
		case 0:
			h = 1;
			break;
		case 1:
			h = -x / (sigma * sigma);
			break;
		case 2:
			h = (x*x - sigma*sigma) / (sigma*sigma*sigma*sigma);
			break;
		default:
			std::cout << "Only gaussian derivatives of order 0, 1, and 2 implemented. Aborting..." << std::endl;
			abort();
	}
	return h;
}

template <typename point_type>
double gaussian_1D(const point_type & x, const double mu, const double sigma)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	
	double sum = (x - mu) * (x - mu) / (sigma*sigma);;
	double normalization_factor = 1 / (sqrt2pi * sigma);

	return normalization_factor * exp(-0.5 * sum);
}

template <typename point_type>
double gaussian(const point_type & x, const double mu, const double sigma)
{
	double g = 1;
	for(int d=0; d<point_type::dims; d++)
	{
		g *= gaussian_1D(x.get(d), mu, sigma);
	}
	return  g;
}


#endif //OPENFPM_NUMERICS_GAUSSIAN_HPP
