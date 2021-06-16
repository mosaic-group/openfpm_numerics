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
			h = 0;
			std::cout << "Only gaussian derivatives of order 0, 1, and 2 implemented. Aborting..." << std::endl;
			abort();
			break;
	}
	return h;
}

template <typename point_type>
double gaussian(const point_type & x, const double mu, const double sigma, const int order=0)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	
	double sum = 0;
	double normalization_factor = 1;
	double h = 1;
	for(int d=0; d<point_type::dims; d++)
	{
		sum += (x.get(d) - mu) * (x.get(d) - mu) / (sigma*sigma);;
		normalization_factor *= 1 / (sqrt2pi * sigma);
		h *= hermite_polynomial(x.get(d), sigma, order);
	}
	return  h * normalization_factor * exp(-0.5 * sum);
}


double gaussian(const double x, const double mu, const double sigma, const int order=0)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	
	double sum = (x - mu) * (x - mu) / (sigma*sigma);;
	double normalization_factor = 1 / (sqrt2pi * sigma);
	double h = hermite_polynomial(x, sigma, order);

	return  h * normalization_factor * exp(-0.5 * sum);
}




#endif //OPENFPM_NUMERICS_GAUSSIAN_HPP
