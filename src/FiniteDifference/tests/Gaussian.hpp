//
// Created by jstark on 14.06.21.
//

#ifndef OPENFPM_NUMERICS_GAUSSIAN_HPP
#define OPENFPM_NUMERICS_GAUSSIAN_HPP

#include "cmath"

template <typename point_type>
double gaussian(const point_type & x, const double mu, const double sigma)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	const double sigmapow2 = sigma * sigma;
	
	double sum = 0;
	double normalization_factor = 1;
	for(int d=0; d<point_type::dims; d++)
	{
		sum += (x.get(d) - mu) * (x.get(d) - mu) / (sigmapow2);;
		normalization_factor *= 1 / (sqrt2pi * sigma);
	}
	
	return  normalization_factor * exp(-0.5 * sum);
}

// 1D
double gaussian(const double x, const double mu, const double sigma)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	const double sigmapow2 = sigma * sigma;
	
	double sum = 0;
	double normalization_factor = 1;

	sum += (x - mu) * (x - mu) / (sigmapow2);;
	normalization_factor *= 1 / (sqrt2pi * sigma);
	
	return  normalization_factor * exp(-0.5 * sum);
}

template <typename point_type>
double gaussian_dx(const point_type & x, const double mu, const double sigma)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	
	const double sigmapow2 = sigma * sigma;
	
	double sum = (x - mu) * (x - mu) / (sigmapow2);
	double normalization_factor = 1 /( sqrt2pi * sigma);
	double prefactor = - x / sigmapow2;
	
	return prefactor *  normalization_factor * exp(-0.5 * sum);
	
}

template <typename point_type>
double gaussian_ddx(const point_type & x, const double mu, const double sigma)
{
	const double pi = 3.14159265358979323846;
	const double sqrt2pi = sqrt(2*pi);
	
	const double sigmapow2 = sigma * sigma;
	const double sigmapow4 = sigmapow2 * sigmapow2;
	
	double sum = (x - mu) * (x - mu) / (sigmapow2);
	double normalization_factor = 1 /( sqrt2pi * sigma);
	double prefactor = - (sigmapow2 - x*x) / sigmapow4;
	
	return prefactor *  normalization_factor * exp(-0.5 * sum);
}


#endif //OPENFPM_NUMERICS_GAUSSIAN_HPP
