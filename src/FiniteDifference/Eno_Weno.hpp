/*
 * 
 * Based on the description in Osher and Fedkiw, Level Set Methods and Dynamic Implicit Surfaces
 * For HJ ENO method - Sec. 3.3
 * For HJ WENO method - Sec. 3.4
*/

// Created by Sachin, modified by Justina

#ifndef __ENO_WENO_HPP__
#define __ENO_WENO_HPP__

#include "Grid/grid_dist_id.hpp"

static double adjustWeights(double v1, double v2, double v3, double v4, double v5)
{
	double phix1, phix2, phix3;
	double s1, s2, s3;
	double a1, a2, a3;
	double eps = 1e-6;
	double w1, w2, w3;
	
	// Eqs. (3.25), (3.26), (3.27)
	phix1 = (v1 / 3.0) - (7.0 * v2 / 6.0) + (11.0 * v3 / 6.0);
	phix2 = -(v2 / 6.0) + (5.0 * v3 / 6.0) + (v4 / 3.0);
	phix3 = (v3 / 3.0) + (5.0 * v4 / 6.0) - (v5 / 6.0);
	
	// Eqs. (3.32), (3.33), (3.34)
	s1 = (13.0 / 12.0) * (v1 - 2*v2 + v3) * (v1 - 2*v2 + v3)
			+ 0.25 * (v1 - 4*v2 + 3*v3) * (v1 - 4*v2 + 3*v3);
	s2 = (13.0 / 12.0) * (v2 - 2*v3 + v4) * (v2 - 2*v3 + v4)
			+ 0.25 * (v2 - v4) * (v2 - v4);
	s3 = (13.0 / 12.0) * (v3 - 2*v4 + v5) * (v3 - 2*v4 + v5)
			+ 0.25 * (3*v3 - 4*v4 + v5) * (3*v3 - 4*v4 + v5);
	
	// Eqs. (3.35), (3.36), (3.37)
	a1 = 0.1 / ((s1 + eps)*(s1 + eps));
	a2 = 0.6 / ((s2 + eps)*(s2 + eps));
	a3 = 0.3 / ((s3 + eps)*(s3 + eps));
	
	// Eqs. (3.39), (3.40), (3.41)
	w1 = a1 / (a1 + a2 + a3);
	w2 = a2 / (a1 + a2 + a3);
	w3 = a3 / (a1 + a2 + a3);
	
	// Eq. (3.28)
	return (w1 * phix1 + w2 * phix2 + w3 * phix3);
}

template<size_t Field, typename grid_type, typename key_type>
double WENO_5_Plus(grid_type & grid, key_type key, size_t x)
{
	double c2, c3, c4, c5, c6;
	double coeff = 1.0 / grid.spacing(x);
	double df;
	
	c2 = (grid.template get<Field>(key.move(x, -1)) - grid.template get<Field>(key.move(x, -2))) * coeff;
	c3 = (grid.template get<Field>(key) - grid.template get<Field>(key.move(x, -1))) * coeff;
	c4 = (grid.template get<Field>(key.move(x,1)) - grid.template get<Field>(key)) * coeff;
	c5 = (grid.template get<Field>(key.move(x, 2)) - grid.template get<Field>(key.move(x, 1))) * coeff;
	c6 = (grid.template get<Field>(key.move(x, 3)) - grid.template get<Field>(key.move(x, 2))) * coeff;
	df = adjustWeights(c6, c5, c4, c3, c2);
	
	return df;
}

template<size_t Field, typename grid_type, typename key_type>
double WENO_5_Minus(grid_type & grid, key_type key, size_t x)
{
	double c1, c2, c3, c4, c5;
	double coeff = 1.0 / grid.spacing(x);
	double df;
	
	c1 = (grid.template get<Field>(key.move(x, -2)) - grid.template get<Field>(key.move(x, -3))) * coeff;
	c2 = (grid.template get<Field>(key.move(x, -1)) - grid.template get<Field>(key.move(x, -2))) * coeff;
	c3 = (grid.template get<Field>(key) - grid.template get<Field>(key.move(x, -1))) * coeff;
	c4 = (grid.template get<Field>(key.move(x,1)) - grid.template get<Field>(key)) * coeff;
	c5 = (grid.template get<Field>(key.move(x, 2)) - grid.template get<Field>(key.move(x, 1))) * coeff;
	
	df = adjustWeights(c1, c2, c3, c4, c5);
	
	return df;
}

template<size_t Field, typename grid_type, typename key_type>
double ENO_3_Plus(grid_type & grid, key_type key, size_t x)
{
	double q1x, q2x, q3x;
	double gridsize = grid.spacing(x);
	double coeff = 1.0 / gridsize;
	
	q1x = (grid.template get<Field>(key.move(x,1)) - grid.template get<Field>(key)) * coeff;
	double d2ix = (grid.template get<Field>(key.move(x,1)) - 2*grid.template get<Field>(key) + grid.template get<Field>(key.move(x,-1))) * 0.5 * coeff * coeff;
	double d2ip1x = (grid.template get<Field>(key.move(x,2)) - 2*grid.template get<Field>(key.move(x,1)) + grid.template get<Field>(key)) * 0.5 * coeff * coeff;
	if(abs(d2ix) <= abs(d2ip1x))
	{
		q2x = -d2ix * gridsize; // (3.22)
		double d3p1hx = (grid.template get<Field>(key.move(x,1)) - 3*grid.template get<Field>(key) + 3*grid.template get<Field>(key.move(x,-1)) - grid.template get<Field>(key.move(x,-2))) * coeff * coeff * coeff / 6.0;
		double d3p3hx = (grid.template get<Field>(key.move(x,2)) - 3*grid.template get<Field>(key.move(x,1)) + 3*grid.template get<Field>(key) - grid.template get<Field>(key.move(x,-1))) * coeff * coeff * coeff / 6.0;
		if(abs(d3p1hx) <= abs(d3p3hx))
			q3x = -d3p1hx * gridsize * gridsize;
		else
			q3x = -d3p3hx * gridsize * gridsize;
	}
	else
	{
		q2x = -d2ip1x * gridsize;
		double d3p1hx = (grid.template get<Field>(key.move(x,2)) - 3*grid.template get<Field>(key.move(x,1)) + 3*grid.template get<Field>(key) - grid.template get<Field>(key.move(x,-1))) * coeff * coeff * coeff / 6.0;
		double d3p3hx = (grid.template get<Field>(key.move(x,3)) - 3*grid.template get<Field>(key.move(x,2)) + 3*grid.template get<Field>(key.move(x,1)) - grid.template get<Field>(key)) * coeff * coeff * coeff / 6.0;
		if(abs(d3p1hx) <= abs(d3p3hx))
			q3x = d3p1hx * 2 * gridsize * gridsize;
		else
			q3x = d3p3hx * 2 * gridsize * gridsize;
		
	}
	return (q1x + q2x + q3x);
}

template<size_t Field, typename grid_type, typename key_type>
double ENO_3_Minus(grid_type & grid, key_type key, size_t x)
{
	double q1x, q2x, q3x;
	double gridsize = grid.spacing(x);
	double coeff = 1.0 / gridsize;
	
	q1x = (grid.template get<Field>(key) - grid.template get<Field>(key.move(x, -1))) * coeff;
	double d2im1x = (grid.template get<Field>(key) - 2*grid.template get<Field>(key.move(x, -1)) + grid.template get<Field>(key.move(x, -2))) * 0.5 * coeff * coeff;
	double d2ix = (grid.template get<Field>(key.move(x,1)) - 2*grid.template get<Field>(key) + grid.template get<Field>(key.move(x,-1))) * 0.5 * coeff * coeff;
	
	if(abs(d2im1x) <= abs(d2ix))
	{
		q2x = d2im1x * gridsize;
		double d3p1hx = (grid.template get<Field>(key) - 3*grid.template get<Field>(key.move(x,-1)) + 3*grid.template get<Field>(key.move(x,-2)) - grid.template get<Field>(key.move(x,-3))) * coeff * coeff * coeff / 6.0;
		double d3p3hx = (grid.template get<Field>(key.move(x,1)) - 3*grid.template get<Field>(key) + 3*grid.template get<Field>(key.move(x,-1)) - grid.template get<Field>(key.move(x,-2))) * coeff * coeff * coeff / 6.0;
		if(abs(d3p1hx) <= abs(d3p3hx))
			q3x = d3p1hx * 2 * gridsize * gridsize;
		else
			q3x = d3p3hx * 2 * gridsize * gridsize;
	}
	else
	{
		q2x = d2ix * gridsize;
		double d3p1hx = (grid.template get<Field>(key.move(x,1)) - 3*grid.template get<Field>(key) + 3*grid.template get<Field>(key.move(x,-1)) - grid.template get<Field>(key.move(x,-2))) * coeff * coeff * coeff / 6.0;
		double d3p3hx = (grid.template get<Field>(key.move(x,2)) - 3*grid.template get<Field>(key.move(x,1)) + 3*grid.template get<Field>(key) - grid.template get<Field>(key.move(x,-1))) * coeff * coeff * coeff / 6.0;
		if(abs(d3p1hx) <= abs(d3p3hx))
			q3x = -d3p1hx * gridsize * gridsize;
		else
			q3x = -d3p3hx * gridsize * gridsize;
	}
	
	return (q1x + q2x + q3x);
}

#endif

