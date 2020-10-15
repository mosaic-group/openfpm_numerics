//
// Created by jstark on 2019-11-11.
//
/**
 * @file GaussFilter.hpp
 *
 * @brief Header file containing function which applies Gaussian smoothing (Gaussian blur) to an n-dim. grid.
 *
 * @details In cases of very steep gradients of the initial (pre-redistancing) phi at the interface, it can be
 * beneficial to smoothen phi prior to redistancing. If the initial function is a heaviside-like function, Gaussian
 * smoothing transforms it into a rather sigmoid-like functions.
 *
 * @author Justina Stark
 * @date November 2019
 */
#ifndef REDISTANCING_SUSSMAN_GAUSSFILTER_HPP
#define REDISTANCING_SUSSMAN_GAUSSFILTER_HPP


#include <iostream>
#include <typeinfo>
#include <cmath>

#include "Grid/grid_dist_id.hpp"
#include "HelpFunctions.hpp"
#include "HelpFunctionsForGrid.hpp"

/**@brief Discrete Gaussian filter with sigma=1.
 *
 * @details If Gaussian blur of sigma>1 is desired, just repeat smoothing.
 *
 * @tparam Phi_0_in Index of property that contains the values which should be Gaussian blurred.
 * @tparam Phi_smoothed_out Index of property where the output of smoothing should be stored.
 * @tparam grid_in_type Type of input grid.
 * @param grid_in Input-OpenFPM-grid on which Gaussian blur should be applied.
 */
template <size_t Phi_0_in, size_t Phi_smoothed_out, typename grid_in_type>
void gauss_filter(grid_in_type & grid_in)
{
	grid_in.template ghost_get<Phi_0_in>();
	
	// create a temporary grid with same size on which smoothing is perfomed to not override Phi_0 of the input grid
	typedef aggregate<double, double> props_temp;
	typedef grid_dist_id<grid_in_type::dims, typename grid_in_type::stype, props_temp> g_temp_type;
	g_temp_type g_temp(grid_in.getDecomposition(), grid_in.getGridInfoVoid().getSize(), Ghost<grid_in_type::dims, long int>(3));
	//	Some indices for better readability
	static constexpr size_t Phi_0_temp          = 0;
	static constexpr size_t Phi_smoothed_temp   = 1;
	
	copy_gridTogrid<Phi_0_in, Phi_0_temp>(grid_in, g_temp, true); // copy incl. ghost
	
	double gauss1d[] = {0.006, 0.061, 0.242, 0.383}; // 1D Gauss kernel for sigma = 1.0

//	loop over grid points. Gaussian kernel separated in its dimensions and first applied in x, then in y and z
	for (int d = 0; d < g_temp_type::dims; d++)
	{
		g_temp.template ghost_get<Phi_0_temp, Phi_smoothed_temp>();
		auto dom1 = g_temp.getDomainIterator();
		while (dom1.isNext())
		{
			auto key = dom1.get();
			auto key_g = g_temp.getGKey(key);
			
			g_temp.template get<Phi_smoothed_temp>(key) =
					gauss1d[0] * g_temp.template get<Phi_0_temp>(key.move(d, -3)) +
					gauss1d[1] * g_temp.template get<Phi_0_temp>(key.move(d, -2)) +
					gauss1d[2] * g_temp.template get<Phi_0_temp>(key.move(d, -1)) +
					gauss1d[3] * g_temp.template get<Phi_0_temp>(key) +
					gauss1d[2] * g_temp.template get<Phi_0_temp>(key.move(d, +1)) +
					gauss1d[1] * g_temp.template get<Phi_0_temp>(key.move(d, +2)) +
					gauss1d[0] * g_temp.template get<Phi_0_temp>(key.move(d, +3));
			++dom1;
		}
		copy_gridTogrid<Phi_smoothed_temp, Phi_0_temp>(g_temp, g_temp);
		g_temp.template ghost_get<Phi_0_temp, Phi_smoothed_temp>();
	}
	copy_gridTogrid<Phi_smoothed_temp, Phi_smoothed_out>(g_temp, grid_in, true); // copy incl. ghost
}


#endif //REDISTANCING_SUSSMAN_GAUSSFILTER_HPP
