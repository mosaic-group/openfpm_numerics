//
// Created by jstark on 2020-05-08.
//
/**
 * @file ComputeGradient.hpp
 *
 * @brief Header file containing the functions to compute an upwind gradient and the magnitude of gradient.
 *
 * @details Upwind gradient is computed according to <a href="https://www.researchgate.net/publication/2642502_An_Efficient_Interface_Preserving_Level_Set_Re
 * -Distancing_Algorithm_And_Its_Application_To_Interfacial_Incompressible_Fluid_Flow.html">M. Sussman and E. Fatemi,
 * “Efficient, interface-preserving level set redistancing algorithm and its application to interfacial
 * incompressible fluid flow” (1999)</a> ).
 *
 * @author Justina Stark
 * @date May 2020
 */
#ifndef REDISTANCING_SUSSMAN_COMPUTEGRADIENT_HPP
#define REDISTANCING_SUSSMAN_COMPUTEGRADIENT_HPP

// Include standard library header files
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cmath>
#include <array>
// Include OpenFPM header files
#include "Grid/grid_dist_id.hpp"

/**@brief Computes the forward finite difference of a scalar property on the current grid node.
 *
 * @tparam Phi Size_t index of property for which the gradient should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @return Forward finite difference for the property under index Phi on the current node with index key.
 */
template <size_t Phi, typename gridtype, typename keytype>
double forward_FD(gridtype & grid, keytype & key, size_t d)
{
	return (grid.template get<Phi> (key.move(d, 1)) - grid.template get<Phi> (key)) / grid.getSpacing()[d];
}
/**@brief Computes the backward finite difference of a scalar property on the current grid node.
 *
 * @tparam Phi Size_t index of property for which the gradient should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @return Backward finite difference for the property under index Phi on the current node with index key.
 */
template <size_t Phi, typename gridtype, typename keytype>
double backward_FD(gridtype & grid, keytype & key, size_t d)
{
	return (grid.template get<Phi> (key) - grid.template get<Phi> (key.move(d, -1))) / grid.getSpacing()[d];
}
/**@brief Get the first order upwind finite difference of a scalar property on the current grid node.
 *
 * @details Calls #forward_FD() and #backward_FD() and chooses between dplus and dminus by upwinding.
 *
 * @tparam Phi Size_t index of property for which the gradient should be computed.
 * @tparam Phi_0_sign Size_t index of property that contains the initial sign of that property for which the upwind FD
 *                    should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @return First order upwind finite difference for the property under index Phi on the current node with index key.
 */
template <size_t Phi, size_t Phi_0_sign, typename gridtype, typename keytype>
double upwind_FD(gridtype &grid, keytype &key, size_t d)
{
	double dplus = 0, dminus = 0, dphi = 0;
	
	dplus = forward_FD<Phi>(grid, key, d);
	dminus = backward_FD<Phi>(grid, key, d);
	int sign_phi0 = grid.template get<Phi_0_sign> (key);
	
	if (dplus * sign_phi0 < 0
	    && (dminus + dplus) * sign_phi0 < 0) dphi = dplus;
	
	else if (dminus * sign_phi0 > 0
	         && (dminus + dplus) * sign_phi0 > 0) dphi = dminus;
	
	else if (dminus * sign_phi0 < 0
	         && dplus * sign_phi0 > 0) dphi = 0;
	
	//else if (-10e-12 <= dminus <= 10e-12 && -10e-12 <= dplus <= 10e-12) dphi = 0;
	
	return dphi;
}
/**@brief Computes 1st order upwind finite difference inside and forward/backward finite difference at the grid borders.
 *
 * @details Checks if a point lays within the grid or at the border. For the internal grid points it calls upwind
 * finite difference #upwind_FD(). For the border points, simply #forward_FD() and #backward_FD() is used,
 * respectively,
 * depending on the side of the border.
 *
 * @tparam Phi Size_t index of property for which the gradient should be computed.
 * @tparam Phi_0_sign Size_t index of property that contains the initial sign of that property for which the upwind FD
 *                    should be computed.
 * @tparam Phi_grad Size_t index of property where the gradient result should be stored.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the gradient should be computed.
 */
// Use upwinding for inner grid points and one sided backward / forward stencil at border grid points
template <size_t Phi, size_t Phi_0_sign, size_t Phi_grad, typename gridtype>
void get_first_order_gradient_depending_on_point_position(gridtype &grid)
{
	grid.template ghost_get<Phi>();
	auto dom = grid.getDomainIterator();
	while (dom.isNext())
	{
		auto key = dom.get();
		auto key_g = grid.getGKey(key);
		
		for(size_t d = 0; d < gridtype::dims; d++ )
		{
			if (key_g.get(d) > 0 && key_g.get(d) < grid.size(d) - 1) // if point lays not at the border of the grid
			{
				grid.template get<Phi_grad> (key) [d] = upwind_FD<Phi, Phi_0_sign>(grid, key, d);
			}
			
			else if (key_g.get(d) == 0) // if point lays at left border, use right sided kernel
			{
				grid.template get<Phi_grad> (key) [d] = forward_FD<Phi>(grid, key, d);
			}
			
			else if (key_g.get(d) >= grid.size(d) - 1) // if point lays at right border, use left sided kernel
			{
				grid.template get<Phi_grad> (key) [d] = backward_FD<Phi>(grid, key, d);
			}
		}
		++dom;
	}
}
/**@brief Calls #get_first_order_gradient_depending_on_point_position(). Computes 1st order upwind finite difference.
 *
* @tparam Phi Size_t index of property for which the gradient should be computed.
 * @tparam Phi_0_sign Size_t index of property that contains the initial sign of that property for which the upwind FD
 *                    should be computed.
 * @tparam Phi_grad Size_t index of property where the gradient result should be stored.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the gradient should be computed.
 */
/*Approximate the initial gradients between neighbor particles using fwd and bwd differences.
 These gradients are needed for the re-distancing step according to Sussman & Fatemi 1999.*/
template <size_t Phi_in, size_t Phi_0_sign, size_t Phi_grad_out, typename gridtype>
void get_upwind_gradient(gridtype & grid)
{
	get_first_order_gradient_depending_on_point_position<Phi_in, Phi_0_sign, Phi_grad_out>(grid);
}

/**@brief Computes the magnitude of the gradient (L2-norm of gradient vector).
 *
 * @tparam Phi_grad_in Size_t index of property that contains the gradient.
 * @tparam Phi_magnOfGrad_out Size_t index of property where the magnitude of gradient should be stored.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the magnitude of gradient should be computed.
 */
template <size_t Phi_grad_in, size_t Phi_magnOfGrad_out, typename gridtype>
void get_gradient_magnitude(gridtype & grid)
{
	grid.template ghost_get<Phi_grad_in>();
	auto dom = grid.getDomainGhostIterator();
	while(dom.isNext())
	{
		double sum = 0;
		auto key = dom.get();
		for(size_t d = 0; d < gridtype::dims; d++)
		{
			sum += grid.template get<Phi_grad_in> (key)[d] * grid.template get<Phi_grad_in> (key)[d];
		}
		grid.template get<Phi_magnOfGrad_out> (key) = sqrt(sum);
		++dom;
	}
}



#endif //REDISTANCING_SUSSMAN_COMPUTEGRADIENT_HPP
