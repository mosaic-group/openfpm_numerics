//
// Created by jstark on 17.05.21.
//
/**
 * @file Upwind_gradient.hpp
 *
 * @brief Approximating upwind gradients on a grid with the following options for the order of accuracy: 1, 3 or 5
 *
 * @details For solving hyperbolic PDEs using differencing that is biased by the direction of the wave front. Upwinding
 * is done according to <a href="https://www.researchgate
 * .net/publication/2642502_An_Efficient_Interface_Preserving_Level_Set_Re
 * -Distancing_Algorithm_And_Its_Application_To_Interfacial_Incompressible_Fluid_Flow.html">M. Sussman and E. Fatemi,
 * “Efficient, interface-preserving level set redistancing algorithm and its application to interfacial
 * incompressible fluid flow” (1999)</a> ), paragraph 4.1, 2 b). Order 1 is achieved by upwinding of forward and
 * backward finite difference, order 3 and 5 by upwinding the ENO and WENO scheme, respectively.
 *
 * @author Justina Stark
 * @date May 2021
 */
#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_UPWIND_GRADIENT_HPP
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_UPWIND_GRADIENT_HPP

// Include standard libraries
#include <cmath>

// Include OpenFPM header files
#include "Grid/grid_dist_id.hpp"
#include "FD_order1.hpp"
#include "Eno_Weno.hpp"

/**@brief Upwinding: For a specific dimension, from the forward and backward gradient find the upwind side.
 *
 * @param dplus: Gradient approximated using RHS neighbors.
 * @param dminus: Gradient approximated using LHS neighbors.
 * @param sign: Sign of the velocity with which the wave front is moving.
 * @return Scalar double upwind gradient approximation in the dimension given.
 */
double upwinding(double dplus, double dminus, int sign)
{
	double grad_upwind = 0;
	if (dplus * sign < 0
			&& (dminus + dplus) * sign < 0) grad_upwind = dplus;
	
	else if (dminus * sign > 0
			&& (dminus + dplus) * sign > 0) grad_upwind = dminus;
	
	else if (dminus * sign < 0
			&& dplus * sign > 0) grad_upwind = 0;
	return grad_upwind;
}


/**@brief Get the upwind finite difference of a scalar property on the current grid node.
 *
 * @details Order of accuracy can be 1, 3 or 5.
 *
 * @tparam Field Size_t index of property for which the gradient should be computed.
 * @tparam Sign Size_t index of property that contains the initial sign of that property for which the upwind FD
 *                    should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @param order Order of accuracy of the difference scheme. Can be 1, 3 or 5.
 * @return Upwind finite difference in one dimension of the property under index Field on the current node with index key.
 */
template <size_t Field, size_t Sign, typename gridtype, typename keytype>
double FD_upwind(gridtype &grid, keytype &key, size_t d, size_t order)
{
	double dplus = 0, dminus = 0;
	int sign_phi0 = grid.template get<Sign> (key);
	
	switch(order)
	{
		case 1:
			dplus  = FD_forward<Field>(grid, key, d);
			dminus = FD_backward<Field>(grid, key, d);
			break;
		case 3:
			dplus  = ENO_3_Plus<gridtype, keytype, Field>(key, grid, d);
			dminus = ENO_3_Minus<gridtype, keytype,Field>(key, grid, d);
			break;
		case 5:
			dplus  = WENO_5_Plus<gridtype, keytype, Field>(key, grid, d);
			dminus = WENO_5_Minus<gridtype, keytype, Field>(key, grid, d);
			break;
		default:
			auto &v_cl = create_vcluster();
			if (v_cl.rank() == 0) std::cout << "Order of accuracy chosen not valid. Using default order 1." <<
						std::endl;
			dplus  = FD_forward<Field>(grid, key, d);
			dminus = FD_backward<Field>(grid, key, d);
			break;
	}
	
	return upwinding(dplus, dminus, sign_phi0);
}


/**@brief Computes upwind gradient with order of accuracy 1, 3 or 5.
 *
 * @details Checks if a point lays within the grid or at the boundary. For the internal grid points, it calls upwind
 * finite difference #FD_upwind(). For the border points, simply #FD_forward() and #FD_backward() is used,
 * respectively, depending on the side of the border.
 *
 * @tparam Field Size_t index of property for which the gradient should be computed.
 * @tparam Sign Size_t index of property that contains the initial sign of that property for which the upwind FD
 *                    should be computed.
 * @tparam Gradient Size_t index of property where the gradient result should be stored.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the gradient should be computed.
 * @param one_sided_BC Bool variable, if true, use one-sided kernel for boundary-nodes. If false, extend stencil onto
 * ghost nodes.
 * @param order Size_t variable, order of accuracy the upwind FD scheme should have. Can be 1, 3 or 5.
 */
// Use upwinding for inner grid points and one sided backward / forward stencil at border (if one_sided_BC=true)
template <size_t Field, size_t Sign, size_t Gradient, typename gridtype>
void upwind_gradient(gridtype & grid, const bool one_sided_BC, size_t order)
{
	grid.template ghost_get<Field>(KEEP_PROPERTIES);
	grid.template ghost_get<Sign>(KEEP_PROPERTIES);
	
	auto dom = grid.getDomainIterator();
	
	if (one_sided_BC)
	{
		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = grid.getGKey(key);
			
			for(size_t d = 0; d < gridtype::dims; d++ )
			{
				// Grid nodes inside and distance from boundary > stencil-width
				if (key_g.get(d) > 2 && key_g.get(d) < grid.size(d) - 3) // if point lays with min. 3 nodes distance to
					// boundary
				{
					grid.template get<Gradient> (key) [d] = FD_upwind<Field, Sign>(grid, key, d, order);
				}
					
					// Grid nodes in stencil-wide vicinity to boundary
				else if (key_g.get(d) > 0 && key_g.get(d) < grid.size(d) - 1) // if point lays not on the grid boundary
				{
					grid.template get<Gradient> (key) [d] = FD_upwind<Field, Sign>(grid, key, d, 1);
				}
				
				else if (key_g.get(d) == 0) // if point lays at left boundary, use right sided kernel
				{
					grid.template get<Gradient> (key) [d] = FD_forward<Field>(grid, key, d);
				}
				
				else if (key_g.get(d) >= grid.size(d) - 1) // if point lays at right boundary, use left sided kernel
				{
					grid.template get<Gradient> (key) [d] = FD_backward<Field>(grid, key, d);
				}
			}
			++dom;
		}
	}
	
	else
	{
		while (dom.isNext())
		{
			auto key = dom.get();
			for(size_t d = 0; d < gridtype::dims; d++ )
			{
				grid.template get<Gradient> (key) [d] = FD_upwind<Field, Sign>(grid, key, d, order);
			}
			++dom;
		}
	}
}

/**@brief Checks if ghost layer is thick enough for a given stencil-width.
 *
 * @tparam gridtype Type of input grid.
 * @param grid Grid, for which ghost layer width should be checked.
 * @param required_width Size_t variable, number of ghost nodes that are required in one direction.
 * @return Bool, true if ghost layer sufficiently large, false otherwise.
 */
template <typename gridtype>
static bool ghost_width_is_sufficient(gridtype & grid, size_t required_width)
{
	auto ghost = grid.getDecomposition().getGhost();
	int np_ghost[gridtype::dims];
	
	for (int d = 0 ; d < gridtype::dims ; d++)
	{
		np_ghost[d] = ghost.getHigh(d) / grid.getSpacing()[d];
	}
	
	int min_width = np_ghost[0];
	for (int d = 0 ; d < gridtype::dims ; d++)
	{
		if(np_ghost[d] < min_width) min_width = np_ghost[d];
	}
	
	return (required_width >= min_width);
}

/**@brief Calls #upwind_gradient. Computes upwind gradient of desired order {1, 3, 5} for the whole n-dim grid.
 *
* @tparam Field Size_t index of property for which the gradient should be computed.
 * @tparam Sign Size_t index of property that contains the initial sign of that property for which the upwind FD
 *                    should be computed / sign of the velocity.
 * @tparam Gradient Size_t index of property where the upwind gradient result should be stored.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the gradient should be computed.
 */
template <size_t Field_in, size_t Sign, size_t Gradient_out, typename gridtype>
void get_upwind_gradient(gridtype & grid, const size_t order=5, const bool one_sided_BC=true)
{
	grid.template ghost_get<Field_in>();
	grid.template ghost_get<Sign>();
	
	if (!one_sided_BC)
	{
		auto &v_cl = create_vcluster();
		if (v_cl.rank() == 0)
		{
			// Check if ghost layer thick enough for stencil
			size_t stencil_width;
			(order > 1) ? stencil_width = 3 : stencil_width = 1;
			if (!ghost_width_is_sufficient(grid, stencil_width))
			{
				std::cout << "Error: Ghost layer not big enough. Either run with one_sided_BC=true or create a ghost "
				             "layer that has a width of least " << stencil_width << " grid node(s)" << std::endl;
				abort();
			}
		}
	}
	
	switch(order)
	{
		case 1:
		case 3:
		case 5:
			upwind_gradient<Field_in, Sign, Gradient_out>(grid, one_sided_BC, order);
			break;
		default:
			auto &v_cl = create_vcluster();
			if (v_cl.rank() == 0) std::cout << "Order of accuracy chosen not valid. Using default order 1." <<
						std::endl;
			upwind_gradient<Field_in, Sign, Gradient_out>(grid, one_sided_BC, 1);
			break;
	}
	
	
}


#endif //OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_UPWIND_GRADIENT_HPP

