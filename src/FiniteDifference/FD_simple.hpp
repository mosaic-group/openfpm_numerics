//
// Created by jstark on 17.05.21.
//
#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FD_SIMPLE_HPP
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FD_SIMPLE_HPP

// Include OpenFPM header files
#include "Grid/grid_dist_id.hpp"

/**@brief Computes the forward finite difference of a scalar property on the current grid node.
 *
 * @tparam Field Size_t index of property for which the gradient should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @return Forward finite difference for the property under index Field on the current node with index key.
 */
template <size_t Field, typename gridtype, typename keytype>
auto FD_forward(gridtype & grid, keytype & key, size_t d)
{
	return (grid.template get<Field> (key.move(d, 1)) - grid.template get<Field> (key)) / grid.getSpacing()[d];
}
/**@brief Computes the backward finite difference of a scalar property on the current grid node.
 *
 * @tparam Field Size_t index of property for which the gradient should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @return Backward finite difference for the property under index Field on the current node with index key.
 */
template <size_t Field, typename gridtype, typename keytype>
auto FD_backward(gridtype & grid, keytype & key, size_t d)
{
	return (grid.template get<Field> (key) - grid.template get<Field> (key.move(d, -1))) / grid.getSpacing()[d];
}


/**@brief Computes the central finite difference of a scalar field on the current grid node.
 *
 * @tparam Field Size_t index of property for which the gradient should be computed.
 * @tparam gridtype Type of input grid.
 * @tparam keytype Type of key variable.
 * @param grid Grid, on which the gradient should be computed.
 * @param key Key that contains the index of the current grid node.
 * @param d Variable (size_t) that contains the dimension.
 * @return Partial central finite difference for dimension d of Field on the current node with index key.
 */
template <size_t Field, typename gridtype, typename keytype>
auto FD_central(gridtype & grid, keytype & key, size_t d)
{
	return (grid.template get<Field>(key.move(d, 1))
			- grid.template get<Field>(key.move(d, -1)))
			/ (2 * grid.getSpacing()[d]);
}

/**@brief Computes the central finite difference of a scalar field on the full grid.
 * 
 * @tparam Field Size_t index of input property for which the gradient should be computed (scalar field).
 * @tparam Gradient Size_t index of output property (vector field).
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the gradient should be computed.
 * @param one_sided_BC Boolian. If true, stencil becoming 1st order one sided at the grid boundary (fwd/bwd FD). If 
 * false, stencil remains centered 2nd order at the boundary, using the ghost layer.
 */
template <size_t Field, size_t Gradient, typename gridtype>
void get_central_FD_grid(gridtype & grid, const bool one_sided_BC)
{
	grid.template ghost_get<Field>(KEEP_PROPERTIES);
	auto dom = grid.getDomainIterator();
	if (one_sided_BC)
	{
		while (dom.isNext())
		{
			auto key = dom.get();
			auto key_g = grid.getGKey(key);
			for(size_t d = 0; d < gridtype::dims; d++)
			{
				// Grid nodes inside and distance from boundary > stencil-width
				if (key_g.get(d) > 0 && key_g.get(d) < grid.size(d) - 1) // if point lays with min. 1 nodes distance to
					// boundary
				{
					grid.template get<Gradient> (key) [d] = FD_central<Field>(grid, key, d);
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
			for(size_t d = 0; d < gridtype::dims; d++)
			{
				grid.template get<Gradient> (key) [d] = FD_central<Field>(grid, key, d);
			}
			++dom;
		}
	}
}


#endif //OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FD_SIMPLE_HPP

