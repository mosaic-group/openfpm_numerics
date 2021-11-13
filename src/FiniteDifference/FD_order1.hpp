//
// Created by jstark on 17.05.21.
//

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FD_ORDER1_HPP
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FD_ORDER1_HPP

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
double FD_forward(gridtype & grid, keytype & key, size_t d)
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
double FD_backward(gridtype & grid, keytype & key, size_t d)
{
	return (grid.template get<Field> (key) - grid.template get<Field> (key.move(d, -1))) / grid.getSpacing()[d];
}







#endif //OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FD_ORDER1_HPP

