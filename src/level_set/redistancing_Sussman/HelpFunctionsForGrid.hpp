//
// Created by jstark on 2020-05-12.
//
/**
 * @file HelpFunctionsForGrid.hpp
 *
 * @brief Header file containing help-functions that perform on OpenFPM-grids.
 *
 * @author Justina Stark
 * @date May 2020
 */
#ifndef REDISTANCING_SUSSMAN_HELPFUNCTIONSFORGRID_HPP
#define REDISTANCING_SUSSMAN_HELPFUNCTIONSFORGRID_HPP

#include <iostream>
#include <limits>

#include "HelpFunctions.hpp"
#include "VCluster/VCluster.hpp"
#include "Grid/grid_dist_id.hpp"

/**@brief Computes the time step for the iterative re-distancing fulfilling CFL condition.
 *
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 * @return Time step.
 */
template <typename grid_type>
typename grid_type::stype get_time_step_CFL(grid_type &grid)
{
	typename grid_type::stype sum = 0.0;
	for (size_t d = 0; d < grid_type::dims; d++)
	{
		sum += 1.0 / (grid.spacing(d) * grid.spacing(d));
	}
	return 0.5 / sum;
}

#if 0
/**@brief Computes the time step size fulfilling CFL condition according to https://www.cfd-online
 * .com/Wiki/Courant–Friedrichs–Lewy_condition for arbitrary dimensionality.
 *
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 * @param u Array of size grid_type::dims containing the velocity in each dimension.
 * @param Cmax Courant number.
 * @return Time step.
 */
template <typename grid_type>
typename grid_type::stype get_time_step_CFL(grid_type & grid, typename grid_type::stype u [grid_type::dims], typename grid_type::stype C)
{
	typename grid_type::stype sum = 0;
	for (size_t d = 0; d < grid_type::dims; d++)
	{
		sum += u[d] / grid.spacing(d);
	}
	return C / sum;
}
/**@brief Computes the time step size fulfilling CFL condition according to https://www.cfd-online
 * .com/Wiki/Courant–Friedrichs–Lewy_condition for arbitrary dimensionality.
 *
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 * @param u Velocity of propagating wave if isotropic for each direction.
 * @param Cmax Courant number.
 * @return Time step.
 */
template <typename grid_type>
typename grid_type::stype get_time_step_CFL(grid_type & grid, typename grid_type::stype u, typename grid_type::stype C)
{
	typename grid_type::stype sum = 0;
	for (size_t d = 0; d < grid_type::dims; d++)
	{
		sum += u / grid.spacing(d);
	}
	return C / sum;
}
#endif
/**@brief Initializes given property \p Prop of an OpenFPM grid including ghost layer with a given value from \p
 * init_value.
 *
 * @tparam Prop Index of property that should be initialized with the value from \p init_value.
 * @tparam grid_type Template type of input OpenFPM grid.
 * @tparam T Template type of the variable containing the initialization value \p init_value.
 * @param grid OpenFPM grid whose property \p Prop should be initialized, including its ghost layer.
 * @param init_value Variable that contains the value that should be copied to \p Prop.
 */
template <size_t Prop, typename grid_type, typename T>
void init_grid_and_ghost(grid_type & grid, T init_value)
{
	auto dom_ghost = grid.getDomainGhostIterator();
	while(dom_ghost.isNext())
	{
		auto key = dom_ghost.get();
		grid.template get<Prop>(key) = init_value;
		++dom_ghost;
	}
}

/**@brief Copies the value stored in a given property from a given source grid to a given destination grid.
 *
 * @tparam attr_sc Index of property that contains the value that should be copied.
 * @tparam attr_ds Index of property to which the value should be copied to.
 * @tparam grid_source_type Template type of the source grid \p grid_sc.
 * @tparam grid_dest_type Template type of the destination grid \p grid_ds.
 * @param grid_sc OpenFPM grid from which value is copied (source).
 * @param grid_ds OpenFPM grid to which value is copied to (destination).
 * @param include_ghost Bool variable that defines if the copying should include the ghost layer. False, if not; true,
 *                      if yes. Default: false.
 */
template <size_t attr_sc, size_t attr_ds, typename grid_source_type, typename grid_dest_type>
void copy_gridTogrid(const grid_source_type & grid_sc, grid_dest_type & grid_ds, bool include_ghost=false)
{
	assert(grid_source_type::dims == grid_dest_type::dims);
	assert(grid_sc.size() == grid_ds.size());
	if (include_ghost)
	{
		auto dom_sc = grid_sc.getDomainGhostIterator();
		auto dom_ds = grid_ds.getDomainGhostIterator();
		while (dom_sc.isNext())
		{
			auto key_sc = dom_sc.get();
			auto key_ds = dom_ds.get();
			grid_ds.template get<attr_ds>(key_ds) = grid_sc.template get<attr_sc>(key_sc);
			++dom_sc;
			++dom_ds;
		}
	}
	else
	{
		auto dom_sc = grid_sc.getDomainIterator();
		auto dom_ds = grid_ds.getDomainIterator();
		while (dom_sc.isNext())
		{
			auto key_sc = dom_sc.get();
			auto key_ds = dom_ds.get();
			grid_ds.template get<attr_ds>(key_ds) = grid_sc.template get<attr_sc>(key_sc);
			++dom_sc;
			++dom_ds;
		}
	}
}

/**@brief Determines the biggest spacing of a grid which is potentially anisotropic when comparing x, y (and z)
 * coordinate axis.
 *
 * @tparam grid_type Template type of input grid.
 * @param grid Input OpenFPM grid.
 * @return Biggest spacing of same type as grid-space.
 */
template <typename grid_type>
typename grid_type::stype get_biggest_spacing(grid_type & grid)
{
	typename grid_type::stype h_max = 0;
	for (size_t d = 0; d < grid_type::dims; d++)
	{
		if (grid.spacing(d) > h_max) h_max = grid.spacing(d);
	}
	return h_max;
}

/**@brief Determines the smallest spacing of a grid which is potentially anisotropic when comparing x, y (and z)
 * coordinate axis.
 *
 * @tparam grid_type Template type of input grid.
 * @param grid Input OpenFPM grid.
 * @return Smallest spacing of same type as grid-space.
 */
template <typename grid_type>
typename grid_type::stype get_smallest_spacing(grid_type & grid)
{
	typename grid_type::stype spacing [grid_type::dims];
	for (size_t d = 0; d < grid_type::dims; d++)
	{
		spacing[d] = grid.spacing(d);
	}
	std::sort(std::begin(spacing), std::end(spacing));
	return spacing[0];
}

/**@brief Computes the average difference between the values stored at two different properties of the same grid, that
 * is, the total difference of these values summed over all the grid nodes divided by the gridsize.
 *
 * @tparam Prop1 Index of the first property.
 * @tparam Prop2 Index of the second property.
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 * @return Sum over all grid nodes of the difference between the value stored at \p Prop1 and the value stored at \p
 * Prop2.
 */
template <size_t Prop1, size_t Prop2, typename grid_type>
auto average_difference(grid_type & grid)
{
	auto dom = grid.getDomainIterator();
	typedef typename std::decay_t<decltype(grid.template get<Prop1>(dom.get()))> prop_type;
	prop_type total_diff = 0;
	while (dom.isNext())
	{
		auto key = dom.get();
		total_diff += abs( grid.template get<Prop1>(key) - grid.template get<Prop2>(key) );
		++dom;
	}
	return total_diff / grid.size();
}

/**@brief Determines the maximum value stored on a given grid at a given property.
 *
 * @tparam Prop Index of property for which maximum should be found.
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 * @return Maximum value of property \p Prop in \p grid.
 */
template <size_t Prop, typename grid_type>
auto get_max_val(grid_type & grid)
{
	auto dom = grid.getDomainIterator();
	typedef typename std::decay_t<decltype(grid.template get<Prop>(dom.get()))> prop_type;
	prop_type max_value = std::numeric_limits<prop_type>::lowest();
	while(dom.isNext())
	{
		auto key = dom.get();
		if (grid.template get<Prop>(key) > max_value) max_value = grid.template get<Prop>(key);
		++dom;
	}
	auto & v_cl = create_vcluster();
	v_cl.max(max_value);
	v_cl.execute();
	return max_value;
}

/**@brief Determines the minimum value stored on a given grid at a given property.
 *
 * @tparam Prop Index of property for which minimum should be found.
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 * @return Minimum value of property \p Prop in \p grid.
 */
template <size_t Prop, typename grid_type>
auto get_min_val(grid_type & grid)
{
	auto dom = grid.getDomainIterator();
	typedef typename std::decay_t<decltype(grid.template get<Prop>(dom.get()))> prop_type;
	prop_type min_value = std::numeric_limits<prop_type>::max();
	while(dom.isNext())
	{
		auto key = dom.get();
		if (grid.template get<Prop>(key) < min_value) min_value = grid.template get<Prop>(key);
		++dom;
	}
	auto & v_cl = create_vcluster();
	v_cl.min(min_value);
	v_cl.execute();
	return min_value;
}

/**@brief Determines the sign of a value stored at a given property and stores it in another property.
 *
 * @tparam Prop_in Index of property that contains the value whose sign should be determined.
 * @tparam Prop_out Index of property where the sign should be written to.
 * @tparam grid_type Template type of the input grid.
 * @param grid Input OpenFPM grid.
 */
template <size_t Prop_in, size_t Prop_out, typename grid_type>
void init_sign_prop(grid_type & grid)
{
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		grid.template get<Prop_out>(key) = sgn(grid.template get<Prop_in>(key));
		++dom;
	}
}


/**@brief Computes the magnitude of the gradient (L2-norm of gradient vector).
 *
 * @tparam Phi_grad_in Size_t index of property that contains the gradient.
 * @tparam Phi_magnOfGrad_out Size_t index of property where the magnitude of gradient should be stored.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the magnitude of gradient should be computed.
 */
template <size_t Vector_in, size_t Magnitude_out, typename magnitude_type, typename gridtype>
void get_vector_magnitude(gridtype & grid)
{
	grid.template ghost_get<Vector_in>();
	auto dom = grid.getDomainGhostIterator();
	while(dom.isNext())
	{
		magnitude_type sum = 0;
		auto key = dom.get();
		for(size_t d = 0; d < gridtype::dims; d++)
		{
			sum += grid.template get<Vector_in> (key)[d] * grid.template get<Vector_in> (key)[d];
		}
		grid.template get<Magnitude_out> (key) = sqrt(sum);
		++dom;
	}
}

/**@brief Computes the magnitude of the gradient (L2-norm of gradient vector).
 *
 * @tparam Phi_grad_in Size_t index of property that contains the gradient.
 * @tparam gridtype Type of input grid.
 * @param grid Grid, on which the magnitude of gradient should be computed.
 */
template <size_t Vector_in, typename key_type, typename gridtype>
auto get_vector_magnitude(gridtype & grid, key_type & key)
{
	typedef typename std::remove_const_t<std::remove_reference_t<decltype(grid.template get<Vector_in>(key)[0])>>
	prop_type;
	prop_type sum = 0;
	for(size_t d = 0; d < gridtype::dims; d++)
	{
		sum += grid.template get<Vector_in> (key)[d] * grid.template get<Vector_in> (key)[d];
	}
	return sqrt(sum);
}


#endif //REDISTANCING_SUSSMAN_HELPFUNCTIONSFORGRID_HPP
