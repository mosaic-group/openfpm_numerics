//
// Created by jstark on 2020-08-18.
//
/**
 * @file Circle.hpp
 *
 * @brief Header file containing functions for creating a filled 2D circle of defined radius.
 *
 * @author Justina Stark
 * @date August 2020
 */
#ifndef ACCURACY_TESTS_CIRCLE_HPP
#define ACCURACY_TESTS_CIRCLE_HPP

#include <iostream>

/**@brief Checks if a point lays inside a circle of given radius and given center coordinates.
 *
 * @tparam coord_type Inferred type of the point coordinates x and y.
 * @tparam radius_type Inferred type of radius the circle is supposed to have.
 * @param x_coord X-coordinate of the point for that laying inside the circle should be checked.
 * @param y_coord Y-coordinate of the point for that laying inside the circle should be checked.
 * @param radius Radius of the filled circle.
 * @param center_x X-coordinate of the circle center.
 * @param center_y Y-coordinate of the circle center.
 * @return True, if point lays inside the circle; false, if outside.
 */
template <typename coord_type, typename radius_type>
bool inside_circle(coord_type x_coord, coord_type y_coord, radius_type radius, float center_x=0, float center_y=0)
{
	return (x_coord - center_x) * (x_coord - center_x) + (y_coord - center_y) * (y_coord - center_y) <= radius * radius;
}

/**@brief Creates a filled circle of given radius an given center coordinates on a 2D OpenFPM grid.
 *
 * @details The circle is represented in one of the grid properties via the following indicator function:
 * @f[ \phi_{\text{indicator}} = \begin{cases}
 *    +1 & \text{point lies inside the circle} \\
 *    -1 & \text{point lies outside the circle} \\
 * \end{cases} @f]
 *
 * @tparam Phi_0 Index of grid property that should store the indicator function representing the filled circle.
 * @tparam grid_type Inferred type of the input grid.
 * @tparam radius_type Inferred type of radius the circle is supposed to have.
 * @param grid Input OpenFPM grid.
 * @param radius Radius of the filled circle.
 * @param center_x X-coordinate of the circle center.
 * @param center_y Y-coordinate of the circle center.
 */
template <size_t Phi_0, typename grid_type, typename radius_type>
void init_grid_with_circle(grid_type & grid, radius_type radius, float center_x=0, float center_y=0)
{
	const size_t x = 0;
	const size_t y = 1;
	// assign pixel values to domain. For each pixel get factor_refinement number of grid points with corresponding value
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		auto gkey = grid.getGKey(key);
		auto spacing = grid.getSpacing();
		auto i = gkey.get(x);
		auto j = gkey.get(y);
		double x_coord = i * spacing[x];
		double y_coord = j * spacing[y];
		
		if (inside_circle(x_coord, y_coord, radius, center_x, center_y))
		{
			grid.template get<Phi_0> (key) = 1;
		}
		else
		{
			grid.template get<Phi_0> (key) = -1;
		}
		++dom;
	}
}


#endif //ACCURACY_TESTS_CIRCLE_HPP
