//
// Created by jstark on 2020-05-17.
//
/**
 * @file Sphere.hpp
 *
 * @brief Header file containing functions for creating a filled 3D sphere of defined radius.
 *
 * @author Justina Stark
 * @date May 2020
 */
#ifndef ACCURACY_TESTS_SPHERE_HPP
#define ACCURACY_TESTS_SPHERE_HPP

#include <iostream>
/**@brief Checks if a point lays inside a sphere of given radius and given center coordinates.
 *
 * @tparam coord_type Inferred type of the point coordinates x, y and z.
 * @tparam radius_type Inferred type of radius the sphere is supposed to have.
 * @param x_coord X-coordinate of the point for that laying inside the sphere should be checked.
 * @param y_coord Y-coordinate of the point for that laying inside the sphere should be checked.
 * @param z_coord Z-coordinate of the point for that laying inside the sphere should be checked.
 * @param radius Radius of the filled sphere.
 * @param center_x X-coordinate of the sphere center.
 * @param center_y Y-coordinate of the sphere center.
 * @param center_z Z-coordinate of the sphere center.
 * @return True, if point lays inside the sphere; false, if outside.
 */
template <typename coord_type, typename radius_type>
bool inside_sphere(coord_type x_coord, coord_type y_coord, coord_type z_coord, radius_type radius, float center_x=0, float center_y=0, float center_z=0)
{
	return (x_coord - center_x) * (x_coord - center_x) + (y_coord - center_y) * (y_coord - center_y) + (z_coord - center_z) * (z_coord - center_z) <= radius * radius;
}

/**@brief Creates a filled sphere of given radius an given center coordinates on a 3D OpenFPM grid.
 *
 * @details The sphere is represented in one of the grid properties via the following indicator function:
 * @f[ \phi_{\text{indicator}} = \begin{cases}
 *    +1 & \text{point lies inside the sphere} \\
 *    -1 & \text{point lies outside the sphere} \\
 * \end{cases} @f]
 *
 * @tparam Phi_0 Index of grid property that should store the indicator function representing the filled sphere.
 * @tparam grid_type Inferred type of the input grid.
 * @tparam radius_type Inferred type of radius the sphere is supposed to have.
 * @param grid Input OpenFPM grid.
 * @param radius Radius of the filled sphere.
 * @param center_x X-coordinate of the sphere center.
 * @param center_y Y-coordinate of the sphere center.
 * @param center_z Z-coordinate of the sphere center.
 */
template <size_t Phi_0, typename grid_type, typename radius_type>
void init_grid_with_sphere(grid_type & grid, radius_type radius, float center_x=0, float center_y=0, float center_z=0)
{
	const size_t x = 0;
	const size_t y = 1;
	const size_t z = 2;
	// assign pixel values to domain. For each pixel get factor_refinement number of grid points with corresponding value
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		auto gkey = grid.getGKey(key);
		auto spacing = grid.getSpacing();
		auto i = gkey.get(x);
		auto j = gkey.get(y);
		auto k = gkey.get(z);
		double x_coord = i * spacing[x];
		double y_coord = j * spacing[y];
		double z_coord = k * spacing[z];
		
		if (inside_sphere(x_coord, y_coord, z_coord, radius, center_x, center_y, center_z))
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


#endif //ACCURACY_TESTS_SPHERE_HPP
