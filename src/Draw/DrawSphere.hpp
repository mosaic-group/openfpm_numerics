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
 * @tparam point_type Template type of OpenFPM Point<dimension, type>.
 * @tparam radius_type Template type of radius.
 * @param coords Point_type coordinates of point.
 * @param radius Radius_type radius of the sphere.
 * @param center_x Double x-coordinate of sphere center.
 * @param center_y Double y-coordinate of sphere center.
 * @param center_z Double z-coordinate of sphere center.
 * @return
 */
template <typename point_type, typename radius_type>
bool inside_sphere(point_type coords, radius_type radius, double center_x=0, double center_y=0, double center_z=0)
{
	typedef typename std::remove_const_t<std::remove_reference_t<decltype(coords.get(0))>> space_type;
	if(!(std::is_same<space_type, radius_type>::value))
	{
		std::cout << "Radius type and space type of grid must be the same! Aborting..." << std::endl;
		abort();
	}
	const space_type EPSILON = std::numeric_limits<space_type>::epsilon();
	const space_type X = coords.get(0), Y = coords.get(1), Z = coords.get(2);
	return (X - center_x) * (X - center_x)
	+ (Y - center_y) * (Y - center_y)
	+ (Z - center_z) * (Z - center_z)
	<= radius * radius + EPSILON;
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
 * @param center_x Double x-coordinate of sphere center.
 * @param center_y Double y-coordinate of sphere center.
 * @param center_z Double z-coordinate of sphere center.
 */
template <size_t Phi_0, typename grid_type, typename radius_type>
void init_grid_with_sphere(grid_type & grid, radius_type radius, double center_x=0, double center_y=0, double center_z=0)
{
	if(!(std::is_same<typename grid_type::stype, radius_type>::value))
	{
		std::cout << "Radius type and space type of grid must be the same! Aborting..." << std::endl;
		abort();
	}
	// assign pixel values to domain. For each pixel get factor_refinement number of grid points with corresponding value
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
	
		Point<grid_type::dims, typename grid_type::stype> coords = grid.getPos(key);

		if (inside_sphere(coords, radius, center_x, center_y, center_z))
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
