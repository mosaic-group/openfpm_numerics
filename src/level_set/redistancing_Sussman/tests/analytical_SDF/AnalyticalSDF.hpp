//
// Created by jstark on 2020-10-05.
//
/**
 * @file AnalyticalSDF.hpp
 *
 * @brief Header file containing functions that compute the analytical solution of the signed distance function (SDF)
 *        for a 3D sphere and a 2D circle of defined radius on a grid.
 *
 * @author Justina Stark
 * @date October 2020
 */
#ifndef ANALYTICAL_SDF_HPP
#define ANALYTICAL_SDF_HPP

#include <iostream>
/**@brief Computes the analytical signed distance function of a sphere for a given point in space.
 *
 * @details At the center of the sphere, \a &phi;_SDF_analytic = \a Radius. Moving outwards from the center on, the
 * value for the SDF decreases steadily, eventually becomes 0 at the sphere surface and negative beyond the surface. The
 * analytic SDF for a sphere of radius \a R and center (\a a, \a b, \a c) is thus:
 *
 * @f[ \phi_{SDF}(x, y, z) = R - \sqrt{((x-a)^2 + (y-b)^2 + (z-c)^2)} @f]
 *
 * @tparam coord_type Inferred type of the point coordinates x, y and z.
 * @tparam radius_type Inferred type of radius the sphere is supposed to have.
 * @tparam center_type Inferred type of the center coordinates.
 * @param x_coord X-coordinate of the point for that laying inside the sphere should be checked.
 * @param y_coord Y-coordinate of the point for that laying inside the sphere should be checked.
 * @param z_coord Z-coordinate of the point for that laying inside the sphere should be checked.
 * @param radius Radius of the sphere.
 * @param center_x X-coordinate of the sphere center.
 * @param center_y Y-coordinate of the sphere center.
 * @param center_z Z-coordinate of the sphere center.
 * @return Double variable that contains the exact solution for the signed distance function of a given point in a
 *         sphere of given radius, where the SDF has positive values inside and negative values outside the sphere.
 */
template <typename coord_type, typename radius_type, typename center_type>
double get_analytic_sdf_sphere(coord_type x_coord, coord_type y_coord, coord_type z_coord, radius_type radius,
		center_type center_x=0, center_type center_y=0, center_type center_z=0)
{
	return (radius - sqrt((x_coord - center_x) * (x_coord  - center_x) + (y_coord  - center_y) * (y_coord  -
	center_y) + (z_coord - center_z) * (z_coord - center_z)));
}

/**@brief Initializes the exact solution of the signed distance function of a sphere on an OpenFPM grid.
 *
 * @details Solves the exact SDF for each grid nodes and writes the solution to a given property.
 *
 * @tparam SDF_exact Index of property where the anaytical solution for the SDF should be written to.
 * @tparam grid_type Inferred type of the input grid.
 * @tparam radius_type Inferred type of radius the sphere is supposed to have.
 * @tparam center_type Inferred type of the center coordinates.
 * @param grid Input OpenFPM grid.
 * @param radius Radius of the filled sphere.
 * @param center_x X-coordinate of the sphere center.
 * @param center_y Y-coordinate of the sphere center.
 * @param center_z Z-coordinate of the sphere center.
 */
template <size_t SDF_exact, typename grid_type, typename radius_type, typename center_type>
void init_analytic_sdf_sphere(grid_type & grid, radius_type radius, center_type center_x=0, center_type center_y=0,
		center_type center_z=0)
{
	// some indices
	const size_t x                      = 0;
	const size_t y                      = 1;
	const size_t z                      = 2;
	
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto spacing = grid.getSpacing();
		auto key = dom.get();
		auto gkey = grid.getGKey(key);
		auto i = gkey.get(x);
		auto j = gkey.get(y);
		auto k = gkey.get(z);
		double x_coord = i * spacing[x];
		double y_coord = j * spacing[y];
		double z_coord = k * spacing[z];
		grid.template getProp<SDF_exact>(key) = get_analytic_sdf_sphere(x_coord, y_coord, z_coord, radius, center_x,
				center_y, center_z);
		++dom;
	}
}

/**@brief Computes the analytical signed distance function of a circle for a given point in space.
 *
 * @details At the center of the circle, \a &phi;_SDF_analytic = \a Radius. Moving outwards from the center on, the
 * value for the SDF decreases steadily, eventually becomes 0 at the circle surface and negative beyond the surface. The
 * analytic SDF for a circle of radius \a R and center (\a a, \a b, \a c) is thus:
 *
 * @f[ \phi_{SDF}(x, y, z) = R - \sqrt{((x-a)^2 + (y-b)^2 + (z-c)^2)} @f]
 *
 * @tparam coord_type Inferred type of the point coordinates x, y and z.
 * @tparam radius_type Inferred type of radius the circle is supposed to have.
 * @tparam center_type Inferred type of the center coordinates.
 * @param x_coord X-coordinate of the point for that laying inside the circle should be checked.
 * @param y_coord Y-coordinate of the point for that laying inside the circle should be checked.
 * @param radius Radius of the circle.
 * @param center_x X-coordinate of the circle center.
 * @param center_y Y-coordinate of the circle center.
 * @return Double variable that contains the exact solution for the signed distance function of a given point in a
 *         circle of given radius, where the SDF has positive values inside and negative values outside the circle.
 */
template <typename coord_type, typename radius_type, typename center_type>
double get_analytic_sdf_circle(coord_type x_coord, coord_type y_coord, radius_type radius,
                               center_type center_x=0, center_type center_y=0)
{
	return (radius - sqrt((x_coord - center_x) * (x_coord  - center_x) + (y_coord  - center_y) * (y_coord  -
	center_y)));
}

/**@brief Initializes the exact solution of the signed distance function of a circle on an OpenFPM grid.
 *
 * @details Solves the exact SDF for each grid nodes and writes the solution to a given property.
 *
 * @tparam SDF_exact Index of property where the anaytical solution for the SDF should be written to.
 * @tparam grid_type Inferred type of the input grid.
 * @tparam radius_type Inferred type of radius the circle is supposed to have.
 * @tparam center_type Inferred type of the center coordinates.
 * @param grid Input OpenFPM grid.
 * @param radius Radius of the filled circle.
 * @param center_x X-coordinate of the circle center.
 * @param center_y Y-coordinate of the circle center.
 */
template <size_t SDF_exact, typename grid_type, typename radius_type, typename center_type>
void init_analytic_sdf_circle(grid_type & grid, radius_type radius, center_type center_x=0, center_type center_y=0)
{
	// some indices
	const size_t x                      = 0;
	const size_t y                      = 1;
	
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto spacing = grid.getSpacing();
		auto key = dom.get();
		auto gkey = grid.getGKey(key);
		auto i = gkey.get(x);
		auto j = gkey.get(y);
		double x_coord = i * spacing[x];
		double y_coord = j * spacing[y];
		grid.template getProp<SDF_exact>(key) = get_analytic_sdf_circle(x_coord, y_coord, radius, center_x, center_y);
		++dom;
	}
}

#endif //ANALYTICAL_SDF_HPP