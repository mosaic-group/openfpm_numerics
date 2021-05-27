//
// Created by jstark on 27.05.21.
//

#ifndef OPENFPM_PDATA_ANALYTICALSDF_HPP
#define OPENFPM_PDATA_ANALYTICALSDF_HPP

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
 * @tparam point_type Template type of OpenFPM Point<dimension, type>.
 * @tparam radius_type Template type of radius.
 * @param coords Point_type coordinates of point.
 * @param radius Radius_type radius of the sphere.
 * @param center_x Double x-coordinate of sphere center.
 * @param center_y Double y-coordinate of sphere center.
 * @param center_z Double z-coordinate of sphere center.
 * @return Double variable that contains the exact solution for the signed distance function of a given point in a
 *         sphere of given radius, where the SDF has positive values inside and negative values outside the sphere.
 */
template <typename point_type, typename radius_type, typename center_type>
double get_analytic_sdf_sphere(point_type coords, radius_type radius,
                               center_type center_x=0, center_type center_y=0, center_type center_z=0)
{
	const double X = coords.get(0), Y = coords.get(1), Z = coords.get(2);
	return (radius -
	sqrt((X - center_x) * (X - center_x)
	+ (Y - center_y) * (Y - center_y)
	+ (Z - center_z) * (Z - center_z)));
}

/**@brief Initializes the exact solution of the signed distance function of a sphere on an OpenFPM grid.
 *
 * @details Solves the exact SDF for each grid nodes and writes the solution to a given property.
 *
 * @tparam SDF_exact Index of property where the anaytical solution for the SDF should be written to.
 * @tparam grid_type Inferred type of the input grid.
 * @tparam radius_type Inferred type of radius the sphere is supposed to have.
 * @param grid Input OpenFPM grid.
 * @param radius Radius of the filled sphere.
 * @param center_x Double x-coordinate of sphere center.
 * @param center_y Double y-coordinate of sphere center.
 * @param center_z Double z-coordinate of sphere center.

 */
template <size_t SDF_exact, typename grid_type, typename radius_type, typename center_type>
void init_analytic_sdf_sphere(grid_type & grid, radius_type radius, center_type center_x=0, center_type center_y=0,
                              center_type center_z=0)
{
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		Point<grid_type::dims, double> coords = grid.getPos(key);
		grid.template getProp<SDF_exact>(key) = get_analytic_sdf_sphere(coords, radius, center_x,
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
 * @tparam point_type Template type of OpenFPM Point<dimension, type>.
 * @tparam radius_type Template type of radius.
 * @param coords Point_type coordinates of point.
 * @param radius Radius_type radius of the sphere.
 * @param center_x Double x-coordinate of sphere center.
 * @param center_y Double y-coordinate of sphere center.
 * @param center_z Double z-coordinate of sphere center.
 * @return Double variable that contains the exact solution for the signed distance function of a given point in a
 *         sphere of given radius, where the SDF has positive values inside and negative values outside the sphere.
 */
template <typename point_type, typename radius_type, typename center_type>
double get_analytic_sdf_circle(point_type coords, radius_type radius,
                               center_type center_x=0, center_type center_y=0)
{
	const double X = coords.get(0), Y = coords.get(1);
	return (radius -
			sqrt((X - center_x) * (X - center_x)
					     + (Y - center_y) * (Y - center_y)));
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
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		Point<grid_type::dims, double> coords = grid.getPos(key);
		grid.template getProp<SDF_exact>(key) = get_analytic_sdf_sphere(coords, radius, center_x,
		                                                                center_y, center_z);
		++dom;
	}
}

#endif //ANALYTICAL_SDF_HPP


#endif //OPENFPM_PDATA_ANALYTICALSDF_HPP
