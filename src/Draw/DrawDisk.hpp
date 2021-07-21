//
// Created by jstark on 2020-08-18.
//
/**
 * @file Circle.hpp
 *
 * @brief Header file containing functions for creating a filled 2D disk of defined radius.
 *
 * @author Justina Stark
 * @date August 2020
 */
#ifndef ACCURACY_TESTS_CIRCLE_HPP
#define ACCURACY_TESTS_CIRCLE_HPP

#include <iostream>

/**@brief Checks if a point lays inside a disk of given radius and given center coordinates.
 *
 * @tparam point_type Template type of OpenFPM Point<dimension, type>.
 * @tparam radius_type Template type of radius.
 * @param coords Point_type coordinates of point.
 * @param radius Radius_type radius of the disk.
 * @param center_x Double x-coordinate of disk center.
 * @param center_y Double y-coordinate of disk center.
 * @return
 */
template <typename point_type, typename radius_type>
bool inside_disk(point_type coords, radius_type radius, double center_x=0, double center_y=0)
{
        const double EPSILON = std::numeric_limits<double>::epsilon();
        const double X = coords.get(0), Y = coords.get(1);
        return (X - center_x) * (X - center_x)
        + (Y - center_y) * (Y - center_y)
        <= radius * radius + EPSILON;
}

/**@brief Creates a filled disk of given radius an given center coordinates on a 3D OpenFPM grid.
 *
 * @details The disk is represented in one of the grid properties via the following indicator function:
 * @f[ \phi_{\text{indicator}} = \begin{cases}
 *    +1 & \text{point lies inside the disk} \\
 *    -1 & \text{point lies outside the disk} \\
 * \end{cases} @f]
 *
 * @tparam Phi_0 Index of grid property that should store the indicator function representing the filled disk.
 * @tparam grid_type Inferred type of the input grid.
 * @tparam radius_type Inferred type of radius the disk is supposed to have.
 * @param grid Input OpenFPM grid.
 * @param radius Radius of the filled disk.
 * @param center_x Double x-coordinate of disk center.
 * @param center_y Double y-coordinate of disk center.
 */
template <size_t Phi_0, typename grid_type, typename radius_type>
void init_grid_with_disk(grid_type & grid, radius_type radius, double center_x=0, double center_y=0)
{
        // assign pixel values to domain. For each pixel get factor_refinement number of grid points with corresponding value
        auto dom = grid.getDomainIterator();
        while(dom.isNext())
        {
                auto key = dom.get();

                Point<grid_type::dims, double> coords = grid.getPos(key);

                if (inside_disk(coords, radius, center_x, center_y))
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
