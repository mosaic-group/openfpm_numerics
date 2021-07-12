//
// Created by jstark on 2020-05-17.
//
/**
 * @file LNorms.hpp
 *
 * @brief Header file containing functions for computing the error and the L_2 / L_infinity norm.
 *
 * @details The L_2 and L_inf norm are useful for verification by convergence testings / plots.
 *
 * @author Justina Stark
 * @date May 2020
 */
#ifndef ACCURACY_TESTS_LNORMS_HPP
#define ACCURACY_TESTS_LNORMS_HPP

#include <iostream>
#include <typeinfo>
#include <cmath>
#include <cstdio>
#include <stdlib.h>

// Include OpenFPM header files
#include "Vector/vector_dist.hpp"
#include "Grid/grid_dist_id.hpp"
#include "util/PathsAndFiles.hpp"

/**@brief Structure that bundles the two variables for L_2 and L_infinity norm.
 */
struct L_norms
{
	double l2;   ///< Double variable that is supposed to contain the L_2 norm.
	double linf; ///< Double variable that is supposed to contain the L_infinity norm.
};
/**@brief At each grid node, the absolute error is computed and stored in another property.
 *
 * @details The absolute error is computed as:
 *
 * @f[ Error_{abs} = | numerical solution - analytical solution | @f]
 *
 * @tparam PropNumeric Index of grid property that contains the numerical value.
 * @tparam PropAnalytic Index of grid property that contains the analytical (exact) value.
 * @tparam Error Index of grid property where the computed error should be written to.
 * @tparam gridtype Inferred type of the input grid.
 * @param grid Input OpenFPM grid. Can be of any dimension.
 */
template <size_t PropNumeric, size_t PropAnalytic, size_t Error, typename gridtype>
void get_absolute_error(gridtype & grid)
{
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		grid.template getProp<Error> (key) = abs(grid.template getProp<PropAnalytic> (key) - (grid.template getProp<PropNumeric> (key)));
		++dom;
	}
}
/**@brief At each grid node, the relative error is computed and stored in another property.
 *
 * @details The relative error is computed as:
 *
 * @f[ Error_{rel} = | 1 - \frac{numerical solution}{analytical solution} | @f]
 *
 * @tparam PropNumeric Index of grid property that contains the numerical value.
 * @tparam PropAnalytic Index of grid property that contains the analytical (exact) value.
 * @tparam Error Index of grid property where the computed error should be written to.
 * @tparam gridtype Inferred type of the input grid.
 * @param grid Input OpenFPM grid. Can be of any dimension.
 */
template <size_t PropNumeric, size_t PropAnalytic, size_t Error, typename gridtype>
void get_relative_error(gridtype & grid)
{
	const double EPSILON = std::numeric_limits<double>::epsilon();
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		grid.template getProp<Error> (key) = abs( 1 - (grid.template getProp<PropNumeric> (key) / (grid.template
				getProp<PropAnalytic> (key))) );
		++dom;
	}
}
/**@brief Computes the L_2 and L_infinity norm on the basis of the precomputed error on a grid.
 *
 * @tparam Error Index of grid property that contains the error.
 * @tparam gridtype Inferred type of the input grid.
 * @param grid Input OpenFPM grid. Can be of any dimension.
 * @return Object of type L_norms that contains #L_norms::l2 and #L_norms::linf.
 */
template <size_t Error, typename gridtype>
L_norms get_l_norms_grid(gridtype & grid)
{
	double maxError = 0;
	double sumErrorSq = 0;
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		sumErrorSq += grid.template getProp<Error> (key) * grid.template getProp<Error> (key);
		if (grid.template getProp<Error> (key) > maxError) maxError = grid.template getProp<Error> (key); // update maxError
		++dom;
	}
	auto &v_cl = create_vcluster();
	v_cl.sum(sumErrorSq);
	v_cl.max(maxError);
	v_cl.execute();
	double l2 = sqrt( sumErrorSq / (double)grid.size());
	double linf = maxError;
	return {l2, linf};
}
/**@brief Computes the L_2 and L_infinity norm on the basis of the precomputed error on a particle vector.
 *
 * @tparam Error Index of grid property that contains the error.
 * @tparam vectortype Inferred type of the input particle vector.
 * @param vd Input particle vector.
 * @return Object of type L_norms that contains #L_norms::l2 and #L_norms::linf.
 */
template <size_t Error, typename vectortype>
L_norms get_l_norms_vector(vectortype & vd)
{
	double maxError = 0;
	double sumErrorSq = 0;
	auto dom = vd.getDomainIterator();
	double count = 0;
	while(dom.isNext())
	{
		auto key = dom.get();
		sumErrorSq += vd.template getProp<Error> (key) * vd.template getProp<Error> (key);
		if (vd.template getProp<Error> (key) > maxError) maxError = vd.template getProp<Error> (key); // update maxError
		++dom;
		++count;
	}
	auto &v_cl = create_vcluster();
	v_cl.sum(sumErrorSq);
	v_cl.sum(count);
	v_cl.max(maxError);
	v_cl.execute();
	double l2 = sqrt( sumErrorSq / (double)count);
	double linf = maxError;
	return {l2, linf};
}
/**@brief Converts value into string maintaining a desired precision.
 *
 * @tparam T Template type of vlaue.
 * @param myValue Value of type T.
 * @param n Number of digits after the point the string of the value should have
 * @return String containing myValue with precision n.
 */
template <typename T>
std::string to_string_with_precision(const T myValue, const size_t n = 6)
{
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << myValue;
	return out.str();
}
/**@brief Writes the N (number of grid points on a square grid) and L-norms as strings to a csv-file.
 *
 * @param N Size_t variable that contains the grid size in number of grid points in one dimension for an NxN(xN) grid
 * @param l_norms Object of type L_norms that contains #L_norms::l2 and #L_norms::linf.
 * @param filename Std::string containing the name of the csv file (without the .csv) to which the l-norms should be
 *                 written to.
 * @param path_output Std::string containing the path where the output csv file should be saved.
 */
template <typename T>
static void write_lnorms_to_file(T N, L_norms l_norms, std::string filename, std::string path_output)
{
	auto &v_cl = create_vcluster();
	if (v_cl.rank() == 0)
	{
		std::string path_output_lnorm = path_output + "/" + filename + ".csv";
		create_file_if_not_exist(path_output_lnorm);
		
		std::ofstream l_out;
		l_out.open(path_output_lnorm, std::ios_base::app); // append instead of overwrite
		l_out << to_string_with_precision(N, 16) << ',' << to_string_with_precision(l_norms.l2, 16)
		<< ',' << to_string_with_precision(l_norms.linf, 16) << std::endl;
		l_out.close();
	}
}




#endif //ACCURACY_TESTS_LNORMS_HPP

