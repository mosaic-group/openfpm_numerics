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

#include "level_set/redistancing_Sussman/HelpFunctions.hpp" // for printing to_string_with_precision

/**@brief At each grid node, the absolute error is computed and stored in another property.
 *
 * @details The absolute error is computed as:
 *
 * @f[ Error_{abs} = | numerical solution - analytical solution | @f]
 *
 * @tparam PropNumeric Index of grid property that contains the numerical value.
 * @tparam PropAnalytic Index of grid property that contains the analytical (exact) value.
 * @tparam Error Index of grid property where the computed error should be written to.
 * @tparam gridtype Template type of the input grid.
 * @param grid Input OpenFPM grid. Can be of any dimension.
 */
template <size_t PropNumeric, size_t PropAnalytic, size_t Error, typename gridtype>
void get_absolute_error(gridtype & grid)
{
	auto dom = grid.getDomainIterator();
	typedef typename std::remove_const_t<std::remove_reference_t<decltype(
	grid.template get<PropNumeric>(dom.get()))>> numeric_type;
	typedef typename std::remove_const_t<std::remove_reference_t<decltype(
	grid.template get<PropAnalytic>(dom.get()))>> analytic_type;
	
	if(!(std::is_same<numeric_type, analytic_type>::value))
	{
		std::cout << "Type of numerical and analytical solution must be the same! Aborting..." << std::endl;
		abort();
	}	while(dom.isNext())
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
 * @tparam gridtype Template type of the input grid.
 * @param grid Input OpenFPM grid. Can be of any dimension.
 */
template <size_t PropNumeric, size_t PropAnalytic, size_t Error, typename gridtype>
void get_relative_error(gridtype & grid)
{
	auto dom = grid.getDomainIterator();
	typedef typename std::remove_const_t<std::remove_reference_t<decltype(
	grid.template get<PropNumeric>(dom.get()))>> numeric_type;
	typedef typename std::remove_const_t<std::remove_reference_t<decltype(
	grid.template get<PropAnalytic>(dom.get()))>> analytic_type;
	
	if(!(std::is_same<numeric_type, analytic_type>::value))
	{
		std::cout << "Type of numerical and analytical solution must be the same! Aborting..." << std::endl;
		abort();
	}
	
	while(dom.isNext())
	{
		auto key = dom.get();
		grid.template getProp<Error> (key) = abs( 1 - (grid.template getProp<PropNumeric> (key) / (grid.template
				getProp<PropAnalytic> (key))) );
		++dom;
	}
}

/**@brief Class for computing the l2/l_infinity norm for distributed grids and vectors based on given errors.
 *
 * @tparam lnorm_type Return type for l-norm.
 */
template <typename lnorm_type>
class LNorms
{
public:
	LNorms() = default;
	// Member variables
	lnorm_type l2; // L2 norm
	lnorm_type linf; // L_infinity norm
	
	// Member functions
	/**@brief Computes the L_2 and L_infinity norm on the basis of the precomputed error on a grid.
	 *
	 * @tparam Error Index of grid property that contains the error.
	 * @tparam gridtype Template type of the input grid.
	 * @param grid Input OpenFPM grid. Can be of any dimension.
	 */
	template <size_t Error, typename gridtype>
	void get_l_norms_grid(gridtype & grid)
	{
		auto dom = grid.getDomainIterator();
		typedef typename std::remove_const_t<std::remove_reference_t<decltype(
		grid.template get<Error>(dom.get()))>> error_type;
		
		error_type maxError = 0;
		error_type sumErrorSq = 0;
		
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
		l2 = (lnorm_type) sqrt( sumErrorSq / (error_type)grid.size());
		linf = (lnorm_type) maxError;
	}
	/**@brief Computes the L_2 and L_infinity norm on the basis of the precomputed error on a particle vector.
	 *
	 * @tparam Error Index of grid property that contains the error.
	 * @tparam vectortype Template type of the input particle vector.
	 * @param vd Input particle vector.
	 */
	template <size_t Error, typename vectortype>
	void get_l_norms_vector(vectortype & vd)
	{
		auto dom = vd.getDomainIterator();
		typedef typename std::remove_const_t<std::remove_reference_t<decltype(
		vd.template getProp<Error>(dom.get()))>> error_type;
		
		error_type maxError = 0;
		error_type sumErrorSq = 0;
		int count = 0;
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
		l2 = (lnorm_type) sqrt( sumErrorSq / (error_type)count);
		linf = (lnorm_type) maxError;
	}
	/**@brief Writes the N (number of grid points on a square grid) and L-norms as strings to a csv-file.
	 *
	 * @param N Size_t variable that contains the grid size in number of grid points in one dimension for an NxN(xN) grid
	 * @param precision Precision in number of digits after the points for writing of the l-norms to file.
	 * @param filename Std::string containing the name of the csv file (without the .csv) to which the l-norms should be
	 *                 written to.
	 * @param path_output Std::string containing the path where the output csv file should be saved.
	 */
	 void write_to_file(const size_t N,
						const int precision,
						const std::string & filename,
						const std::string & path_output)
	{
		auto &v_cl = create_vcluster();
		if (v_cl.rank() == 0)
		{
			std::string path_output_lnorm = path_output + "/" + filename + ".csv";
			create_file_if_not_exist(path_output_lnorm);
			
			std::ofstream l_out;
			l_out.open(path_output_lnorm, std::ios_base::app); // append instead of overwrite
			l_out << std::to_string(N)
					<< ',' << to_string_with_precision(l2, precision)
					<< ',' << to_string_with_precision(linf, precision) << std::endl;
			l_out.close();
		}
	}

private:


};









#endif //ACCURACY_TESTS_LNORMS_HPP

