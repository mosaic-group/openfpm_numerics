/*
 * grid_dist_testing.hpp
 *
 *  Created on: Oct 28, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_UTIL_GRID_DIST_TESTING_HPP_
#define OPENFPM_NUMERICS_SRC_UTIL_GRID_DIST_TESTING_HPP_

#include "Grid/grid_dist_key.hpp"
#include "Grid/grid_dist_key.hpp"
#include "Grid/map_grid.hpp"

template<unsigned int dim>
class grid_dist_testing
{
	grid_cpu<dim,aggregate<size_t>> grid_test;

public:

	/*! \brief It create a test map suitable for testing
	 *
	 * \param size of the grid with padding
	 *
	 */
	grid_dist_testing(const size_t (& g_sz)[dim])
	:grid_test(g_sz)
	{
		grid_test.setMemory();

		// We fill the testing map

		auto dom = grid_test.getIterator();
		const auto & gs = grid_test.getGrid();

		while (dom.isNext())
		{
			auto key = dom.get();

			grid_test.template get<0>(key) = gs.LinId(key);

			++dom;
		}
	}

	/*! \brief Get the element at position key
	 *
	 * \warning it is a stub object, it just return number in the underline grid
	 *
	 */
	template <unsigned int p> size_t get(const grid_dist_key_dx<dim> & key) const
	{
		if (p != 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " mapping grid is suppose to have only one scalar size_t \n";

		return grid_test.template get<0>(key.getKey());
	}

	/*! \brief Return the global key
	 *
	 * \warning it is a stub object, it just return the key
	 *
	 */
	grid_key_dx<dim> getGKey(const grid_dist_key_dx<dim> & key) const
	{
		return key.getKey();
	}
};


#endif /* OPENFPM_NUMERICS_SRC_UTIL_GRID_DIST_TESTING_HPP_ */
