/*
 * FiniteDifferences.hpp
 *
 *  Created on: Sep 17, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_

#include "Grid/grid_dist_id.hpp"
#include "Matrix/Matrix.hpp"
#include "eq.hpp"

/*! \brief Finite Differences
 *
 * \tparam dim Dimensionality of the finite differences scheme
 *
 */

template<typename Sys_eqs>
class FDScheme
{
	openfpm::vector<triplet<typename Sys_eqs::stype>> trpl;

	/*! \brief Impose an operator
	 *
	 *
	 *
	 */
/*	template<typename T> void impose(T & op, grid_key_dx_dist_iterator_sub & it)
	{
		size_t lin = 0;

		std::unordered_map<long int,float> cols;

		// iterate all the grid points
		while (it.isNext())
		{
			// get the position
			auto key = it.get();

			// convert into global coordinate the position
			auto keyg = g.getGKey(key);

			// Convert local key into global key
			T::value(op,cols);

			// create the triplet

			for ( auto it = cols.begin(); it != cols.end(); ++it )
			{
				trpl.add();
				trpl.last().i = lin;
				trpl.last().j = it->first;
				trpl.last().value = it->second;
			}

			std::cout << " " << it->first << ":" << it->second;

			cols.clear();

			++loc_cnt;
			++it;
		}
	}*/

	/*! \brief produce the Matrix
	 *
	 *  \tparam Syst_eq System of equations, or a single equation
	 *
	 */
	template<typename Grid> SparseMatrix<typename Sys_eqs::stype> getMatrix(const Grid & g, const ConstField(& fld)[Sys_eqs::num_cfields::value] )
	{
		// iterate across the full domain

		auto it = g.getDomainIterator();
		auto g_sm = g.getGrid();

		Sys_eqs eqs();

		// iterate all the grid points
		while (it.isNext())
		{
			// get the position
			auto key = it.get();

			// convert into global coordinate the position
			auto keyg = g.getGKey(key);

			// Convert local key into global key
			size_t row = g_sm.LidId(keyg);

			eqs.value(keyg,trpl);

			++it;
		}
	}
};

#define EQS_FIELDS 0
#define EQS_SPACE 1


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_ */
