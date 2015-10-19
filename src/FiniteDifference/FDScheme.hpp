/*
 * FiniteDifferences.hpp
 *
 *  Created on: Sep 17, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_

#include "Grid/grid_dist_id.hpp"
#include "Grid/grid_dist_id_iterator_sub.hpp"
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

public:

	/*! \brief Impose an operator
	 *
	 *
	 *
	 */
	template<typename T> void impose(T & op, const grid_sm<Sys_eqs::dims,void> & gs , grid_dist_iterator_sub<Sys_eqs::dims,typename Sys_eqs::b_grid::d_grid> it)
	{
		size_t lin = 0;

		std::unordered_map<long int,float> cols;

		// iterate all the grid points
		while (it.isNext())
		{
			// get the position
			auto key = it.get();

			// convert into global coordinate the position
			auto keyg = it.getGKey(key);

			// Convert local key into global key
			T::value(keyg,gs,cols,1.0);

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

			++it;
		}
	}

	/*! \brief produce the Matrix
	 *
	 *  \tparam Syst_eq System of equations, or a single equation
	 *
	 */
	template<typename Grid> SparseMatrix<typename Sys_eqs::stype> getMatrix(const Grid & g)
	{
#ifdef SE_CLASS1
		if (Sys_eqs::num_cfields)
			std::cerr << __FILE__ << ":" << __LINE__ << " if you do not provide ConstantFields in getMatrix, the system should not contain any constant field, while" << Sys_eqs::num_cfields << "\n";
#endif

		const ConstField fld[Sys_eqs::num_cfields];

		getMatrix(g,fld);
	}

	/*! \brief produce the Matrix
	 *
	 *  \tparam Syst_eq System of equations, or a single equation
	 *
	 */
	template<typename Grid> SparseMatrix<typename Sys_eqs::stype> getMatrix(const Grid & g, const ConstField(& fld)[Sys_eqs::num_cfields] )
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
