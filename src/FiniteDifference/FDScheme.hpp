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
	// Padding
	Padding<Sys_eqs::dims> pd;

	// Sparse Matrix
	openfpm::vector<triplet<typename Sys_eqs::stype>> trpl;

	openfpm::vector<typename Sys_eqs::stype> b;

	// Domain Grid informations
	const grid_sm<Sys_eqs::dims,void> & gs;

	// mapping grid
	grid_dist_id<Sys_eqs::dims,Sys_eqs::stype,scalar<size_t>,Sys_eqs::grid_type::decomposition> g_map;

	// row of the matrix
	size_t row;

	// row on b
	size_t row_b;

	/*! \brief Check if the Matrix is consistent
	 *
	 */
	void consistency()
	{
		// A and B must have the same rows
		if (row != row_b)
			std::cerr << "Error " << __FILE__ << ":" << __LINE__ << "the term B and the Matrix A for Ax=B must contain the same number of rows\n";

		// Indicate all the non zero rows
		openfpm::vector<bool> nz_rows;
		nz_rows.resize(row_b);

		for (size_t i = 0 ; i < trpl.size() ; i++)
		{
			nz_rows.get(trpl.get(i).i) = true;
		}

		// Indicate all the non zero colums
		openfpm::vector<bool> nz_cols;
		for (size_t i = 0 ; i < trpl.size() ; i++)
		{
			nz_cols.get(trpl.get(i).j) = true;
		}

		// all the rows must have a non zero element
		for (size_t i = 0 ; i < nz_rows.size() ; i++)
		{
			if (nz_rows.get(i) == false)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Ill posed matrix not all the row are filled\n";
		}

		// all the colums must have a non zero element
		for (size_t i = 0 ; i < nz_cols.size() ; i++)
		{
			if (nz_cols.get(i) == false)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Ill posed matrix not all the row are filled\n";
		}
	}

public:

	/*! \brief Constructor
	 *
	 * \param pd Padding
	 * \param gs grid infos where Finite differences work
	 *
	 */
	FDScheme(Padding<Sys_eqs::dims> & pd, const grid_sm<Sys_eqs::dims,void> & gs)
	:pd(pd),gs(gs)
	{
		// Counter
		size_t cnt = 0;

		// Resize the b term
		b.resize(gs.size());

		// Create the re-mapping-grid
		auto it = g_map.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			g_map.get(key) = cnt;

			++cnt;
			++it;
		}


	}

	/*! \brief Impose an operator
	 *
	 *
	 *
	 */
	template<typename T> void imposeA(const T & op , grid_dist_iterator_sub<Sys_eqs::dims,typename Sys_eqs::b_grid::d_grid> it)
	{
		std::unordered_map<long int,float> cols;

		// iterate all the grid points
		while (it.isNext())
		{
			// get the position
			auto key = it.get();

			// convert into global coordinate the position
			auto keyg = it.getGKey(key);

			// Calculate the non-zero colums
			T::value(keyg,gs,cols,1.0);

			// create the triplet

			for ( auto it = cols.begin(); it != cols.end(); ++it )
			{
				trpl.add();
				trpl.last().i = row;
				trpl.last().j = it->first;
				trpl.last().value = it->second;

				std::cout << "(" << trpl.last().i << "," << trpl.last().j << "," << trpl.last().value << ")" << "\n";
			}

			cols.clear();
			std::cout << "\n";

			++row;
			++it;
		}
	}

	/*! \brief Impose an operator
	 *
	 *
	 *
	 */
	void imposeB(typename Sys_eqs::stype num , grid_dist_iterator_sub<Sys_eqs::dims,typename Sys_eqs::b_grid::d_grid> it)
	{
		std::unordered_map<long int,float> cols;

		// iterate all the grid points
		while (it.isNext())
		{

			b.get(row_b) = num;

			++row_b;
			++it;
		}
	}

	/*! \brief produce the Matrix
	 *
	 *  \tparam Syst_eq System of equations, or a single equation
	 *
	 */
	template<typename Grid> SparseMatrix<typename Sys_eqs::stype> getMatrix()
	{
#ifdef SE_CLASS1
		consistency();
#endif


	}
};

#define EQS_FIELDS 0
#define EQS_SPACE 1


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_FDSCHEME_HPP_ */
