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
#include "data_type/scalar.hpp"

/*! \brief Finite Differences
 *
 * This class is able to discreatize on a Matrix any operator. In order to create a consistent
 * Matrix it is required that each processor must contain a contiguos range on grid points without
 * holes. In order to ensure this, each processor produce a contiguos local labelling of its local
 * points. Each processor also add an offset equal to the number of local
 * points of the processors with id smaller than him, to produce a global and non overlapping
 * labelling. An example is shown in the figures down, here we have
 * a grid 8x6 divided across two processor each processor label locally its grid points
 *
 * \verbatim
 *
+--------------------------+
| 1   2   3   4| 1  2  3  4|
|              |           |
| 5   6   7   8| 5  6  7  8|
|              |           |
| 9  10  11  12| 9 10 11 12|
+--------------------------+
|13  14  15| 13 14 15 16 17|
|          |               |
|16  17  18| 18 19 20 21 22|
|          |               |
|19  20  21| 23 24 25 26 27|
+--------------------------+

 *
 *
 * \endverbatim
 *
 * To the local relabelling is added an offset to make the local id global and non overlapping
 *
 *
 * \verbatim
 *
+--------------------------+
| 1   2   3   4|23 24 25 26|
|              |           |
| 5   6   7   8|27 28 29 30|
|              |           |
| 9  10  12  13|31 32 33 34|
+--------------------------+
|14  15  16| 35 36 37 38 39|
|          |               |
|17  18  19| 40 41 42 43 44|
|          |               |
|20  21  22| 45 46 47 48 49|
+--------------------------+
 *
 *
 * \endverbatim
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

	typedef grid_dist_id<Sys_eqs::dims,typename Sys_eqs::stype,scalar<size_t>,typename Sys_eqs::b_grid::decomposition> g_map_type;

	// mapping grid
	g_map_type g_map;

	// row of the matrix
	size_t row;

	// row on b
	size_t row_b;

	// Grid points that has each processor
	openfpm::vector<size_t> pnt;

	// Each point in the grid has a global id, to decompose correctly the Matrix each processor contain a
	// contiguos range of global id, example processor 0 can have from 0 to 234 and processor 1 from 235 to 512
	// no processors can have holes in the sequence, this number indicate where the sequence start for this
	// processor
	size_t s_pnt;

	/* \brief calculate the mapping grid size with padding
	 *
	 * \param gs original grid size
	 *
	 * \return padded grid size
	 *
	 */
	inline const std::vector<size_t> padding( const size_t (& sz)[Sys_eqs::dims], Padding<Sys_eqs::dims> & pd)
	{
		std::vector<size_t> g_sz_pad(Sys_eqs::dims);

		for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
			g_sz_pad[i] = sz[i] + pd.getLow(i) + pd.getHigh(i);

		return g_sz_pad;
	}

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

	/* \brief Convert discrete ghost into continous ghost
	 *
	 * \param Ghost
	 *
	 */
	inline Ghost<Sys_eqs::dims,typename Sys_eqs::stype> convert_dg_cg(const Ghost<Sys_eqs::dims,typename Sys_eqs::stype> & g)
	{
		return g_map_type::convert_ghost(g,g_map.getCellDecomposer());
	}

public:

	/*! \brief Constructor
	 *
	 * \param pd Padding
	 * \param gs grid infos where Finite differences work
	 *
	 */
	FDScheme(Padding<Sys_eqs::dims> & pd, const Box<Sys_eqs::dims,typename Sys_eqs::stype> & domain, const grid_sm<Sys_eqs::dims,void> & gs, const typename Sys_eqs::b_grid::decomposition & dec, Vcluster & v_cl)
	:pd(pd),gs(gs),g_map(dec,padding(gs.getSize(),pd),domain,Ghost<Sys_eqs::dims,typename Sys_eqs::stype>(1)),row(0),row_b(0)
	{
		// Calculate the size of the local domain
		size_t sz = g_map.getLocalDomainSize();

		// Get the total size of the local grids on each processors
		v_cl.allGather(sz,pnt);
		v_cl.execute();
		s_pnt = 0;

		// calculate the starting point for this processor
		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
			s_pnt += pnt.get(i);

		// Calculate the starting point

		// Counter
		size_t cnt = 0;

		// Resize the b term
		b.resize(gs.size());

		// Create the re-mapping-grid
		auto it = g_map.getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			g_map.template get<0>(key) = cnt + s_pnt;

			++cnt;
			++it;
		}

		// sync the ghost

		g_map.template ghost_get<0>();
	}

	/*! \brief Impose an operator
	 *
	 *
	 *
	 */
	template<typename T> void imposeA(const T & op , grid_dist_iterator_sub<Sys_eqs::dims,typename Sys_eqs::b_grid::d_grid> it_d)
	{
		// key one
		grid_key_dx<Sys_eqs::dims> gk_one;
		gk_one.one();

		// add padding to start and stop
		grid_key_dx<Sys_eqs::dims> start = it_d.getStart() + pd.getKP1();
		grid_key_dx<Sys_eqs::dims> stop = it_d.getStop() + pd.getKP1();

		auto it = g_map.getSubDomainIterator(start,stop);

		std::unordered_map<long int,float> cols;

		// iterate all the grid points
		while (it.isNext())
		{
			// get the position
			auto key = it.get();

			// Calculate the non-zero colums
			T::value(g_map,key,gs,cols,1.0);

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
