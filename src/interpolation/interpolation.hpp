/*
 * mp4_interpolation.hpp
 *
 *  Created on: May 4, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_HPP_
#define OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_HPP_

#include "NN/CellList/MemFast.hpp"
#include "NN/CellList/CellList.hpp"
#include "Grid/grid_dist_key.hpp"

#define INTERPOLATION_ERROR_OBJECT std::runtime_error("Runtime interpolation error");

constexpr int inte_m2p = 0;
constexpr int inte_p2m = 1;

/*! \brief It store the offsets of the interpolation points
 *
 * \tparam n_ele number of interpolation points
 * \tparam T type in general an
 *
 */
template<unsigned int n_ele>
struct agg_arr
{
	//! offset of the interpolation points
	size_t ele[n_ele];
};

/*! \brief multiply the src  by coeff for several types T
 *
 * \tparam type T
 *
 */
template<typename T>
struct mul_inte
{
	/*! \brief multiply the src  by coeff for several types T
	 *
	 * \param result the result of the multiplication
	 * \param coeff coefficent to use for of the multiplication
	 * \param src source
	 *
	 */
	inline static void value(T & result, const T & coeff, const T & src)
	{
		result += coeff * src;
	}
};

/*! \brief multiply the src  by coeff for several types T
 *
 * \tparam type T
 *
 */
template<typename T, unsigned int N1>
struct mul_inte<T[N1]>
{
	/*! \brief multiply the src  by coeff for several types T
	 *
	 * \param result the result of the multiplication
	 * \param coeff coefficent to use for of the multiplication
	 * \param src source
	 *
	 */
	inline static void value(T (& result)[N1], const T & coeff, const T (& src)[N1])
	{
		for (size_t i = 0 ; i < N1 ; i++)
			result[i] += coeff * src[i];
	}
};

/*! \brief multiply the src  by coeff for several types T
 *
 * \tparam type T
 *
 */
template<typename T, unsigned int N1, unsigned int N2>
struct mul_inte<T[N1][N2]>
{
	/*! \brief multiply the src  by coeff for several types T
	 *
	 * \param result the result of the multiplication
	 * \param coeff coefficent to use for of the multiplication
	 * \param src source
	 *
	 */
	inline static void value(T (& result)[N1][N2], const T & coeff, const T (& src)[N1][N2])
	{
		for (size_t i = 0 ; i < N1 ; i++)
			for (size_t j = 0 ; j < N2 ; j++)
				result[i][j] += coeff * src[i][j];
	}
};

/*! \brief Class that select the operation to do differently if we are doing Mesh to particle (m2p) or particle to mesh (p2m)
 *
 * \tparam prp_g property of the grid to interpolate
 * \tparam prp_v property of the vector to interpolate
 * \param M2P or P2M
 *
 */
template<unsigned int np, unsigned int prp_g, unsigned int prp_v,unsigned int m2p_or_p2m>
struct inte_template
{
	/*! \brief Evaluate the interpolation
	 *
	 * \tparam np number of kernel points in one direction
	 * \tparam grid type of grid
	 * \tparam vector type of vector
	 * \tparam iterator type of the iterator
	 *
	 * \param gd grid for interpolation
	 * \param vd vector for interpolation
	 * \param k_dist grid key grid point for interpolation
	 * \param key_p particle for interpolation
	 * \param a_int interpolation coefficients pre-calculated
	 * \param key indicate which pre-calculated coefficient we have to use
	 *
	 */
	template<unsigned int np_a_int, typename grid, typename vector, typename iterator> inline static void value(grid & gd,
			                                                          vector & vd,
																	  const grid_dist_lin_dx & k_dist,
																	  iterator & key_p,
																	  typename vector::stype (& a_int)[np_a_int],
																	  const size_t & key)
	{
		mul_inte<typename std::remove_reference<decltype(gd.template get<prp_g>(k_dist))>::type>::value(gd.template get<prp_g>(k_dist),a_int[key],vd.template getProp<prp_v>(key_p));
	}
};

/*! \brief Class that select the operation to do differently if we are doing Mesh to particle (m2p) or particle to mesh (p2m)
 *
 * \tparam prp_g property of the grid to interpolate
 * \tparam prp_v property of the vector to interpolate
 * \param M2P or P2M
 *
 */
template<unsigned int np, unsigned int prp_g, unsigned int prp_v>
struct inte_template<np,prp_g,prp_v,inte_m2p>
{
	/*! \brief Evaluate the interpolation
	 *
	 * \tparam grid type of grid
	 * \tparam vector type of vector
	 * \tparam iterator type of the iterator
	 * \tparam key_type type of the key
	 *
	 * \param gd grid for interpolation
	 * \param vd vector for interpolation
	 * \param k_dist grid key grid point for interpolation
	 * \param key_p particle for interpolation
	 * \param a_int interpolation coefficients pre-calculated
	 * \param key indicate which pre-calculated coefficient we have to use
	 *
	 */
	template<unsigned int np_a_int, typename grid, typename vector, typename iterator> inline static void value(grid & gd,
			                                                          vector & vd,
																	  const grid_dist_lin_dx & k_dist,
																	  iterator & key_p,
																	  typename vector::stype (& a_int)[np_a_int],
																	  const size_t & key)
	{
		mul_inte<typename std::remove_reference<decltype(gd.template get<prp_g>(k_dist))>::type>::value(vd.template getProp<prp_v>(key_p),a_int[key],gd.template get<prp_g>(k_dist));
	}
};

/*! \brief Calculate aint
 *
 * This class store
 *
 */
template<unsigned int dims, typename vector, unsigned int np>
struct calculate_aint
{
	/*! \brief Calculate the coefficients of the interpolation a_int for one particles
	 *         having the 1D kernel values
	 *
	 * \param sz size of stencil for the interpolation
	 * \param a_int coefficients on the stencil points
	 * \param a coefficients of the kernel for each direction, consider
	 *          that for 3 dimensions the kernel is the multiplication
	 *          the 1D kernel on each direction. The "a" array store the
	 *          calculated coefficient of the 1D kernel on each direction
	 *
	 */
	static inline void value(size_t (& sz)[vector::dims],
			                 typename vector::stype a_int[openfpm::math::pow(np,vector::dims)],
							 typename vector::stype (& a)[vector::dims][np])
	{
		grid_sm<vector::dims,void> gs(sz);
		grid_key_dx_iterator<vector::dims> kit(gs);

		size_t s = 0;
		while (kit.isNext())
		{
			auto key = kit.get();

			a_int[s] = 1;

			for (size_t i = 0 ; i < vector::dims ; i++)
				a_int[s] *= a[i][key.get(i)];

			s++;
			++kit;
		}
	}
};

/*! \brief Calculate aint 2D
 *
 *
 */
template<typename vector, unsigned int np>
struct calculate_aint<2,vector,np>
{

	/*! \brief Calculate the coefficients of the interpolation a_int for one particles
	 *         having the 1D kernel values
	 *
	 * \param sz size of stencil for the interpolation
	 * \param a_int coefficients on the stencil points
	 * \param a coefficients of the kernel for each direction, consider
	 *          that for 3 dimensions the kernel is the multiplication
	 *          the 1D kernel on each direction. The "a" array store the
	 *          calculated coefficient of the 1D kernel on each direction
	 *
	 */
	static inline void value(size_t (& sz)[vector::dims],
			                 typename vector::stype a_int[openfpm::math::pow(np,vector::dims)],
							 typename vector::stype (& a)[vector::dims][np])
	{
		size_t s = 0;
		for (size_t i = 0 ; i < np ; i++)
		{
			for (size_t j = 0 ; j < np ; j++)
			{
				a_int[s] = a[0][j]*a[1][i];

				s++;
			}
		}
	}
};

/*! \brief Calculate aint
 *
 *
 */
template<typename vector, unsigned int np>
struct calculate_aint<3,vector,np>
{
	/*! \brief Calculate the coefficients of the interpolation a_int for one particles
	 *         having the 1D kernel values
	 *
	 * \param sz size of stencil for the interpolation
	 * \param a_int coefficients on the stencil points
	 * \param a coefficients of the kernel for each direction, consider
	 *          that for 3 dimensions the kernel is the multiplication
	 *          the 1D kernel on each direction. The "a" array store the
	 *          calculated coefficient of the 1D kernel on each direction
	 *
	 */
	static inline void value(size_t (& sz)[vector::dims],
			                 typename vector::stype a_int[openfpm::math::pow(np,vector::dims)],
							 typename vector::stype (& a)[vector::dims][np])
	{
		size_t s = 0;
		for (size_t i = 0 ; i < np ; i++)
		{
			for (size_t j = 0 ; j < np ; j++)
			{
				for (size_t k = 0 ; k < np ; k++)
				{
					a_int[s] = a[0][k]*a[1][j]*a[2][i];

					s++;
				}
			}
		}
	}
};

/*! \brief return the sub-domain where this particle must be interpolated
 *
 * \param p particle position
 *
 * \return the sub-domain id
 *
 */
template<typename vector, typename grid>
inline size_t getSub(Point<vector::dims,typename vector::stype> & p,
		             const CellList<vector::dims,typename vector::stype,Mem_fast<>,shift<vector::dims,typename vector::stype>> & geo_cell,
					 grid & gd)
{
	size_t cell = geo_cell.getCell(p);

	for (size_t i = 0 ; i < geo_cell.getNelements(cell) ; i++)
	{
		size_t ns = geo_cell.get(cell,i);

		if (gd.getDecomposition().getSubDomain(ns).isInside(p))
			return ns;
	}

    // If we end up here it mean that the particle decomposition and the grid decomposition are at machine precision level
    // different. To recover we shift the particle of a machine precision correction inside the domain.

    typename vector::stype dx[vector::dims];
    typename vector::stype best_dx[vector::dims];
    typename vector::stype best_tot_dx;
    long int best_sub;

	Box<vector::dims,typename vector::stype> sub_dom = gd.getDecomposition().getSubDomain(0);

	for (size_t j = 0 ; j < vector::dims ; j++)
	{
			if (sub_dom.getLow(j) > p.get(j))
			{dx[j] = 2*(sub_dom.getLow(j) - p.get(j));}
			else if (sub_dom.getHigh(j) < p.get(j))
			{dx[j] = 2*(sub_dom.getHigh(j) - p.get(j));}
			else {dx[j] = 0;}
	}

	typename vector::stype tot_dx = 0.0;

	for (size_t j = 0 ; j < vector::dims ; j++)
	{tot_dx += fabs(dx[j]);}

	best_tot_dx = tot_dx;
	for (size_t j = 0 ; j < vector::dims ; j++)
	{best_dx[j] = dx[j];}
	best_sub = 0;

    for (size_t i = 1 ; i < gd.getDecomposition().getNSubDomain() ; i++)
    {
		Box<vector::dims,typename vector::stype> sub_dom = gd.getDecomposition().getSubDomain(i);

		for (size_t j = 0 ; j < vector::dims ; j++)
		{
			if (sub_dom.getLow(j) > p.get(j))
			{dx[j] = 2*(sub_dom.getLow(j) - p.get(j));}
			else if (sub_dom.getHigh(j) < p.get(j))
			{dx[j] = 2*(sub_dom.getHigh(j) - p.get(j));}
			else {dx[j] = 0;}
		}

		typename vector::stype tot_dx = 0.0;

		for (size_t j = 0 ; j < vector::dims ; j++)
		{tot_dx += fabs(dx[j]);}

		if (tot_dx < best_tot_dx)
		{
			best_tot_dx = tot_dx;
			for (size_t j = 0 ; j < vector::dims ; j++)
			{best_dx[j] = dx[j];}
			best_sub = i;
		}
	}

	// shift xp

	for (size_t j = 0 ; j < vector::dims ; j++)
	{p.get(j) += best_dx[j];}

	return best_sub;
}

/*! \brief calculate the interpolation for one point
 *
 * \tparam vector of particles
 * \tparam kernel type
 *
 */
template<typename vector,typename kernel>
struct inte_calc_impl
{
	//! Type of the calculations
	typedef typename vector::stype arr_type;

	/*! \brief M2P or P2M for one particle
	 *
	 * \tparam prp_g property to interpolate from(M2P) or to(P2M) for grid
	 * \tparam prp_v property to interpolate to(M2P) or from(P2M) for vector
	 * \tparam m2p_or_p2m mesh to particle or mesh to particle interpolation
	 *
	 * \param it iterator used to retrieve the particle p for interpolation
	 * \param vd vector of particles
	 * \param domain simulation domain
	 * \param ip index of the grid on each direction (1D) used for interpolation
	 * \param gd interpolation grid
	 * \param dx spacing on each direction
	 * \param xp Point that store the position of xp
	 * \param a_int coefficients on the stencil points
	 * \param a coefficients of the kernel for each direction, consider
	 *          that for 3 dimensions the kernel is the multiplication
	 *          the 1D kernel on each direction. The "a" array store the
	 *          calculated coefficient of the 1D kernel on each direction
	 * \param x position of
	 * \param sz grid size
	 * \param geo_cell cell list to convert particle position into sub-domain id
	 * \param offsets array where are stored the linearized offset of the
	 *        kernel stencil for each local grid (sub-domain)
	 *
	 */
	template<unsigned int prp_g, unsigned int prp_v, unsigned int m2p_or_p2m, unsigned int np_a_int, typename iterator, typename grid>
																	 static inline void inte_calc(iterator & it,
																	 vector & vd,
																	 const Box<vector::dims,typename vector::stype> & domain,
																	 int (& ip)[vector::dims][kernel::np],
																	 grid & gd,
																	 const typename vector::stype (& dx)[vector::dims],
																	 typename vector::stype (& xp)[vector::dims],
																	 typename vector::stype (& a_int)[np_a_int],
																	 typename vector::stype (& a)[vector::dims][kernel::np],
																	 typename vector::stype (& x)[vector::dims][kernel::np],
																	 size_t (& sz)[vector::dims],
																	 const CellList<vector::dims,typename vector::stype,Mem_fast<>,shift<vector::dims,typename vector::stype>> & geo_cell,
																	 openfpm::vector<agg_arr<openfpm::math::pow(kernel::np,vector::dims)>> & offsets)
	{
		auto key_p = it.get();

		Point<vector::dims,typename vector::stype> p = vd.getPos(key_p);

		// On which sub-domain we interpolate the particle
		size_t sub = getSub<vector>(p,geo_cell,gd);

		typename vector::stype x0[vector::dims];

		// calculate the position of the particle in the grid
		// coordinated
		for (size_t i = 0 ; i < vector::dims ; i++)
			x0[i] = (p.get(i)-domain.getLow(i))*dx[i];

		// convert into integer
		for (size_t i = 0 ; i < vector::dims ; i++)
			ip[i][0] = (int)x0[i];

		// convert the global grid position into local grid position
		grid_key_dx<vector::dims> base;

		for (size_t i = 0 ; i < vector::dims ; i++)
			base.set_d(i,ip[i][0] - gd.getLocalGridsInfo().get(sub).origin.get(i) - (long int)kernel::np/2 + 1);

		// convenient grid of points

		for (size_t j = 0 ; j < kernel::np-1 ; j++)
		{
			for (size_t i = 0 ; i < vector::dims ; i++)
				ip[i][j+1] = (int)ip[i][j]+1;
		}

		for (size_t i = 0 ; i < vector::dims ; i++)
			xp[i] = x0[i] - ip[i][0];

		for (long int j = 0 ; j < kernel::np ; j++)
		{
			for (size_t i = 0 ; i < vector::dims ; i++)
				x[i][j] = - xp[i] + typename vector::stype((long int)j - (long int)kernel::np/2 + 1);
		}

		for (size_t j = 0 ; j < kernel::np ; j++)
		{
			for (size_t i = 0 ; i < vector::dims ; i++)
				a[i][j] = kernel::value(x[i][j],j);
		}

		calculate_aint<vector::dims,vector,kernel::np>::value(sz,a_int,a);

		grid_dist_lin_dx k_dist_lin;
		k_dist_lin.setSub(sub);

		size_t k = 0;
		auto & gs_info = gd.get_loc_grid(sub).getGrid();

		size_t lin_base = gs_info.LinId(base);

		for (size_t i = 0 ; i < openfpm::math::pow(kernel::np,vector::dims) ; i++)
		{
			size_t lin = offsets.get(sub).ele[k] + lin_base;
			k_dist_lin.getKeyRef() = lin;

			inte_template<kernel::np,prp_g,prp_v,m2p_or_p2m>::value(gd,vd,k_dist_lin,key_p,a_int,k);

			k++;
		}
	}
};

/*! \brief Main class for interpolation Particle to mest p2m and Mesh to particle m2p
 *
 * This function is the main class to interpolate from particle to mesh and mesh to particle
 *
 * \tparam vector type of vector for interpolation
 * \tparam grid type of grid for interpolation
 * \tparam interpolation kernel
 *
 */
template<typename vector,typename grid, typename kernel>
class interpolate
{
	//! Cell list used to convert particles position to sub-domain
	CellList<vector::dims,typename vector::stype,Mem_fast<>,shift<vector::dims,typename vector::stype>> geo_cell;

	/*! Structure to order boxes by volume
	 *
	 *
	 *
	 */
	struct Box_vol
	{
		//! Box
		Box<vector::dims,size_t> bv;

		//! Calculated volume
		size_t vol;

		/*! \brief operator to reorder
		 *
		 * \param bv box to compare with
		 *
		 * \return true if bv has volume bigger volume
		 *
		 */
		bool operator<(const Box_vol & bv)
		{
			return vol < bv.vol;
		}
	};

	//! particles
	vector & vd;

	//! grid
	grid & gd;

	//! Type of the calculations
	typedef typename vector::stype arr_type;

	/*! \brief It calculate the interpolation stencil offsets
	 *
	 * \param offsets array where to store the linearized offset of the
	 *        kernel stencil for each local-grid (sub-domain)
	 * \param sz kernel stencil points in each direction
	 *
	 */
	void calculate_the_offsets(openfpm::vector<agg_arr<openfpm::math::pow(kernel::np,vector::dims)>> & offsets, size_t (& sz)[vector::dims])
	{
		offsets.resize(gd.getN_loc_grid());

		for (size_t i = 0 ; i < offsets.size() ; i++)
		{
			auto & loc_grid = gd.get_loc_grid(i);
			const grid_sm<vector::dims,void> & gs = loc_grid.getGrid();

			grid_sm<vector::dims,void> gs2(sz);
			grid_key_dx_iterator<vector::dims> kit2(gs2);

			size_t k = 0;

			while (kit2.isNext())
			{
				auto key = kit2.get();

				long int lin = gs.LinId(key);

				offsets.get(i).ele[k] = lin;

				++k;
				++kit2;
			}
		}
	}

public:

	/*! \brief construct an interpolation object between a grid and a vector
	 *
	 * When possible and easy to do we suggest to retain the object
	 *
	 * \param vd interpolation vector
	 * \param gd interpolation grid
	 *
	 */
	interpolate(vector & vd, grid & gd)
	:vd(vd),gd(gd)
	{
		// get the processor bounding box in grid units
		Box<vector::dims,typename vector::stype> bb = gd.getDecomposition().getProcessorBounds();
		Box<vector::dims,typename vector::stype> bunit = gd.getDecomposition().getCellDecomposer().getCellBox();

		size_t div[vector::dims];

		for (size_t i = 0 ; i < vector::dims ; i++)
			div[i] = (bb.getHigh(i) - bb.getLow(i)) / bunit.getHigh(i);

		geo_cell.Initialize(bb,div);

		// Now draw the domain into the cell list

		auto & dec = gd.getDecomposition();

		for (size_t i = 0 ; i < dec.getNSubDomain() ; i++)
		{
			const Box<vector::dims,typename vector::stype> & bx = dec.getSubDomain(i);

			// get the cells this box span
			const grid_key_dx<vector::dims> p1 = geo_cell.getCellGrid(bx.getP1());
			const grid_key_dx<vector::dims> p2 = geo_cell.getCellGrid(bx.getP2());

			// Get the grid and the sub-iterator
			auto & gi = geo_cell.getGrid();
			grid_key_dx_iterator_sub<vector::dims> g_sub(gi,p1,p2);

			// add the box-id to the cell list
			while (g_sub.isNext())
			{
				auto key = g_sub.get();
				geo_cell.addCell(gi.LinId(key),i);
				++g_sub;
			}
		}
	};

	/*! \brief Interpolate particles to mesh
	 *
	 * Most of the time the particle set and the mesh are the same
	 * as the one used in the constructor. They also can be different
	 * as soon as they have the same decomposition
	 *
	 * \param gd grid or mesh
	 * \param vd particle set
	 *
	 */
	template<unsigned int prp_v, unsigned int prp_g> void p2m(vector & vd, grid & gd)
	{
#ifdef SE_CLASS1

		if (!vd.getDecomposition().is_equal_ng(gd.getDecomposition()) )
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error: the distribution of the vector of particles" <<
					" and the grid is different. In order to interpolate the two data structure must have the" <<
					" same decomposition" << std::endl;

			ACTION_ON_ERROR(INTERPOLATION_ERROR_OBJECT)
		}

#endif

		Box<vector::dims,typename vector::stype> domain = vd.getDecomposition().getDomain();

		// grid spacing
		typename vector::stype dx[vector::dims];

		for (size_t i = 0 ; i < vector::dims ; i++)
			dx[i] = 1.0/gd.spacing(i);

		size_t sz[vector::dims];

		for (size_t i = 0 ; i < vector::dims ; i++)
			sz[i] = kernel::np;

		// Precalculate the offset for each sub-sub-domain
		openfpm::vector<agg_arr<openfpm::math::pow(kernel::np,vector::dims)>> offsets;

		calculate_the_offsets(offsets,sz);

		// point position
		typename vector::stype xp[vector::dims];

		int ip[vector::dims][kernel::np];
		typename vector::stype x[vector::dims][kernel::np];
		typename vector::stype a[vector::dims][kernel::np];

		typename vector::stype a_int[openfpm::math::pow(kernel::np,vector::dims)];

		auto it = vd.getDomainIterator();

		while (it.isNext() == true)
		{
			inte_calc_impl<vector,kernel>::template inte_calc<prp_g,prp_v,inte_p2m,openfpm::math::pow(kernel::np,vector::dims)>(it,vd,domain,ip,gd,dx,xp,a_int,a,x,sz,geo_cell,offsets);

			++it;
		}
	}

	/*! \brief Interpolate mesh to particle
	 *
	 * Most of the time the particle set and the mesh are the same
	 * as the one used in the constructor. They also can be different
	 * as soon as they have the same decomposition
	 *
	 * \param gd grid or mesh
	 * \param vd particle set
	 *
	 */
	template<unsigned int prp_g, unsigned int prp_v> void m2p(grid & gd, vector & vd)
	{
#ifdef SE_CLASS1

		if (!vd.getDecomposition().is_equal_ng(gd.getDecomposition()) )
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error: the distribution of the vector of particles" <<
					" and the grid is different. In order to interpolate the two data structure must have the" <<
					" same decomposition" << std::endl;

			ACTION_ON_ERROR(INTERPOLATION_ERROR_OBJECT)
		}

#endif

		Box<vector::dims,typename vector::stype> domain = vd.getDecomposition().getDomain();

		// grid spacing
		typename vector::stype dx[vector::dims];

		for (size_t i = 0 ; i < vector::dims ; i++)
			dx[i] = 1.0/gd.spacing(i);

		// point position
		typename vector::stype xp[vector::dims];

		int ip[vector::dims][kernel::np];
		typename vector::stype x[vector::dims][kernel::np];
		typename vector::stype a[vector::dims][kernel::np];

		size_t sz[vector::dims];

		for (size_t i = 0 ; i < vector::dims ; i++)
			sz[i] = kernel::np;

		// Precalculate the offset for each sub-sub-domain
		openfpm::vector<agg_arr<openfpm::math::pow(kernel::np,vector::dims)>> offsets;

		calculate_the_offsets(offsets,sz);

		//		grid_cpu<vector::dims,aggregate<typename vector::stype>> a_int(sz);
		//		a_int.setMemory();
		typename vector::stype a_int[openfpm::math::pow(kernel::np,vector::dims)];

		auto it = vd.getDomainIterator();

		while (it.isNext() == true)
		{
			inte_calc_impl<vector,kernel>::template inte_calc<prp_g,prp_v,inte_m2p,openfpm::math::pow(kernel::np,vector::dims)>(it,vd,domain,ip,gd,dx,xp,a_int,a,x,sz,geo_cell,offsets);

			++it;
		}
	}

};

#endif /* OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_HPP_ */
