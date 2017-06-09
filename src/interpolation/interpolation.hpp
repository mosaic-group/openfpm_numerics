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

#define INTERPOLATION_ERROR_OBJECT std::runtime_error("Runtime interpolation error");

constexpr int inte_m2p = 0;
constexpr int inte_p2m = 1;

template<typename T>
struct mul_inte
{
	inline static void value(T & result, const T & coeff, const T & src)
	{
		result += coeff * src;
	}
};

template<typename T, unsigned int N1>
struct mul_inte<T[N1]>
{
	inline static void value(T (& result)[N1], const T & coeff, const T (& src)[N1])
	{
		for (size_t i = 0 ; i < N1 ; i++)
			result[i] += coeff * src[i];
	}
};

template<typename T, unsigned int N1, unsigned int N2>
struct mul_inte<T[N1][N2]>
{
	inline static void value(T (& result)[N1][N2], const T & coeff, const T (& src)[N1][N2])
	{
		for (size_t i = 0 ; i < N1 ; i++)
			for (size_t j = 0 ; j < N2 ; j++)
				result[i][j] += coeff * src[i][j];
	}
};

template<unsigned int prp_g, unsigned int prp_v,unsigned int m2p_or_p2m>
struct inte_template
{
	template<typename grid, typename vector, typename iterator, typename key_type> inline static void value(grid & gd,
			                                                          vector & vd,
																	  const grid_dist_key_dx<vector::dims> & k_dist,
																	  iterator & key_p,
																	  const grid_cpu<vector::dims,aggregate<typename vector::stype>> & a_int,
																	  const key_type & key)
	{
		mul_inte<typename std::remove_reference<decltype(gd.template get<prp_g>(k_dist))>::type>::value(gd.template get<prp_g>(k_dist),a_int.template get<0>(key),vd.template getProp<prp_v>(key_p));
	}
};

template<unsigned int prp_g, unsigned int prp_v>
struct inte_template<prp_g,prp_v,inte_m2p>
{
	template<typename grid, typename vector, typename iterator, typename key_type> inline static void value(grid & gd,
			                                                          vector & vd,
																	  const grid_dist_key_dx<vector::dims> & k_dist,
																	  iterator & key_p,
																	  const grid_cpu<vector::dims,aggregate<typename vector::stype>> & a_int,
																	  const key_type & key)
	{
		mul_inte<typename std::remove_reference<decltype(gd.template get<prp_g>(k_dist))>::type>::value(vd.template getProp<prp_v>(key_p),a_int.template get<0>(key),gd.template get<prp_g>(k_dist));
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
		             const CellList<vector::dims,typename vector::stype,Mem_fast<vector::dims,typename vector::stype>,shift<vector::dims,typename vector::stype>> & geo_cell,
					 grid & gd)
{
	size_t cell = geo_cell.getCell(p);

	for (size_t i = 0 ; i < geo_cell.getNelements(cell) ; i++)
	{
		size_t ns = geo_cell.get(cell,i);

		if (gd.getDecomposition().getSubDomain(ns).isInside(p))
			return ns;
	}

#ifdef SE_CLASS1
	std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << std::endl;
	ACTION_ON_ERROR(INTERPOLATION_ERROR_OBJECT)
#endif

	return (size_t)-1;
}


template<unsigned int prp_g, unsigned int prp_v, unsigned int m2p_or_p2m, typename kernel, typename iterator, typename vector, typename grid, typename grid_inte>
																 inline void inte_calc(iterator & it,
		                     	 	 	 	 	 	 	 	 	 const vector & vd,
																 const Box<vector::dims,typename vector::stype> & domain,
																 int (& ip)[vector::dims][kernel::np],
																 grid & gd,
																 const typename vector::stype (& dx)[vector::dims],
																 typename vector::stype (& xp)[vector::dims],
																 const grid_inte & a_int,
																 typename vector::stype (& a)[vector::dims][kernel::np],
																 typename vector::stype (& x)[vector::dims][kernel::np],
																 size_t (& sz)[vector::dims],
																 const CellList<vector::dims,typename vector::stype,Mem_fast<vector::dims,typename vector::stype>,shift<vector::dims,typename vector::stype>> & geo_cell)
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

	grid_sm<vector::dims,void> gs(sz);
	grid_key_dx_iterator<vector::dims> kit(gs);

	while (kit.isNext())
	{
		auto key = kit.get();

		a_int.template get<0>(key) = 1;

		for (size_t i = 0 ; i < vector::dims ; i++)
			a_int.template get<0>(key) *= a[i][key.get(i)];

		++kit;
	}

	grid_key_dx_iterator<vector::dims> kit2(gs);
	grid_dist_key_dx<vector::dims> k_dist;
	k_dist.setSub(sub);

	while (kit2.isNext())
	{
		auto key = kit2.get();

		for (size_t i = 0 ; i < vector::dims ; i++)
			k_dist.getKeyRef().set_d(i,key.get(i) + base.get(i));

		inte_template<prp_g,prp_v,m2p_or_p2m>::value(gd,vd,k_dist,key_p,a_int,key);

		++kit2;
	}
}


template<typename vector,typename grid, typename kernel>
class interpolate
{
	CellList<vector::dims,typename vector::stype,Mem_fast<vector::dims,typename vector::stype>,shift<vector::dims,typename vector::stype>> geo_cell;

	struct Box_vol
	{
		Box<vector::dims,size_t> bv;

		size_t vol;

		void operator<(const Box_vol & bv)
		{
			return vol < bv.vol;
		}
	};

	//! particles
	vector & vd;

	//! grid
	grid & gd;


public:

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

		// point position
		typename vector::stype xp[vector::dims];

		int ip[vector::dims][kernel::np];
		typename vector::stype x[vector::dims][kernel::np];
		typename vector::stype a[vector::dims][kernel::np];

		size_t sz[vector::dims];

		for (size_t i = 0 ; i < vector::dims ; i++)
			sz[i] = kernel::np;

		grid_cpu<vector::dims,aggregate<typename vector::stype>> a_int(sz);
		a_int.setMemory();

		auto it = vd.getDomainIterator();

		while (it.isNext() == true)
		{
			inte_calc<prp_g,prp_v,inte_p2m,kernel>(it,vd,domain,ip,gd,dx,xp,a_int,a,x,sz,geo_cell);

			++it;
		}
	}

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

		grid_cpu<vector::dims,aggregate<typename vector::stype>> a_int(sz);
		a_int.setMemory();

		auto it = vd.getDomainIterator();

		while (it.isNext() == true)
		{
			inte_calc<prp_g,prp_v,inte_m2p,kernel>(it,vd,domain,ip,gd,dx,xp,a_int,a,x,sz,geo_cell);

			++it;
		}
	}

};

#endif /* OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_HPP_ */
