/*
 * Average.hpp
 *
 *  Created on: Nov 18, 2015
 *      Author: Pietro Incardona
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_AVERAGE_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_AVERAGE_HPP_


#define CENTRAL 0
#define CENTRAL_B_ONE_SIDE 1
#define FORWARD 2
#define BACKWARD 3

#include "util/mathutil.hpp"
#include "Vector/map_vector.hpp"
#include "Grid/comb.hpp"
#include "FiniteDifference/util/common.hpp"
#include "util/util_num.hpp"

/*! \brief Average
 *
 * \tparam d on which dimension average
 * \tparam Field which field average
 * \tparam impl which implementation
 *
 */
template<unsigned int d, typename Field, typename Sys_eqs, unsigned int impl=CENTRAL>
class Avg
{
	/*! \brief Create the row of the Matrix
	 *
	 * \tparam ord
	 *
	 */
	inline static void value(const grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}

	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case on non staggered case this function just return pos, in case of staggered,
	 *  it calculate where the operator is calculated on a staggered grid
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}
};

/*! \brief Average on direction i
 *
 *
 */
template<unsigned int d, typename arg, typename Sys_eqs>
class Avg<d,arg,Sys_eqs,CENTRAL>
{
	public:

	/*! \brief fill the row
	 *
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap , const grid_sm<Sys_eqs::dims,void> & gs, typename Sys_eqs::stype (& spacing )[Sys_eqs::dims] , std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		// if the system is staggered the CENTRAL derivative is equivalent to a forward derivative
		if (is_grid_staggered<Sys_eqs>::value())
		{
			Avg<d,arg,Sys_eqs,BACKWARD>::value(g_map,kmap,gs,spacing,cols,coeff);
			return;
		}

		long int old_val = kmap.getKeyRef().get(d);
		kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
		arg::value(g_map,kmap,gs,spacing,cols,coeff/2);
		kmap.getKeyRef().set_d(d,old_val);


		old_val = kmap.getKeyRef().get(d);
		kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
		arg::value(g_map,kmap,gs,spacing,cols,coeff/2);
		kmap.getKeyRef().set_d(d,old_val);
	}


	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case on non staggered case this function just return pos, in case of staggered,
	 *  it calculate where the operator is calculated on a staggered grid
	 *
	 *  \param pos from the position
	 *  \param fld Field we are deriving, if not provided the function just return pos
	 *  \param s_pos position of the properties in the staggered grid
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, long int fld = -1, const openfpm::vector<comb<Sys_eqs::dims>> & s_pos = openfpm::vector<comb<Sys_eqs::dims>>())
	{
		if (is_grid_staggered<Sys_eqs>::value())
		{
			if (fld == -1)
				return pos;

			if (s_pos[fld][d] == 1)
			{
				grid_key_dx<Sys_eqs::dims> ret = pos;
				ret.set_d(d,0);
				return pos;
			}
			else
			{
				grid_key_dx<Sys_eqs::dims> ret = pos;
				ret.set_d(d,1);
				return pos;
			}
		}

		return pos;
	}
};

/*! \brief Average FORWARD on direction i
 *
 *
 */
template<unsigned int d, typename arg, typename Sys_eqs>
class Avg<d,arg,Sys_eqs,FORWARD>
{
	public:

	/*! \brief fill the row
	 *
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap , const grid_sm<Sys_eqs::dims,void> & gs, typename Sys_eqs::stype (& spacing )[Sys_eqs::dims] , std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{

		long int old_val = kmap.getKeyRef().get(d);
		kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
		arg::value(g_map,kmap,gs,spacing,cols,coeff/2);
		kmap.getKeyRef().set_d(d,old_val);

		// backward
		arg::value(g_map,kmap,gs,spacing,cols,coeff/2);
	}


	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case on non staggered case this function just return pos, in case of staggered,
	 *  it calculate where the operator is calculated on a staggered grid
	 *
	 *  \param pos from the position
	 *  \param fld Field we are deriving, if not provided the function just return pos
	 *  \param s_pos position of the properties in the staggered grid
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, long int fld = -1, const openfpm::vector<comb<Sys_eqs::dims>> & s_pos = openfpm::vector<comb<Sys_eqs::dims>>())
	{
		return pos;
	}
};

/*! \brief Average BACKWARD on direction i
 *
 *
 */
template<unsigned int d, typename arg, typename Sys_eqs>
class Avg<d,arg,Sys_eqs,BACKWARD>
{
	public:

	/*! \brief fill the row
	 *
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap , const grid_sm<Sys_eqs::dims,void> & gs, typename Sys_eqs::stype (& spacing )[Sys_eqs::dims], std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		long int old_val = kmap.getKeyRef().get(d);
		kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
		arg::value(g_map,kmap,gs,spacing,cols,coeff/2);
		kmap.getKeyRef().set_d(d,old_val);

		// forward
		arg::value(g_map,kmap,gs,spacing,cols,coeff/2);
	}


	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case on non staggered case this function just return pos, in case of staggered,
	 *  it calculate where the operator is calculated on a staggered grid
	 *
	 *  \param pos from the position
	 *  \param fld Field we are deriving, if not provided the function just return pos
	 *  \param s_pos position of the properties in the staggered grid
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, long int fld = -1, const openfpm::vector<comb<Sys_eqs::dims>> & s_pos = openfpm::vector<comb<Sys_eqs::dims>>())
	{
		return pos;
	}
};

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_AVERAGE_HPP_ */
