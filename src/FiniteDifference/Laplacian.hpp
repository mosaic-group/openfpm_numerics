/*
 * Laplacian.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_LAPLACIAN_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_LAPLACIAN_HPP_


/*! \brief Laplacian second order on h (spacing)
 *
 * \tparam Field which field derive
 * \tparam impl which implementation
 *
 */
template<typename arg, typename Sys_eqs, unsigned int impl=CENTRAL>
class Lap
{
	/*! \brief Create the row of the Matrix
	 *
	 * \tparam ord
	 *
	 */
	inline static std::unordered_map<long int,typename Sys_eqs::stype> value(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}

	/*! \brief Calculate the position where the laplacian is calculated
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

/*! \brief Derivative on direction i
 *
 *
 */
template<typename arg, typename Sys_eqs>
class Lap<arg,Sys_eqs,CENTRAL>
{
public:

	/*! \brief fill the row
	 *
	 *
	 */
	inline static void value(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		// for each dimension
		for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
		{
			// forward
			if (Sys_eqs::boundary[i] == PERIODIC )
			{
				long int old_val = pos.get(i);
				pos.set_d(i, openfpm::math::positive_modulo(pos.get(i) + 1, gs.size(i)));
				arg::value(pos,gs,cols,coeff);
				pos.set_d(i,old_val);
			}
			else
			{
				long int old_val = pos.get(i);
				pos.set_d(i, pos.get(i) + 1);
				arg::value(pos,gs,cols,coeff);
				pos.set_d(i,old_val);
			}

			// backward
			if (Sys_eqs::boundary[i] == PERIODIC )
			{
				long int old_val = pos.get(i);
				pos.set_d(i, openfpm::math::positive_modulo(pos.get(i) - 1, gs.size(i)));
				arg::value(pos,gs,cols,coeff);
				pos.set_d(i,old_val);
			}
			else
			{
				long int old_val = pos.get(i);
				pos.set_d(i, pos.get(i) - 1);
				arg::value(pos,gs,cols,coeff);
				pos.set_d(i,old_val);
			}
		}

		arg::value(pos,gs,cols, - 2.0 * Sys_eqs::dims * coeff);
	}


	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case of non staggered case this function just return pos, in case of staggered,
	 * it calculate where the operator is calculated on a staggered grid
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const openfpm::vector<comb<Sys_eqs::dims>> & s_pos, long int & fld)
	{
		return pos;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_LAPLACIAN_HPP_ */
