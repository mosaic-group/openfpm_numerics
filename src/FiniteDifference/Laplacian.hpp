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
	/*! \brief Calculate which colums of the Matrix are non zero
	 *
	 * stub_or_real it is just for change the argument type when testing, in normal
	 * conditions it is a distributed map
	 *
	 * \param pos position where the laplacian is calculated
	 * \param gs Grid info
	 * \param cols non-zero colums calculated by the function
	 * \param coeff coefficent (constant in front of the derivative)
	 *
	 * ### Example
	 *
	 * \snippet FDScheme_unit_tests.hpp Laplacian usage
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap , const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}

	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * \param position where we are calculating the derivative
	 * \param gs Grid info
	 * \param s_pos staggered position of the properties
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs)
	{
		std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
	}
};

/*! \brief Laplacian second order approximation CENTRAL Scheme
 *
 * \verbatim

          1.0
           *
           | -4.0
   1.0 *---+---* 1.0
           |
           *
          1.0

 * \endverbatim
 *
 *
 */
template<typename arg, typename Sys_eqs>
class Lap<arg,Sys_eqs,CENTRAL>
{
public:


	/*! \brief Calculate which colums of the Matrix are non zero
	 *
	 * stub_or_real it is just for change the argument type when testing, in normal
	 * conditions it is a distributed map
	 *
	 * \param pos position where the laplacian is calculated
	 * \param gs Grid info
	 * \param cols non-zero colums calculated by the function
	 * \param coeff coefficent (constant in front of the derivative)
	 *
	 * ### Example
	 *
	 * \snippet FDScheme_unit_tests.hpp Laplacian usage
	 *
	 */
	inline static void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition>::type & g_map, grid_dist_key_dx<Sys_eqs::dims> & kmap , const grid_sm<Sys_eqs::dims,void> & gs, typename Sys_eqs::stype (& spacing )[Sys_eqs::dims] , std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff)
	{
		// for each dimension
		for (size_t i = 0 ; i < Sys_eqs::dims ; i++)
		{
			long int old_val = kmap.getKeyRef().get(i);
			kmap.getKeyRef().set_d(i, kmap.getKeyRef().get(i) + 1);
			arg::value(g_map,kmap,gs,spacing,cols,coeff/spacing[i]/spacing[i]);
			kmap.getKeyRef().set_d(i,old_val);

			old_val = kmap.getKeyRef().get(i);
			kmap.getKeyRef().set_d(i, kmap.getKeyRef().get(i) - 1);
			arg::value(g_map,kmap,gs,spacing,cols,coeff/spacing[i]/spacing[i]);
			kmap.getKeyRef().set_d(i,old_val);

			arg::value(g_map,kmap,gs,spacing,cols, - 2.0 * coeff/spacing[i]/spacing[i]);
		}
	}


	/*! \brief Calculate the position where the derivative is calculated
	 *
	 * In case of non staggered case this function just return a null grid_key, in case of staggered,
	 * the CENTRAL Laplacian scheme return the position of the staggered property
	 *
	 * \param position where we are calculating the derivative
	 * \param gs Grid info
	 * \param s_pos staggered position of the properties
	 *
	 */
	inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
	{
		return arg::position(pos,gs,s_pos);
	}
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_LAPLACIAN_HPP_ */
