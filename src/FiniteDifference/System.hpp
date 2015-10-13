/*
 * System.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SYSTEM_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SYSTEM_HPP_


/*! \brief System of equations
 *
 * This class model a system of equations
 *
 * \tparam dim Dimensionality of the system
 * \tparam nvf number of variable fields
 * \tparam ncf number of constants fields
 * \tparam eqs sets of equations
 *
 */
template<unsigned int dim, unsigned int nvf, unsigned int ncf, typename ... eqs>
class System
{
	// Define the number of constant fields
	typedef num_cfields boost::mpl::int_<nf>;

	// Define the number of variable fields
	typedef num_vfields boost::mpl::int_<nf>;

	// set of equations as boost::mpl::vector
	typedef eqs_v make_vactor<eqs>;

	/*! \brief Create the row of the Matrix
	 *
	 * \tparam ord
	 *
	 */
	template<unsigned int ord=EQS_FIELDS> void value(const grid_key_dx_dist<dim> & pos)
	{
		if (EQS_FIELDS)
			value_f(pos);
		else
			value_s(pos);
	}

	/*! \brief fill the row
	 *
	 *
	 */
	template<unsigned int eq_id> void value_s(grid_key_dx_dist<dim> & it)
	{
		boost::mpl::at<eqs_v,boost::mpl::int_<eq_id>>::type eq;

		eq.value_s(it);
	}

	/*! \brief fill the row
	 *
	 *
	 */
	template<unsigned int eq_id> void value_f(grid_key_dx_dist<dim> & it)
	{
		boost::mpl::at<eqs_v,boost::mpl::int_<eq_id>>::type eq;

		eq.value_f(it);
	}
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SYSTEM_HPP_ */
