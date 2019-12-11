/*
 * sum.hpp
 *
 *  Created on: Oct 13, 2015
 *      Author: i-bird
 *  Modified on: Dec 09, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_

#include <boost/mpl/vector.hpp>
#include "config.h"
#include <unordered_map>
#include "util/for_each_ref.hpp"

/*! \brief It model an expression expr1 + ... exprn
 *
 * \tparam expr1_type Type of first expression.
 * \tparam expr2_type Type of second expression.
 * \tparam expr1_type::sys_nn System of equations.
 *
 * ## Example
 *
 * \snippet FDScheme_unit_tests.hpp sum example
 *
 */
template<typename expr1_type, typename expr2_type>
class sum
{

public:

  typedef typename expr1_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;
  
  expr1_type expr1; /**< First expression to sum. */
  expr2_type expr2; /**< Second expression to sum. */


  // Constructors
  /**
   * @fn sum
   * @brief Default constructor.
   */
  sum(): expr1{expr1_type{}}, expr2{expr2_type{}} {}
  
  /**
   * @fn sum(expr1_type&, expr2_type&)
   * @brief Constructor.
   * @param[in] _expr1 First expression to add.
   * @param[in] _expr2 Second expression to add.
   */
  sum(const expr1_type & _expr1, const expr2_type & _expr2) : expr1{_expr1}, expr2{_expr2} {}
  
  
  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map Grid mapping, it convert grid position into vector/Matrix row
   * \param kmap grid position
   * \param gs grid information
   * \param spacing grid spacing
   * \param cols unordered map contain the map colum -> value
   * \param coeff it contain an additional actual coefficients in front of the values
   * \param imp_pos Position in the cell where to compute the value. (Important for staggered grids).
   *
   * ### Example
   *
   * \snippet FDScheme_unit_tests.hpp sum example
   *
   */
  inline void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition::extended_type>::type & g_map,
		    grid_dist_key_dx<Sys_eqs::dims> & kmap,
		    const grid_sm<Sys_eqs::dims,void> & gs,
		    typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
		    std::unordered_map<long int,typename Sys_eqs::stype > & cols,
		    typename Sys_eqs::stype coeff,
		    comb<Sys_eqs::dims> imp_pos) const
  {
    expr1.value(g_map,kmap,gs,spacing,cols,coeff,imp_pos);
    expr2.value(g_map,kmap,gs,spacing,cols,coeff,imp_pos);
  }


};

/*! \brief It ancapsulate the minus operation
 *
 * \tparam expr1_type Type of expression.
 * \tparam Sys_eqs System of equations.
 *
 */
template<typename expr1_type>
class minus
{

public:

  typedef typename expr1_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;
  
  expr1_type expr; /**< Expression to substract. */

  // Constructors
  /**
   * @fn minus
   * @brief Default constructor.
   */
  minus(): expr{expr1_type{}} {}
  
  /**
   * @fn minus(expr1_type&)
   * @brief Constructor.
   * @param[in] env Environment object.
   * @param[in] total_mag Total magnetization of the system.
   */
  minus(const expr1_type & _expr) : expr{_expr} {}

  /*! \brief Create the row of the Matrix
   *
   * \tparam ord
   *
   * \snippet FDScheme_unit_tests.hpp minus example
   *
   * \param g_map Grid mapping, it convert grid position into vector/Matrix row
   * \param kmap grid position
   * \param gs grid information
   * \param spacing grid spacing
   * \param cols unordered map contain the map colum -> value
   * \param coeff it contain an additional actual coefficients in front of the values
   *
   */
  inline void value(const typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition::extended_type>::type & g_map,
		    grid_dist_key_dx<Sys_eqs::dims> & kmap,
		    const grid_sm<Sys_eqs::dims,void> & gs,
		    typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
		    std::unordered_map<long int,typename Sys_eqs::stype > & cols,
		    typename Sys_eqs::stype coeff,
		    comb<Sys_eqs::dims> imp_pos) const
  {
    expr.value(g_map,kmap,gs,spacing,cols,-coeff,imp_pos);
  }


};

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_SUM_HPP_ */
