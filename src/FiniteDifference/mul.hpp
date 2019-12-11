/*
 * mul.hpp
 *
 *  Created on: Oct 22, 2015
 *      Author: i-bird
 *  Modified on: Dec 06, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_MUL_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_MUL_HPP_

#include <typeinfo>
#include "util/util_debug.hpp"
#include "util/util_num.hpp"

#define HAS_VAL 1
#define HAS_POS_VAL 2

//! Evaluate the constant field function
template <unsigned int, typename T>
struct has_val
{
  //! evaluate
  static float call_val()
  {
    std::cerr << "Error the type " << demangle(typeid(T).name()) << "interpreted as constant field, does not have a function val() or val_pos(), please see the numeric examples in Finite Differences for more information\n";
    return 0;
  }
};

//! Evaluate the constant field function
template <typename T>
struct has_val<HAS_VAL,T>
{
  //! evaluate
  static decltype(T::val()) call_val()
  {
    return T::val();
  }
};

/*! \brief It model an expression expr1 * expr2
 *
 * \warning expr1 MUST be a constant expression only expr2 depend form variable, this requirement ensure linearity in the solving variable of the equations
 *
 * \tparam expr1_type Type of first expression.
 * \tparam expr2_type Type of second expression.
 * \tparam expr1_type::sys_nn System of equations.
 *
 */
template<typename expr1_type, typename expr2_type>
class mul
{
public:

  typedef typename expr1_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;
  
  expr1_type expr1; /**< First expression to sum. */
  expr2_type expr2; /**< Second expression to sum. */


  // Constructors
  /**
   * @fn mul
   * @brief Default constructor.
   */
  mul(): expr1{expr1_type{}}, expr2{expr2_type{}} {}
  
  /**
   * @fn mul(expr1_type&, expr2_type&)
   * @brief Constructor.
   * @param[in] _expr1 First expression to add.
   * @param[in] _expr2 Second expression to add.
   */
  mul(const expr1_type & _expr1, const expr2_type & _expr2) : expr1{_expr1}, expr2{_expr2} {}
  
  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map mapping grid
   * \param kmap position where the multiplication is calculated
   * \param gs Grid info
   * \param spacing of the grid
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
   * \param imp_pos Position in the cell where to compute the value. (Important for staggered grids).
   *
   */
  inline void value(const grid_dist_id<Sys_eqs::dims,typename Sys_eqs::stype,aggregate<size_t>,typename Sys_eqs::b_grid::decomposition::extended_type> & g_map,
		    grid_dist_key_dx<Sys_eqs::dims> & kmap,
		    const grid_sm<Sys_eqs::dims,void> & gs,
		    typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
		    std::unordered_map<long int,typename Sys_eqs::stype > & cols,
		    typename Sys_eqs::stype coeff,
		    comb<Sys_eqs::dims> imp_pos) const
  {
    coeff *= expr1.get(kmap,imp_pos);
    expr2.value(g_map,kmap,gs,spacing,cols,coeff,imp_pos);
  }
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_MUL_HPP_ */
