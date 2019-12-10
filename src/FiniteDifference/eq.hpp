/*
 * eq.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: i-bird
 *  Modified on: Dec 05, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_HPP_

#define EQS_FIELD 0
#define EQS_POS 1

//#define PERIODIC true
//#define NON_PERIODIC false

#include "util/util_num.hpp"
#include "FiniteDifference/util/common.hpp"
#include "Matrix/SparseMatrix.hpp"

// ---------------------------------------------------------------------------------------------------------------------------------------------

/*! \brief Equation
 *
 * It model an equation like expr1 = expr2
 *
 * \tparam expr1
 * \tparam expr2
 *
 */
template<typename expr1,typename expr2,typename Sys_eqs>
class Eq
{
  /*! \brief Create the row of the Matrix
   *
   * \tparam ord
   *
   */
  template<unsigned int ord=EQS_FIELD> static void value(const grid_key_dx<Sys_eqs::dims> & pos)
  {
    if (EQS_FIELD)
      value_f(pos);
    else
      value_s(pos);
  }

  /*! \brief fill the row
   *
   *
   */
  static openfpm::vector<cval<typename Sys_eqs::stype>> value_s(grid_key_dx<Sys_eqs::dims> & it)
  {
    return expr1::value_s(it) - expr2::value_s(it);
  }

  /*! \brief fill the row
   *
   *
   */
  static void value_f(grid_key_dx<Sys_eqs::dims> & it)
  {
    return expr1::value_s(it) - expr2::value_s(it);
  }
};

// ---------------------------------------------------------------------------------------------------------------------------------------------

// spatial position + value

template<unsigned int dim,typename T>
struct pos_val
{
  /*! \brief Initialize to zero the value
   *
   */
  pos_val()
  {
    value = 0.0;
  }

  grid_key_dx<dim> pos;
  T value;
};

// ---------------------------------------------------------------------------------------------------------------------------------------------

template<unsigned int f, typename Sys_eqs>
class Field
{
  typedef typename stub_or_real<Sys_eqs,Sys_eqs::dims,typename Sys_eqs::stype,typename Sys_eqs::b_grid::decomposition::extended_type::extended_type>::type map_grid;

public:

  typedef Sys_eqs sys_eqs_type;

  Field() {};

  /*! \brief fill the row
   *
   *
   */
  static void value(const map_grid & g_map,
		    grid_dist_key_dx<Sys_eqs::dims> & kmap,
		    const grid_sm<Sys_eqs::dims,void> & gs,
		    typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
		    std::unordered_map<long int,typename Sys_eqs::stype > & cols,
		    typename Sys_eqs::stype coeff)
  {
    cols[g_map.template get<0>(kmap)*Sys_eqs::nvar + f] += coeff;
  }

  /*! \brief
   *
   *
   */
  static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    return grid_key_dx<Sys_eqs::dims>(s_pos[f]);
  }
};

// ---------------------------------------------------------------------------------------------------------------------------------------------

inline size_t mat_factor(size_t nvar, size_t sz, const size_t ord)
{
  return nvar;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------


/**
 * \class coeff
 * \brief Creates a coefficient object to be used in the equations.
 * \tparam coeff_type Type of coefficient. It could be just a double or a grid_dist_id, for example.
 * \tparam Sys_eqs Struct with information regarding the system of equations.
 */
template<typename coeff_type, typename Sys_eqs>
class coeff {

public:

  const typename std::conditional<has_get<coeff_type,Sys_eqs::dims>::value,coeff_type &,coeff_type>::type c;

  typedef Sys_eqs sys_eqs_type; /**< Extra helper type. */

  /**
   * \fn coeff(const coeff_type &)
   * \brief Constructor.
   */
  coeff(const coeff_type & c_) : c{c_} {}

  /**
   * \fn get(grid_dist_key_dx<Sys_eqs::dims> &)
   * \brief Compute the value of the coeff. This function is called when the coefficient changes from point to point in the grid.
   * \return Value of the coefficient on a specific point in the grid (key).
   */
  template<typename U = coeff_type, typename std::enable_if<has_get<U,Sys_eqs::dims>::value,int>::type = 0>
  typename Sys_eqs::stype get(grid_dist_key_dx<Sys_eqs::dims> & key) const {
    return c.template get(key);
  }

  /**
   * \fn get(grid_dist_key_dx<Sys_eqs::dims> &)
   * \brief Compute the value of the coeff. This function is called when the coefficient is just a number that does not change from point to point.
   * \return Value of the coefficient.
   */
  template<typename U = coeff_type, typename std::enable_if<!has_get<U,Sys_eqs::dims>::value,int>::type = 0>
  typename Sys_eqs::stype get(grid_dist_key_dx<Sys_eqs::dims> & key) const {
    return c;
  }
};

// ---------------------------------------------------------------------------------------------------------------------------------------------

#include "mul.hpp"
#include "Average.hpp"
#include "Derivative.hpp"
#include "sum.hpp"
#include "Laplacian.hpp"

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_EQ_HPP_ */
