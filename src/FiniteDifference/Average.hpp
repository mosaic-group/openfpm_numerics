/*
 * Average.hpp
 *
 *  Created on: Nov 18, 2015
 *      Author: Pietro Incardona
 *  Modified on: Dec 09, 2019
 *      Author: amfoggia
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
template<unsigned int d, typename expr_type, unsigned int impl=CENTRAL>
class Avg
{

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;
  
  /*! \brief Create the row of the Matrix
   *
   * \tparam ord
   *
   */
  inline void value(const grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff,comb<Sys_eqs::dims> imp_pos) const
  {
    std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
  }

  /*! \brief Calculate the position where the derivative is calculated
   *
   * In case of non staggered case this function just return a null grid_key, in case of staggered,
   *  it calculate how the operator shift the calculation in the cell
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const openfpm::vector<comb<Sys_eqs::dims>> & s_pos = openfpm::vector<comb<Sys_eqs::dims>>())
  {
    std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";
  }
};

/*! \brief Central average scheme on direction i
 *
 * \verbatim
 *
 *  +0.5     +0.5
 *   *---+---*
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class Avg<d,expr_type,CENTRAL>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;

  Avg() : expr{expr_type{}} {}

  Avg(expr_type expr_) : expr{expr_} {}

  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * stub_or_real it is just for change the argument type when testing, in normal
   * conditions it is a distributed map
   *
   * \param g_map It is the map explained in FDScheme
   * \param kmap position where the average is calculated
   * \param gs Grid info
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
   * \param imp_pos Position in the cell where to compute the value. (Important for staggered grids).
   *
   * ### Example
   *
   * \snippet FDScheme_unit_tests.hpp Usage of stencil derivative
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
    // if the system is staggered the CENTRAL derivative is equivalent to a forward derivative
    if (is_grid_staggered<Sys_eqs>::value())
      {
	Avg<d,expr_type,BACKWARD>{expr}.value(g_map,kmap,gs,spacing,cols,coeff,imp_pos);
	return;
      }

    long int old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
    expr.value(g_map,kmap,gs,spacing,cols,coeff/2,imp_pos);
    kmap.getKeyRef().set_d(d,old_val);


    old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
    expr.value(g_map,kmap,gs,spacing,cols,coeff/2,imp_pos);
    kmap.getKeyRef().set_d(d,old_val);
  }


  /*! \brief Calculate the position where the average is calculated
   *
   * It follow the same concept of central derivative
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    auto arg_pos = expr_type::position(pos,gs,s_pos);
    if (is_grid_staggered<Sys_eqs>::value())
      {
	if (arg_pos.get(d) == -1)
	  {
	    arg_pos.set_d(d,0);
	    return arg_pos;
	  }
	else
	  {
	    arg_pos.set_d(d,-1);
	    return arg_pos;
	  }
      }

    return arg_pos;
  }
};

/*! \brief FORWARD average on direction i
 *
 * \verbatim
 *
 *  +0.5    0.5
 *    +------*
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class Avg<d,expr_type,FORWARD>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;

  Avg() : expr{expr_type{}}{}

  Avg(expr_type expr_) : expr{expr_} {}

  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * stub_or_real it is just for change the argument type when testing, in normal
   * conditions it is a distributed map
   *
   * \param g_map It is the map explained in FDScheme
   * \param kmap position where the average is calculated
   * \param gs Grid info
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
   * \param imp_pos Position in the cell where to compute the value. (Important for staggered grids).
   *
   * ### Example
   *
   * \snippet FDScheme_unit_tests.hpp Usage of stencil derivative
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
    
    long int old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
    expr.value(g_map,kmap,gs,spacing,cols,coeff/2,imp_pos);
    kmap.getKeyRef().set_d(d,old_val);

    // backward
    expr.value(g_map,kmap,gs,spacing,cols,coeff/2,imp_pos);
  }


  /*! \brief Calculate the position in the cell where the average is calculated
   *
   * In case of non staggered case this function just return a null grid_key, in case of staggered,
   * the FORWARD scheme return the position of the staggered property
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    return expr_type::position(pos,gs,s_pos);
  }
};

/*! \brief First order BACKWARD derivative on direction i
 *
 * \verbatim
 *
 *  +0.5    0.5
 *    *------+
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class Avg<d,expr_type,BACKWARD>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;

  Avg() : expr{expr_type{}} {}

  Avg(expr_type expr_) : expr{expr_} {}

  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * stub_or_real it is just for change the argument type when testing, in normal
   * conditions it is a distributed map
   *
   * \param g_map It is the map explained in FDScheme
   * \param kmap position where the average is calculated
   * \param gs Grid info
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
   * \param imp_pos Position in the cell where to compute the value. (Important for staggered grids).
   *
   * ### Example
   *
   * \snippet FDScheme_unit_tests.hpp Usage of stencil derivative
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
    long int old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
    expr.value(g_map,kmap,gs,spacing,cols,coeff/2,imp_pos);
    kmap.getKeyRef().set_d(d,old_val);

    // forward
    expr.value(g_map,kmap,gs,spacing,cols,coeff/2,imp_pos);
  }


  /*! \brief Calculate the position in the cell where the average is calculated
   *
   * In case of non staggered case this function just return a null grid_key, in case of staggered,
   * the BACKWARD scheme return the position of the staggered property
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    return expr_type::position(pos,gs,s_pos);
  }
};

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_AVERAGE_HPP_ */
