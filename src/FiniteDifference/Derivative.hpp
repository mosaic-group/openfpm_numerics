/*
 * Derivative.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: Pietro Incardona
 *  Modified on: Dec 09, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_DERIVATIVE_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_DERIVATIVE_HPP_


#include "util/mathutil.hpp"
#include "Vector/map_vector.hpp"
#include "Grid/comb.hpp"
#include "FiniteDifference/util/common.hpp"
#include "util/util_num.hpp"
#include <unordered_map>
#include "FD_util_include.hpp"

/*! \brief Derivative second order on h (spacing)
 *
 * \tparam d on which dimension derive
 * \tparam Field which field derive
 * \tparam impl which implementation
 *
 */
template<unsigned int d, typename expr_type, unsigned int impl=CENTRAL>
class D
{

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;
  
  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param pos position where the derivative is calculated
   * \param gs Grid info
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
   *
   * ### Example
   *
   * \snippet FDScheme_unit_tests.hpp Usage of stencil derivative
   *
   */
  inline void value(const grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, std::unordered_map<long int,typename Sys_eqs::stype > & cols, typename Sys_eqs::stype coeff) const
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
   * \return where (in which cell) the derivative is calculated
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos,
						    const grid_sm<Sys_eqs::dims,void> & gs,
						    const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " only CENTRAL, FORWARD, BACKWARD derivative are defined";

    return pos;
  }
};

/*! \brief Second order central Derivative scheme on direction i
 *
 * \verbatim
 *
 *  -1        +1
 *   *---+---*
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class D<d,expr_type,CENTRAL>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;

  D() {}

  D(expr_type expr_) : expr{expr_} {}

  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map mapping grid
   * \param kmap position where the derivative is calculated
   * \param gs Grid info
   * \param spacing grid spacing
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
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
		    typename Sys_eqs::stype coeff) const
  {
    // if the system is staggered the CENTRAL derivative is equivalent to a forward derivative
    if (is_grid_staggered<Sys_eqs>::value())
      {
	D<d,expr_type,BACKWARD>{}.value(g_map,kmap,gs,spacing,cols,coeff);
	return;
      }

    long int old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
    expr.value(g_map,kmap,gs,spacing,cols,coeff/spacing[d]/2.0 );
    kmap.getKeyRef().set_d(d,old_val);


    old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
    expr.value(g_map,kmap,gs,spacing,cols,-coeff/spacing[d]/2.0 );
    kmap.getKeyRef().set_d(d,old_val);
  }


  /*! \brief Calculate the position where the derivative is calculated
   *
   * In case on non staggered case this function just return a null grid_key, in case of staggered,
   *  it calculate how the operator shift in the cell
   *
   \verbatim

   +--$--+
   |     |
   #  *  #
   |     |
   0--$--+

   # = velocity(y)
   $ = velocity(x)
   * = pressure

   \endverbatim
   *
   * Consider this 2D staggered cell and a second order central derivative scheme, this lead to
   *
   * \f$ \frac{\partial v_y}{\partial x} \f$ is calculated on position (*), so the function return the grid_key {0,0}
   *
   * \f$ \frac{\partial v_y}{\partial y} \f$ is calculated on position (0), so the function return the grid_key {-1,-1}
   *
   * \f$ \frac{\partial v_x}{\partial x} \f$ is calculated on position (0), so the function return the grid_key {-1,-1}
   *
   * \f$ \frac{\partial v_x}{\partial y} \f$ is calculated on position (*), so the function return the grid_key {0,0}
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   * \return where (in which cell grid) the derivative is calculated
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos,
						    const grid_sm<Sys_eqs::dims,void> & gs,
						    const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
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


/*! \brief Second order one sided Derivative scheme on direction i
 *
 * \verbatim
 *
 *  -1.5    2.0   -0.5
 *    +------*------*
 *
 * or
 *
 *  -0.5    2.0   -1.5
 *    *------*------+
 *
 *  in the bulk
 *
 *  -1        +1
 *   *---+---*
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class D<d,expr_type,CENTRAL_B_ONE_SIDE>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;
  
  D() {}

  D(expr_type expr_) : expr{expr_} {}
  
  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map mapping grid points
   * \param kmap position where the derivative is calculated
   * \param gs Grid info
   * \param spacing of the grid
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
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
		    typename Sys_eqs::stype coeff) const
  {
#ifdef SE_CLASS1
    if (Sys_eqs::boundary[d] == PERIODIC)
      std::cerr << __FILE__ << ":" << __LINE__ << " error, it make no sense use one sided derivation with periodic boundary, please use CENTRAL\n";
#endif

    grid_key_dx<Sys_eqs::dims> pos = g_map.getGKey(kmap);

    if (pos.get(d) == (long int)gs.size(d)-1 )
      {
	expr.value(g_map,kmap,gs,spacing,cols,1.5*coeff/spacing[d]);

	long int old_val = kmap.getKeyRef().get(d);
	kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
	expr.value(g_map,kmap,gs,spacing,cols,-2.0*coeff/spacing[d]);
	kmap.getKeyRef().set_d(d,old_val);

	old_val = kmap.getKeyRef().get(d);
	kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 2);
	expr.value(g_map,kmap,gs,spacing,cols,0.5*coeff/spacing[d]);
	kmap.getKeyRef().set_d(d,old_val);
      }
    else if (pos.get(d) == 0)
      {
	expr.value(g_map,kmap,gs,spacing,cols,-1.5*coeff/spacing[d]);

	long int old_val = kmap.getKeyRef().get(d);
	kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
	expr.value(g_map,kmap,gs,spacing,cols,2.0*coeff/spacing[d]);
	kmap.getKeyRef().set_d(d,old_val);

	old_val = kmap.getKeyRef().get(d);
	kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 2);
	expr.value(g_map,kmap,gs,spacing,cols,-0.5*coeff/spacing[d]);
	kmap.getKeyRef().set_d(d,old_val);
      }
    else
      {
	long int old_val = kmap.getKeyRef().get(d);
	kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
	expr.value(g_map,kmap,gs,spacing,cols,coeff/spacing[d]);
	kmap.getKeyRef().set_d(d,old_val);

	old_val = kmap.getKeyRef().get(d);
	kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
	expr.value(g_map,kmap,gs,spacing,cols,-coeff/spacing[d]);
	kmap.getKeyRef().set_d(d,old_val);
      }
  }

  /*! \brief Calculate the position where the derivative is calculated
   *
   * In case on non staggered case this function just return a null grid_key, in case of staggered,
   *  it calculate how the operator shift in the cell
   *
   \verbatim

   +--$--+
   |     |
   #  *  #
   |     |
   0--$--+

   # = velocity(y)
   $ = velocity(x)
   * = pressure

   \endverbatim
   *
   * In the one side stencil the cell position depend if you are or not at the boundary.
   * outside the boundary is simply the central scheme, at the boundary it is simply the
   * staggered position of the property
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   * \return where (in which cell grid) the derivative is calculated
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    return expr_type::position(pos,gs,s_pos);
  }
};


/*! \brief First order FORWARD derivative on direction i
 *
 * \verbatim
 *
 *  -1.0    1.0
 *    +------*
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class D<d,expr_type,FORWARD>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;

  D() {}

  D(expr_type expr_) : expr{expr_} {}

  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map mapping grid
   * \param kmap position where the derivative is calculated
   * \param gs Grid info
   * \param spacing grid spacing
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
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
		    typename Sys_eqs::stype coeff) const
  {

    long int old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) + 1);
    expr.value(g_map,kmap,gs,spacing,cols,coeff/spacing[d]);
    kmap.getKeyRef().set_d(d,old_val);

    // backward
    expr.value(g_map,kmap,gs,spacing,cols,-coeff/spacing[d]);
  }


  /*! \brief Calculate the position where the derivative is calculated
   *
   * In case of non staggered case this function just return a null grid_key, in case of staggered,
   * the FORWARD scheme return the position of the staggered property
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   * \return where (in which cell grid) the derivative is calculated
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos,
						    const grid_sm<Sys_eqs::dims,void> & gs,
						    const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    return expr_type::position(pos,gs,s_pos);
  }
};

/*! \brief First order BACKWARD derivative on direction i
 *
 * \verbatim
 *
 *  -1.0    1.0
 *    *------+
 *
 * \endverbatim
 *
 */
template<unsigned int d, typename expr_type>
class D<d,expr_type,BACKWARD>
{
public:

  typedef typename expr_type::sys_eqs_type Sys_eqs; /**< System of equations. */
  typedef Sys_eqs sys_eqs_type;

  expr_type expr;

  D() {}

  D(expr_type expr_) : expr{expr_} {}

  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map mapping grid
   * \param kmap position where the derivative is calculated
   * \param gs Grid info
   * \param spacing of the grid
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
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
		    typename Sys_eqs::stype coeff) const
  {

    long int old_val = kmap.getKeyRef().get(d);
    kmap.getKeyRef().set_d(d, kmap.getKeyRef().get(d) - 1);
    expr.value(g_map,kmap,gs,spacing,cols,-coeff/spacing[d]);
    kmap.getKeyRef().set_d(d,old_val);

    // forward
    expr.value(g_map,kmap,gs,spacing,cols,coeff/spacing[d]);
  }


  /*! \brief Calculate the position where the derivative is calculated
   *
   * In case of non staggered case this function just return a null grid_key, in case of staggered,
   * the BACKWARD scheme return the position of the staggered property
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   * \return where (in which cell grid) the derivative is calculated
   *
   */
  inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos, const grid_sm<Sys_eqs::dims,void> & gs, const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  {
    return expr_type::position(pos,gs,s_pos);
  }
};

#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_DERIVATIVE_HPP_ */
