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

//! Multiplication expression
// template<typename v_expr>
// struct const_mul_functor_value
// {
//   //! Number of elements in the vector v_expr
//   typedef boost::mpl::size<v_expr> size;

//   //! Last element of sum
//   typedef typename boost::mpl::at<v_expr,boost::mpl::int_<size::value-1> >::type last;

//   //! sum functor
//   std::unordered_map<long int,typename last::stype> & cols;

//   //! grid size
//   const grid_sm<last::dims,void> & gs;

//   //! grid mapping
//   const grid_dist_id<last::dims,typename last::stype,aggregate<size_t>,typename last::b_grid::decomposition::extended_type> & g_map;

//   //! grid position
//   grid_dist_key_dx<last::dims> & kmap;

//   //! coefficent
//   typename last::stype coeff;

//   //! spacing
//   typename last::stype (& spacing)[last::dims];

//   /*! \brief constructor
//    *
//    * \param g_map mapping grid
//    * \param kmap grid point (row) where we evaluate the non-zero colums
//    * \param gs grid size
//    * \param spacing grid spacing
//    * \param cols non zero colums
//    * \param coeff multiplication coefficent
//    *
//    */
//   const_mul_functor_value(const grid_dist_id<last::dims,typename last::stype,aggregate<size_t>,typename last::b_grid::decomposition::extended_type> & g_map,
// 			  grid_dist_key_dx<last::dims> & kmap,
// 			  const grid_sm<last::dims,void> & gs,
// 			  typename last::stype (& spacing)[last::dims],
// 			  std::unordered_map<long int,typename last::stype> & cols,
// 			  typename last::stype coeff)
//     :cols(cols),gs(gs),g_map(g_map),kmap(kmap),coeff(coeff),spacing(spacing)
//   {};



//   //! It call this function for every constant expression in the mul
//   template<typename T>
//   void operator()(T& t)
//   {
//     typedef typename boost::mpl::at<v_expr, boost::mpl::int_<T::value> >::type cfield;

//     coeff *= has_val<is_const_field<cfield>::value * 1,cfield>::call_val();
//   }

//   /*! \brief Get the coefficent
//    *
//    * \return the coefficent
//    *
//    */
//   typename last::stype getCoeff()
//   {
//     return coeff;
//   }
// };

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
  mul(expr1_type & _expr1, expr2_type & _expr2) : expr1{_expr1}, expr2{_expr2} {}
  
  /*! \brief Calculate which colums of the Matrix are non zero
   *
   * \param g_map mapping grid
   * \param kmap position where the multiplication is calculated
   * \param gs Grid info
   * \param spacing of the grid
   * \param cols non-zero colums calculated by the function
   * \param coeff coefficent (constant in front of the derivative)
   *
   */
  inline void value(const grid_dist_id<Sys_eqs::dims,typename Sys_eqs::stype,aggregate<size_t>,typename Sys_eqs::b_grid::decomposition::extended_type> & g_map,
		    grid_dist_key_dx<Sys_eqs::dims> & kmap,
		    const grid_sm<Sys_eqs::dims,void> & gs,
		    typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
		    std::unordered_map<long int,typename Sys_eqs::stype > & cols,
		    typename Sys_eqs::stype coeff) const
  {
    coeff *= expr1.value();
    expr2.value(g_map,kmap,gs,spacing,cols,coeff);

    // const_mul_functor_value<v_expr> mfv(g_map,kmap,gs,spacing,cols,coeff);

    // //
    // boost::mpl::for_each_ref< boost::mpl::range_c<int,0,v_sz::type::value - 2> >(mfv);

    // //! Last element of multiplication
    // typedef typename boost::mpl::at< v_expr ,boost::mpl::int_<v_sz::value-2> >::type last_m;

    // last_m::value(g_map,kmap,gs,spacing,cols,mfv.coeff);
  }

  /*! \brief Calculate the position in the cell where the mul operator is performed
   *
   * it just return the position of the staggered property in the last expression
   *
   * \param pos position where we are calculating the derivative
   * \param gs Grid info
   * \param s_pos staggered position of the properties
   *
   */
  // inline static grid_key_dx<Sys_eqs::dims> position(grid_key_dx<Sys_eqs::dims> & pos,
  // 						    const grid_sm<Sys_eqs::dims,void> & gs,
  // 						    const comb<Sys_eqs::dims> (& s_pos)[Sys_eqs::nvar])
  // {
  //   return boost::mpl::at<v_expr, boost::mpl::int_<v_sz::type::value - 2> >::type::position(pos,gs,s_pos);
  // }
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_MUL_HPP_ */
