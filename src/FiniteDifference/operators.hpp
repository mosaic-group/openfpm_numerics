/*
 * operators.hpp
 *
 *  Created on: Dec 05, 2019
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_HPP_

#include "FiniteDifference/Derivative.hpp"
#include "FiniteDifference/Laplacian.hpp"
#include "FiniteDifference/Average.hpp"
#include "FiniteDifference/sum.hpp"
#include "FiniteDifference/mul.hpp"
#include "eq.hpp"


/**
 * @fn operator+
 * @brief Creates an object of the type sum<>.
 * tparam expr1_type Type of the first expression to add.
 * tparam expr2_type Type of the second expression to add.
 * tparam sfinae Extra parameter set to a default value "typename expr1_type::sys_eqs_type" that makes only some type of expression to match this operator.
 * param[in] expr1 First expression to add.
 * param[in] expr2 Second expression to add.
 * @return sum<> object of the two operands.
 */
template<typename expr1_type, typename expr2_type, typename sfinae = typename expr1_type::sys_eqs_type>
sum<expr1_type, expr2_type> operator+(const expr1_type & expr1, const expr2_type & expr2) {
  return sum<expr1_type, expr2_type>{expr1, expr2};
}

/**
 * @fn operator-
 * @brief Creates an object of the type sum<expr1,minus<expr2>>.
 * tparam expr1_type Type of the first expression for the minus operator.
 * tparam expr2_type Type of the expression to substract.
 * tparam sfinae Extra parameter set to a default value "typename expr1_type::sys_eqs_type" that makes only some type of expression to match this operator.
 * param[in] expr1 Expression to add.
 * param[in] expr2 Expression to substract.
 * @return sum<expr1,minus<expr2>> object of the two operands.
 */
template<typename expr1_type, typename expr2_type, typename sfinae = typename expr1_type::sys_eqs_type>
sum<expr1_type, minus<expr2_type>> operator-(const expr1_type & expr1, const expr2_type & expr2) {
  auto minus_obj = minus<expr2_type>{expr2};
  return sum<expr1_type, minus<expr2_type>>(expr1, minus_obj);
}

/**
 * @fn operator*
 * @brief Creates an object of the type mul<>.
 * tparam expr1_type Type of the first expression to multiply.
 * tparam expr2_type Type of the second expression to multiply.
 * tparam sfinae Extra parameter set to a default value "typename expr1_type::sys_eqs_type" that makes only some type of expression to match this operator.
 * param[in] expr1 First expression to multiply.
 * param[in] expr2 Second expression to multiply.
 * @return mul<> object of the two operands.
 */
template<typename expr1_type, typename expr2_type, typename sfinae = typename expr1_type::sys_eqs_type>
mul<expr1_type, expr2_type> operator*(const expr1_type & expr1, const expr2_type & expr2) {
  return mul<expr1_type, expr2_type>{expr1, expr2};
}

/**
 * @class Der
 * @brief Creates an derivative object to be used with any kind of expression.
 * tparam d Direction of the derivative.
 * tparam expr_type Type of the expression to derive.
 * tparam impl Type of derivative: CENTRAL, BACKWARD, FORWARD.
 */
template<unsigned int d, typename Sys_eqs, unsigned int impl=CENTRAL>
class Der {

public:
  
  typedef Sys_eqs sys_eqs_type; /**< Extra type. Used for recognition of "valid" expressions. */

  /**
   * @fn Der()
   * @brief Default constructor.
   */
  Der() {}

  /**
   * @fn operator() (const expr_type &)
   * @brief Constructs a D (derivative) object.
   * param[in] expr_ Expression to derive.
   * return D<> (derivative) object.
   */
  template<typename expr_type>
  D<d,expr_type,impl> operator() (const expr_type expr_) {
    return D<d,expr_type,impl>(expr_);
  }
};


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_HPP_ */
