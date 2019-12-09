/*
 * operators.hpp
 *
 *  Created on: Dec 05, 2019 (re-written)
 *      Author: amfoggia
 */

#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_HPP_

/**
 * @fn operator+
 * @brief Creates an object of the type sum<>
 * @return sum<> object of the two operands.
 */
sum<expr1_type, expr2_type, expr1_type::sys_nn_type> operator+(expr1_type & expr1, expr2_type & expr2) {
  return sum<expr1_type, expr2_type, expr1_type::sys_nn_type>{expr1, expr2};
}

/**
 * @fn operator-
 * @brief Creates an object of the type sum<expr1,minus<expr2>>
 * @return sum<expr1,minus<expr2>> object of the two operands.
 */
sum<expr1_type, minus<expr2_type>, expr1_type::sys_nn_type> operator+(expr1_type & expr1, expr2_type & expr2) {
  auto minus_obj = minus<expr2_type, expr1_type::sys_nn_type>{expr2};
  return sum<expr1_type, minus<expr2_type, expr1_type::sys_nn_type>, expr1_type::sys_nn_type>(expr1, minus_obj);
}

/**
 * @fn operator*
 * @brief Creates an object of the type mul<>
 * @return mul<> object of the two operands.
 */
mul<expr1_type, expr2_type, expr1_type::sys_nn_type> operator*(expr1_type & expr1, expr2_type & expr2) {
  return mul<expr1_type, expr2_type, expr1_type::sys_nn_type>{expr1, expr2};
}


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_HPP_ */
