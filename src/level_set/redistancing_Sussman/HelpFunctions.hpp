//
// Created by jstark on 2019-11-08.
//
/**
 * @file HelpFunctions.hpp
 *
 * @brief Header file containing small mathematical help-functions that are needed e.g. for the Sussman redistancing.
 *
 * @author Justina Stark
 * @date November 2019
 */
#ifndef REDISTANCING_SUSSMAN_HELPFUNCTIONS_HPP
#define REDISTANCING_SUSSMAN_HELPFUNCTIONS_HPP

#include <iostream>
#include <fstream>
/**@brief Gets the sign of a variable.
 *
 * @tparam T Inferred type of input variable.
 * @param val Scalar variable of arbitrary type for which the sign should be determined.
 * @return Integer variable that contains a -1, if argument negative, 0 if argument=0, and +1 if argument positive.
 */
template <class T>
int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}
/**@brief Gets the smoothed sign of a variable.
 *
 * @details Sign function with smoothing effect on numerical solution (see: <a href="https://www.math.uci
 * .edu/~zhao/publication/mypapers/ps/local_levelset.ps">Peng \a et \a al. "A PDE-Based Fast Local
 * Level Set Method"</a>, equation (36). Peng \a et \a al.: "The choice of approximation to S(d) by (36) solves the
 * problem of the changing of sign of Phi (thus moving the interface across the cell boundary) in the reinitialization
 * step when Phi is steep and speeds up the convergence when Phi is flat at the interface."
 *
 * @f[ sgn_{\epsilon}(\phi) = \frac{\phi}{ \sqrt{\phi^2 + |\nabla\phi|^2 \Delta x^2} }@f]
 *
 * where the \a &phi; is the approximation of the SDF from the last iteration.
 *
 * @tparam T Inferred type of the variable for which the smooth sign should be determined.
 * @param val Scalar variable for which the smooth sign should be determined.
 * @param epsilon Scalar variable containing the smoothing factor. In case of Peng: @f[\epsilon =
 * |\nabla\phi|^2 \Delta x^2@f]
 * @return Smooth sign of the argument variable: @f[sgn_{\epsilon}(\phi)@f].
 */
template <typename T>
T smooth_S(T val, T epsilon)
{
	return (val / sqrt(val * val + epsilon * epsilon));
}


/**@brief Checks, if two values are sufficiently close to each other within a given tolerance range, as to be
 * considered as approximately equal.
 *
 * @tparam T Inferred type of the two variables, for which it should be checked, if they are sufficiently close.
 * @param val1 Variable that contains the first value.
 * @param val2 Variable that contains the second value.
 * @param tolerance Tolerance (or error) by which values are allowed to deviate while still being considered as equal.
 * @return True, if \p val1 and \p val2 are the same +/-  \p tolerance. False, if they differ by more
 * than \p tolerance.
 */
template <class T>
bool isApproxEqual(T val1, T val2, T tolerance)
{
	return (val1 <= val2 + tolerance && val1 >= val2 - tolerance);
}

/**@brief Appends the value of a given variable of any type to a textfile as string.
 *
 * @tparam T Inferred type of the input variable, whose value should be appended to textfile as string.
 * @param textfile Std::string that contains the path and filename of the textfile, to which \p value should be added.
 * @param value Variable that contains the value which should be appended to the \p textfile.
 */
template <typename T>
void append_value_to_textfile(std::string & textfile, T value)
{
	std::ofstream out(textfile);
	out << value;
}

/**@brief Converts value into string maintaining a desired precision.
 *
 * @tparam T Template type of vlaue.
 * @param myValue Value of type T.
 * @param n Number of digits after the point the string of the value should have
 * @return String containing myValue with precision n.
 */
template <typename T>
std::string to_string_with_precision(const T myValue, const size_t n = 6)
{
	std::ostringstream out;
	out.precision(n);
	out << std::fixed << myValue;
	return out.str();
}


#endif //REDISTANCING_SUSSMAN_HELPFUNCTIONS_HPP
