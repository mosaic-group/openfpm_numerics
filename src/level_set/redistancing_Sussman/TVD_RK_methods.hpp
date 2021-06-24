//
// Created by jstark on 22.06.21.
//

#ifndef OPENFPM_NUMERICS_TVD_RK_METHODS_HPP
#define OPENFPM_NUMERICS_TVD_RK_METHODS_HPP
//void func ( void (*f)(int) ) {
//  for ( int ctr = 0 ; ctr < 5 ; ctr++ ) {
//    (*f)(ctr);
//  }
//}
#include "HelpFunctions.hpp"
#include "HelpFunctionsForGrid.hpp"
#include "FiniteDifference/Upwind_gradient.hpp"


//template <size_t Field_in, size_t Sign, size_t Gradient_out, typename gridtype>
//void get_upwind_gradient(gridtype & grid, const size_t order=5, const bool one_sided_BC=true)

template <size_t U, size_t Sign, size_t Gradient, typename grid_type, typename key_type>
void L(grid_type & grid, key_type key, size_t order, bool one_sided_BC=true)
{
	get_upwind_gradient<U, Sign, Gradient>(grid, order, one_sided_BC);
	//const double phi_n = grid.template get<Phi_0_temp>(key);
	//			const double phi_n_magnOfGrad = grid.template get<Phi_grad_temp>(key).norm();
	//			double epsilon = phi_n_magnOfGrad * spacing_x;
	//			grid.template get<Phi_nplus1_temp>(key) = get_phi_nplus1(phi_n, phi_n_magnOfGrad, time_step,
	
	double phi_magnOfGradient = grid.template get<Phi_grad_temp>(key).norm();
	double epsilon = phi_magnOfGradient * grid.spacing(0);
	return smooth_S(grid.template get<U>(key), epsilon) * (1 - phi_magnOfGradient);
}

template <size_t Un, size_t Unplus1, typename grid_type, typename T>
T tvd_runge_kutta_1_stepper(grid_type & grid, T dt, T (*L)(T))
{
	return un + dt * (*L)(un);
}


template <typename T>
T tvd_runge_kutta_1_stepper(const T & un, T dt, T (*L)(T))
{
	return un + dt * (*L)(un);
}

template <typename T>
T tvd_runge_kutta_3_stepper(const T & un, T dt, T (*L)(T))
{
	double u1 = tvd_runge_kutta_1_stepper(un, dt, &L(un));
	double u2 = 0.75 * un + 0.25 * u1 + 0.25 * dt * (*L)(u1);
	return 1.0/3.0 * un + 2.0/3.0 * u2 + 2.0/3.0 * dt * (*L)(u2);
}

#endif //OPENFPM_NUMERICS_TVD_RK_METHODS_HPP
