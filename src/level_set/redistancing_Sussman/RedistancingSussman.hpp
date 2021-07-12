//
// Created by jstark on 2020-04-06.
//
/**
 * @file RedistancingSussman.hpp
 * @class RedistancingSussman
 *
 * @brief Class for reinitializing a level-set function into a signed distance function using Sussman redistancing.
 *
 * @details First-time-only redistancing step (see: <a href="https://www.researchgate
 * .net/publication/2642502_An_Efficient_Interface_Preserving_Level_Set_Re
 * -Distancing_Algorithm_And_Its_Application_To_Interfacial_Incompressible_Fluid_Flow.html">M. Sussman and E. Fatemi,
 * “Efficient, interface-preserving level set redistancing algorithm and its application to interfacial
 * incompressible fluid flow” (1999)</a> ). In the Sussman-redistancing, a level-set-function Phi_0, which can be a
 * simple step-function having positive values inside the object and negative values outside, is reinitialized to a
 * signed distance function Phi_SDF by finding the steady-state solution of the following PDE:
 * @f[ \phi_{t} + sgn_{\epsilon}(\phi)(|\nabla\phi| - 1) = 0  @f]
 *
 * The signed distance function is defined as:
 *
 * @f[ \phi_{\text{SDF}} = \begin{cases}
 *    +d & \text{orthogonal distance to the closest surface of a point lying inside the object} \\
 *     0 & \text{object surface} \\
 *    -d & \text{orthogonal distance to the closest surface of a point lying outside the object} \\
 * \end{cases} @f]
 *
 * @author Justina Stark
 * @date April 2020
 */

#ifndef REDISTANCING_SUSSMAN_REDISTANCINGSUSSMAN_HPP
#define REDISTANCING_SUSSMAN_REDISTANCINGSUSSMAN_HPP

// Include standard library header files
#include <iostream>
#include <string>
#include <fstream>
// Include OpenFPM header files
#include "Vector/vector_dist.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"

// Include other header files
#include "HelpFunctions.hpp"
#include "HelpFunctionsForGrid.hpp"
//#include "ComputeGradient.hpp"
#include "FiniteDifference/Upwind_gradient.hpp"


/** @brief Optional convergence criterium checking the total change.
 *
 * @details The change between the current and the previous iterations is computed as sum over all the If
 * Conv_tol_change
 * .check = true, the total change of Phi
 * between the iterations is considered and the
 * redistancing finishes when the Conv_tol_change.value (or the max-iteration) is reached.
 */
struct Conv_tol_change
{
	bool check = true; 	///< If true, the total change of Phi (see DistFromSol::change)
	///< between the iterations is considered and the redistancing finishes when the Conv_tol_change.value (or the
	///< max-iteration) is reached. If false, the change is not considered as steady-state
	///< criterion. Default: true.
	double value = 1e-6; 	///< Double variable that stores the tolerance value for the change at which the iterative
	///< redistancing process is considered as steady-state. (see also DistFromSol::change)
};

/** @brief Optional convergence criterium checking the residual.
 *
 * @details If Conv_tol_residual.check = true, the residual, that is abs(magnitude gradient of phi - 1), is considered
 * and the redistancing finishes when the Conv_tol_residual.value (or the max-iteration) is reached.
 */
struct Conv_tol_residual
{
	bool check = false; ///< If true, the residual of Phi (see DistFromSol::residual) is considered and the
	///< redistancing finishes when the Conv_tol_residual.value (or the
	///< max-iteration) is reached. If false, the change is not considered as steady-state criterion. Default: false.
	double value = 1e-1; 		///< Double variable that stores the tolerance value for the residual at which the
	///< iterative redistancing process is considered as steady-state. (see also
	///< DistFromSol::residual)
};

/** @brief Structure to bundle options for redistancing.
 * @struct Redist_options
 * @details For the redistancing, we can choose some options. These options will then be passed bundled as a structure to
 * the redistancing function. Setting these options is optional, since they all have a Default value as well.
 *
 * @param min_iter: Minimum number of iterations before steady state in narrow band will be checked (Default: 100).
 * @param max_iter: Maximum number of iterations you want to run the redistancing, even if steady state might not yet
 *                have been reached (Default: 1e6).
 * @param order_space_op: Order of accuracy of the upwind gradient computation when solving the eikonal equation
 *                        during the redistancing. Options are: {1, 3, 5}. Default: 5;
 * @param convTolChange.value: Convolution tolerance for the normalized total change of Phi in the narrow band between
 *                           two consecutive iterations (Default: 1e-6).
 * @param convTolChange.check: Set true, if you want to use the normalized total change between two iterations as
 *                           measure of how close you are to the steady state solution. Redistancing will then stop
 *                           if convTolChange.value is reached or if the current iteration is bigger than max_iter.
 * @param convTolResidual.value: Convolution tolerance for the residual, that is abs(magnitude gradient of phi - 1) of
 *                             Phi in the narrow band (Default 1e-1).
 * @param convTolResidual.check: Set true, if you want to use the residual of the current iteration as measure of how
 *                             close you are to the steady state solution. Redistancing will then stop if
 *                             convTolResidual.value is reached or if the current iteration is bigger than max_iter.
 * @param interval_check_convergence: Interval of number of iterations at which convergence to steady state is checked
 *                                  (Default: 100).
 * @param width_NB_in_grid_points: Width of narrow band in number of grid points. Must be at least 4, in order to
 *                               have at least 2 grid points on each side of the interface. Is automatically set
 *                               to 4, if a value smaller than 4 is chosen (Default: 4).
 * @param print_current_iterChangeResidual: If true, the number of the current iteration, the corresponding change
 *                                        w.r.t the previous iteration and the residual is printed (Default: false).
 * @param print_steadyState_iter: If true, the number of the steady-state-iteration, the corresponding change
 *                              w.r.t the previous iteration and the residual is printed (Default: false).
 */
struct Redist_options
{
	size_t min_iter = 1000;
	size_t max_iter = 1e6;
	
	size_t order_space_op = 5;
	size_t order_timestepper = 3;
	
	Conv_tol_change convTolChange;
	Conv_tol_residual convTolResidual;
	
	size_t interval_check_convergence = 100;
	size_t width_NB_in_grid_points = 4;
	bool print_current_iterChangeResidual = false;
	bool print_steadyState_iter = true;
	bool save_temp_grid = false;
};

/** @brief Bundles total residual and total change over all the grid points.
 *
 * @details In order to evaluate, how close the current solution is away from the convergence criterium.
 *
 * @see RedistancingSussman::get_residual_and_change_NB()
 */
struct DistFromSol
{
	double residual; 	///< Double variable that contains the absolute value of how far the gradient magnitude of
	///< the current \a &phi; of iteration number \a i is away from being equal to 1: @f[ abs(|\nabla\phi_i| - 1 ) @f]
	///< It is computed for all grid points that lie within the narrow band and
	///< normalized to the number of grid points that lie in that narrow band.
	double change;   	///< Double variable that contains the absolute value of the change of \a &phi; between the
	///< current
	///< iteration \a i and the previous iteration \a i-1: @f[ abs( \phi_i - \phi_{i-1} ) @f] It is
	///< computed for all grid
	///< points that lie within the narrow band and normalized to the number of grid points that
	///< lie in that narrow band.
};


/**@brief Class for reinitializing a level-set function into a signed distance function using Sussman redistancing.
 * @file RedistancingSussman.hpp
 * @class RedistancingSussman
 * @tparam grid_in_type Inferred type of input grid, which stores the initial level-set function Phi_0.
 */
template <typename grid_in_type>
class RedistancingSussman
{
public:
	/** @brief Constructor initializing the redistancing options, the temporary internal grid and reference variable
	 * to the input grid.
	 *
	 * @param grid_in Input grid with min. 2 properties: 1.) Phi_0, 2.) Phi_SDF <- will be overwritten with
	 * re-distancing result
	 * @param redistOptions User defined options for the Sussman redistancing process
	 *
	 */
	RedistancingSussman(grid_in_type &grid_in, Redist_options &redistOptions) : redistOptions(redistOptions),
	                                                                            r_grid_in(grid_in),
	                                                                            g_temp(grid_in.getDecomposition(),
	                                                                                   grid_in.getGridInfoVoid().getSize(),
	                                                                                   Ghost<grid_in_type::dims, long int>(3))
	{
		time_step = get_time_step_CFL(g_temp);
#ifdef SE_CLASS1
		if(!(redistOptions.order_timestepper == 1 || redistOptions.order_timestepper == 3))
		{
			std::cout << "You set the order of time discretization to " << redistOptions.order_timestepper <<
						 " but time stepping only implemented for order 1 and 3. Using order 3 instead..." << std::endl;
			redistOptions.order_timestepper = 3;
		}
		if(!(redistOptions.order_space_op == 1
		|| redistOptions.order_space_op == 3
		|| redistOptions.order_space_op == 5))
		{
			std::cout << "You set the order of space discretization to " << redistOptions.order_space_op <<
					" but space discretization only implemented for order 1, 3 and 5. Using order 5 instead..." <<
					std::endl;
			redistOptions.order_space_op = 5;
		}
#endif // SE_CLASS1
	}
	
	/**@brief Aggregated properties for the temporary grid.
	 *
	 * @details The initial (input) Phi_0 (will be updated by Phi_{n+1} after each redistancing step),
	 * Phi_{n+1} (received from redistancing),
	 * gradient of Phi_{n+1},
	 * sign of the original input Phi_0 (for the upwinding).
	 */
	typedef aggregate<double, double, Point<grid_in_type::dims, double>, int, double, double> props_temp;
	/** @brief Type definition for the temporary grid.
	 */
	typedef grid_dist_id<grid_in_type::dims, typename grid_in_type::stype, props_temp> g_temp_type;
	/**
	 * @brief Create temporary grid, which is only used inside the class for the redistancing.
	 *
	 * @details The temporary grid stores the following 4 properties:
	 * the initial (input) Phi_0 (will be updated by Phi_{n+1} after each redistancing step),
	 * Phi_{n+1}(received from redistancing),
	 * gradient of Phi_{n+1},
	 * sign of the original input Phi_0 (for the upwinding).
	 */
	g_temp_type g_temp;
	
	/**@brief Runs the Sussman-redistancing.
	 *
	 * @details Copies Phi_0 from input grid to an internal temporary grid which allows
	 * having more properties. Computes the gradients. Runs the redistancing on the internal temporary grid. Copies
	 * resulting signed distance function to the Phi_SDF_out property of the input grid.
	 */
	template<size_t Phi_0_in, size_t Phi_SDF_out> void run_redistancing()
	{
		init_temp_grid<Phi_0_in>();
		init_sign_prop<Phi_0, Sign_Phi0in>(
				g_temp); // initialize Sign_Phi0in with the sign of the initial (pre-redistancing) Phi_0
		// Get initial gradients
		get_upwind_gradient<Phi_0, Sign_Phi0in, Phi_grad>(g_temp, redistOptions.order_space_op, true);
		std::cout << "timestep before redistancing starts: " << time_step << std::endl;
		iterative_redistancing(g_temp); // Do the redistancing on the temporary grid
		copy_gridTogrid<Phi_nplus1, Phi_SDF_out>(g_temp, r_grid_in); // Copy resulting SDF function to input grid
	}
	
	/** @brief Overwrite the time_step found via CFL condition with an individual time_step.
	 *
	 * @details If the user wants to overwrite the time_step found via CFL condition with an individual time_step.
	 * Should be only used carefully, time_step must not be too large (jump over solution) nor too small (extremely
	 * slow).
	 *
	 * @param dt Artificial time step by which re-distancing should be performed.
	 *
	 */
	template<typename T>
	void set_user_time_step(T dt)
	{
		time_step = dt;
	}
	/** @brief Access the artificial timestep (private member) which will be used for the iterative redistancing.
	 * @see get_time_step_CFL(g_temp_type &grid), set_user_time_step()
	 */
	double get_time_step()
	{
		/// This timestep is computed according to the grid spacing fulfilling the CFL condition.
		return time_step;
	}

private:
	//	Some indices for better readability
	static constexpr size_t Phi_0       = 0; ///< Property index of Phi_0 on the temporary grid.
	static constexpr size_t Phi_nplus1  = 1; ///< Property index of Phi_n+1 on the temporary grid.
	static constexpr size_t Phi_grad    = 2; ///< Property index of gradient of Phi_n on the temporary grid.
	static constexpr size_t Sign_Phi0in = 3; ///< Property index of sign of initial (input) Phi_0 (temp. grid).
	static constexpr size_t L_factor    = 4; ///< Property index of RHS: sign(phi)*(1-|gradPhi|)
	static constexpr size_t Phi_interm  = 5; ///< Property index of intermediate state for RK3.
	
	
	//	Member variables
	Redist_options redistOptions; ///< Instantiate redistancing options.
	grid_in_type &r_grid_in; ///< Define reference to input grid.
	
	double h_max = get_biggest_spacing(g_temp); ///< Grid spacing in less resolved direction.
	double h_min = get_smallest_spacing(g_temp); ///< Grid spacing in higher resolved direction.
	/// Transform the half-bandwidth in no_of_grid_points into physical half-bandwidth kappa.
	double kappa = ceil(redistOptions.width_NB_in_grid_points / 2.0) * h_max;
	/**@brief Artificial timestep for the redistancing iterations.
	 * @see get_time_step_CFL(g_temp_type &grid), get_time_step(), set_user_time_step()
	 */
	double time_step;
	
	//	Member functions
	/** @brief Copies values from input grid to internal temporary grid and initializes ghost layer with minimum
	 * value of input grid.
	 */
	template<size_t Phi_0_in>
	void init_temp_grid()
	{
		double min_value = get_min_val<Phi_0_in>(r_grid_in); // get minimum Phi_0 value on the input grid
		init_grid_and_ghost<Phi_0>(g_temp, min_value); // init. Phi_0 (incl. ghost) with min. Phi_0
		init_grid_and_ghost<Phi_nplus1>(g_temp, min_value); // init. Phi_nplus1 (incl. ghost) with min. Phi_0
		copy_gridTogrid<Phi_0_in, Phi_0>(r_grid_in, g_temp); // Copy Phi_0 from the input grid to Phi_0
	}
	
	template <size_t U, size_t Sign, size_t Gradient, size_t L>
	void get_L()
	{
		g_temp.template ghost_get<U, Sign>();
		get_upwind_gradient<U, Sign, Gradient>(g_temp, redistOptions.order_space_op, true);
		g_temp.template ghost_get<Gradient>(KEEP_PROPERTIES);
		double spacing_x = g_temp.getSpacing()[0];
		auto dom = g_temp.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			const double phi_n = g_temp.template get<U>(key);
			const double phi_n_magnOfGrad = g_temp.template get<Gradient>(key).norm();
			double epsilon = phi_n_magnOfGrad * spacing_x;
			g_temp.template get<L>(key) = smooth_S(phi_n, epsilon) * (1 - phi_n_magnOfGrad);
			++dom;
		}
	}
	
	template <size_t Un, size_t L, size_t U1>
	void get_u1(double dt)
	{
		g_temp.template ghost_get<Un, L>();
		auto dom = g_temp.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			g_temp.template get<U1>(key) = g_temp.template get<Un>(key) + dt * g_temp.template get<L>(key);
			++dom;
		}
	}
	
	template <size_t Un, size_t U1, size_t L, size_t U2>
	void get_u2(double dt)
	{
		g_temp.template ghost_get<Un, U1, L>();
		auto dom = g_temp.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			g_temp.template get<U2>(key) =
			                0.75 * g_temp.template get<Un>(key) + 
	                        0.25 * g_temp.template get<U1>(key) +
                            0.25 * dt * g_temp.template get<L>(key);
			++dom;
		}
	}
	
	template <size_t Un, size_t U2, size_t L, size_t Unplus1>
	void get_u3(double dt)
	{
		g_temp.template ghost_get<Un, U2, L>();
		auto dom = g_temp.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			g_temp.template get<Unplus1>(key) =
							1.0/3.0 * g_temp.template get<Un>(key) +
							2.0/3.0 * g_temp.template get<U2>(key) +
					        2.0/3.0 * dt * g_temp.template get<L>(key);
			++dom;
		}
	}
	
	template <size_t Un, size_t L, size_t Sign, size_t Gradient, size_t Unplus1>
	void tvd_runge_kutta_1_stepper(double dt)
	{
		get_L<Un, Sign, Gradient, L>();
		get_u1<Un, L, Unplus1>(dt);
	}
	
	template <size_t Un, size_t Uintermed, size_t L, size_t Sign, size_t Gradient, size_t Unplus1>
	void tvd_runge_kutta_3_stepper(double dt)
	{
		get_L<Un, Sign, Gradient, L>();
		get_u1<Un, L, Uintermed>(dt);
		get_L<Uintermed, Sign, Gradient, L>();
		//template <size_t Un, size_t U1, size_t L, size_t U2>
		get_u2<Un, Uintermed, L, Uintermed>(dt);
		get_L<Uintermed, Sign, Gradient, L>();
		get_u3<Un, Uintermed, L, Unplus1>(dt);
	}
	

	
	/** @brief Go one re-distancing time-step on the whole grid.
    *
    * @param grid Internal temporary grid.
    */
	void go_one_redistancing_step_whole_grid(g_temp_type &grid)
	{
		switch(redistOptions.order_timestepper)
		{
			case 1:
				//	template <size_t Un, size_t L, size_t Sign, size_t Gradient, size_t Unplus1>
				tvd_runge_kutta_1_stepper<Phi_0, L_factor, Sign_Phi0in, Phi_grad, Phi_nplus1>(time_step);
				break;
			case 3:
			default:
				// 	template <size_t Un, size_t U1, size_t L, size_t Sign, size_t Gradient, size_t Unplus1>
				tvd_runge_kutta_3_stepper<Phi_0, Phi_interm, L_factor, Sign_Phi0in, Phi_grad, Phi_nplus1>(time_step);
				break;
		}

	}
	
	/** @brief Updates Phi_n with the new Phi_n+1 and recomputes the gradients.
	 *
	 * @param grid Internal temporary grid.
	 */
	void update_grid(g_temp_type &grid)
	{
		copy_gridTogrid<Phi_nplus1, Phi_0>(grid, grid); // Update Phi_0
		get_upwind_gradient<Phi_0, Sign_Phi0in, Phi_grad>(grid, redistOptions.order_space_op, true);
	}
	
	/** @brief Checks if a node lays within the narrow band around the interface.
	 *
	 * @param Phi Value of Phi at that specific node.
	 *
	 * @return True, if node lays within nb., false, if the distance to the interface is > kappa.
	 */
	bool lays_inside_NB(double Phi)
	{
		return (abs(Phi) <= kappa);
	}
	
	/** @brief Checks how far current solution is from fulfilling the user-defined convergence criteria.
	 *
	 * @param grid Internal temporary grid.
	 *
	 * @return Total residual (1 - phi_gradient_magnitude) and total change from the last time step,
	 * both normalized by the number of grid nodes in the narrow band.
	 */
	DistFromSol get_residual_and_change_NB(g_temp_type &grid)
	{
		double total_residual = 0;
		double total_change = 0;
		double total_points_in_nb = 0;
		auto dom = grid.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			if (lays_inside_NB(grid.template get<Phi_nplus1>(key)))
			{
				total_points_in_nb += 1.0;
				auto dphi_magn = grid.template get<Phi_grad>(key).norm();
				total_residual += abs(dphi_magn - 1);
				total_change += abs(grid.template get<Phi_nplus1>(key) - grid.template get<Phi_0>(key));
			}
			++dom;
		}
		auto &v_cl = create_vcluster();
		v_cl.sum(total_points_in_nb);
		v_cl.sum(total_residual);
		v_cl.sum(total_change);
		v_cl.execute();
		return {total_residual / total_points_in_nb, total_change / total_points_in_nb};
	}
	
	/** @brief Prints out the iteration number, residual and change of the current re-distancing iteration
	 *
	 * @param grid Internal temporary grid.
	 * @param iter Current re-distancing iteration.
	 *
	 * @return Total residual (1 - phi_gradient_magnitude) and total change from the last time step,
	 * both normalized by the number of grid nodes in the narrow band.
	 */
	void print_out_iteration_change_residual(g_temp_type &grid, size_t iter)
	{
		DistFromSol distFromSol = get_residual_and_change_NB(grid);
		auto &v_cl = create_vcluster();
		if (v_cl.rank() == 0)
		{
			if (iter == 0)
			{
				std::cout << "Iteration,Change,Residual" << std::endl;
			}
			std::cout << iter << "," << distFromSol.change << "," << distFromSol.residual << std::endl;
		}
	}
	
	/** @brief Checks steady-state is reached in the narrow band.
	 *
	 * @details Checks if change and/or residual between 2 iterations smaller than user-defined convolution tolerance.
	 *
	 * @param grid Internal temporary grid.
	 *
	 * @return True, if steady-state reached, else false.
	 *
	 * @see Conv_tol_change, Conv_tol_residual
	 */
	bool steady_state_NB(g_temp_type &grid)
	{
		bool steady_state = false;
		DistFromSol distFromSol = get_residual_and_change_NB(grid);
		if (redistOptions.convTolChange.check && redistOptions.convTolResidual.check)
		{
			steady_state = (distFromSol.change <= redistOptions.convTolChange.value &&
					distFromSol.residual <= redistOptions.convTolResidual.value);
		}
		else
		{
			if (redistOptions.convTolChange.check)
			{
				steady_state = (distFromSol.change <= redistOptions.convTolChange.value);
			}       // Use the normalized total change between two iterations in the narrow bands steady-state criterion
			if (redistOptions.convTolResidual.check)
			{
				steady_state = (distFromSol.residual <= redistOptions.convTolResidual.value);
			} // Use the normalized total residual of phi compared to SDF in the narrow bands steady-state criterion
		}
		return steady_state;
	}
	
	/** @brief Runs Sussman re-distancing on the internal temporary grid.
	 *
	 * @details The number of iterations is minimum redistOptions.min_iter iterations and finishes when either
	 * steady-state or redist_options.max_iter is reached. The steady-state convergence accuracy depends on the
	 * user defined redist_options.convTolChange.value and redist_options.convTolResidual.value, respectively.
	 *
	 * @param grid Internal temporary grid.
	 *
	 */
	void iterative_redistancing(g_temp_type &grid)
	{
		for (size_t i = 0; i <= redistOptions.max_iter; i++)
		{
			go_one_redistancing_step_whole_grid(grid);
			
			if (i % redistOptions.interval_check_convergence == 0) // after some iterations check if steady state
				// is reached in the narrow band
			{
				if (redistOptions.print_current_iterChangeResidual)
				{
					print_out_iteration_change_residual(grid, i);
				}
				
				if (i >= redistOptions.min_iter)
				{
					if (steady_state_NB(grid))
					{
						if (redistOptions.print_steadyState_iter)
						{
							auto &v_cl = create_vcluster();
							if (v_cl.rank() == 0)
							{
								std::cout << "Steady state reached at iteration: " << i << std::endl;
								std::cout << "Final_Iteration,Change,Residual" << std::endl;
							}
							print_out_iteration_change_residual(grid, i);
						}
						update_grid(grid); // Update Phi
						if (redistOptions.save_temp_grid)
						{
							g_temp.setPropNames({"Phi_0", "Phi_nplus1", "Phi_grad", "Sign_Phi0in"});
							g_temp.save("g_temp_redistancing.hdf5"); // HDF5 file}
						}
						break;
					}
				}
			}
			update_grid(grid);
		}
		// If save_temp_grid set true, save the temporary grid to an hdf5 file that can be reloaded onto a grid and
		// reused
		if (redistOptions.save_temp_grid)
		{
			g_temp.setPropNames({"Phi_0", "Phi_nplus1", "Phi_grad", "Sign_Phi0in"});
			g_temp.save("g_temp_redistancing.hdf5"); // HDF5 file}
		}
	}
};

#endif //REDISTANCING_SUSSMAN_REDISTANCINGSUSSMAN_HPP
