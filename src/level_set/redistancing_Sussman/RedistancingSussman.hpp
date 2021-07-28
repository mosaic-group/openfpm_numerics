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
	double value = 1e-15; 	///< Double variable that stores the tolerance value for the change at which the iterative
	///< redistancing process is considered as steady-state. (see also DistFromSol::change)
};

/** @brief Optional convergence criterium checking the residual.
 *
 * @details If Conv_tol_residual.check = true, the residual, that is abs(magnitude gradient of phi - 1), is considered
 * and the redistancing finishes when the Conv_tol_residual.value (or the max-iteration) is reached.
 */
struct Conv_tol_residual
{
	bool check = true; ///< If true, the residual of Phi (see DistFromSol::residual) is considered and the
	///< redistancing finishes when the Conv_tol_residual.value (or the
	///< max-iteration) is reached. If false, the change is not considered as steady-state criterion. Default: true.
	double value = 1e-3; 		///< Double variable that stores the tolerance value for the residual at which the
	///< iterative redistancing process is considered as steady-state. (see also
	///< DistFromSol::residual)
};

/** @brief Structure to bundle options for redistancing.
 * @struct Redist_options
 * @details For the redistancing, we can choose some options. These options will then be passed bundled as a structure to
 * the redistancing function. Setting these options is optional, since they all have a Default value as well.
 *
 * @param min_iter: Minimum number of iterations before steady state in narrow band will be checked (Default: 1e5).
 * @param max_iter: Maximum number of iterations you want to run the redistancing, even if steady state might not yet
 *                have been reached (Default: 1e12).
 * @param order_space_op: Order of accuracy of the upwind gradient computation when solving the eikonal equation
 *                        during the redistancing. Options are: {1, 3, 5} (Default: 5).
 * @param convTolChange.value: Convolution tolerance for the normalized total change of Phi in the narrow band between
 *                           two consecutive iterations (Default: 1e-15).
 * @param convTolChange.check: Set true, if you want to use the normalized total change between two iterations as
 *                           measure of how close you are to the steady state solution. Redistancing will then stop
 *                           if convTolChange.value is reached or if the current iteration is bigger than max_iter.
 * @param convTolResidual.value: Convolution tolerance for the residual, that is abs(magnitude gradient of phi - 1) of
 *                             Phi in the narrow band (Default 1e-3).
 * @param convTolResidual.check: Set true, if you want to use the residual of the current iteration as measure of how
 *                             close you are to the steady state solution. Redistancing will then stop if
 *                             convTolResidual.value is reached or if the current iteration is bigger than max_iter.
 * @param interval_check_convergence: Interval of number of iterations at which convergence to steady state is checked
 *                                  (Default: 100).
 * @param width_NB_in_grid_points: Width of narrow band in number of grid points. Convergence is checked for this
 *                                  area around the interface only, so don't choose too small! (Default: 8).
 * @param print_current_iterChangeResidual: If true, the number of the current iteration, the corresponding change
 *                                        w.r.t the previous iteration and the residual is printed (Default: false).
 * @param print_steadyState_iter: If true, the number of the steady-state-iteration, the corresponding change
 *                              w.r.t the previous iteration and the residual is printed (Default: false).
 * @param save_temp_grid: If true, save the temporary grid as hdf5 that can be reloaded onto a grid
 
 */
struct Redist_options
{
	size_t min_iter = 1e5;
	size_t max_iter = 1e12;
	
	Conv_tol_change convTolChange;
	Conv_tol_residual convTolResidual;
	
	size_t interval_check_convergence = 100;
	size_t width_NB_in_grid_points = 8;
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
	double change;   	///< Double variable that contains the absolute value of the change of \a &phi; between the
	///< current
	///< iteration \a i and the previous iteration \a i-1: @f[ abs( \phi_i - \phi_{i-1} ) @f] It is
	///< computed for all grid
	///< points that lie within the narrow band and normalized to the number of grid points that
	///< lie in that narrow band.
	double residual; 	///< Double variable that contains the absolute value of how far the gradient magnitude of
	///< the current \a &phi; of iteration number \a i is away from being equal to 1: @f[ abs(|\nabla\phi_i| - 1 ) @f]
	///< It is computed for all grid points that lie within the narrow band and
	///< normalized to the number of grid points that lie in that narrow band.
	int count; ///< Integer variable that contains the number of points that could be assigned to the narrow band.
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
		time_step = get_time_step_CFL(grid_in);
		order_upwind_gradient = 1;
#ifdef SE_CLASS1
		assure_minimal_thickness_of_NB();
#endif // SE_CLASS1
	}
	
	/**@brief Aggregated properties for the temporary grid.
	 *
	 * @details The initial (input) Phi_0 (will be updated by Phi_{n+1} after each redistancing step),
	 * Phi_{n+1} (received from redistancing),
	 * gradient of Phi_{n+1},
	 * sign of the original input Phi_0 (for the upwinding).
	 */
	typedef aggregate<double, Point<grid_in_type::dims, double>, int>
	        props_temp;
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
		init_sign_prop<Phi_n_temp, Phi_0_sign_temp>(
				g_temp); // initialize Phi_0_sign_temp with the sign of the initial (pre-redistancing) Phi_0
		// Get initial gradients
		get_upwind_gradient<Phi_n_temp, Phi_0_sign_temp, Phi_grad_temp>(g_temp, order_upwind_gradient, true);
		iterative_redistancing(g_temp); // Do the redistancing on the temporary grid
		copy_gridTogrid<Phi_n_temp, Phi_SDF_out>(g_temp, r_grid_in); // Copy resulting SDF function to input grid
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
	
	int get_finalIteration()
	{
		return final_iter;
	}
	
	double get_finalChange()
	{
		return distFromSol.change;
	}
	
	double get_finalResidual()
	{
		return distFromSol.residual;
	}
	
	int get_finalNumberNbPoints()
	{
		return distFromSol.count;
	}
	
private:
	//	Some indices for better readability
	static constexpr size_t Phi_n_temp          = 0; ///< Property index of Phi_0 on the temporary grid.
	static constexpr size_t Phi_grad_temp       = 1; ///< Property index of gradient of Phi_n on the temporary grid.
	static constexpr size_t Phi_0_sign_temp     = 2; ///< Property index of sign of initial (input) Phi_0 (temp. grid).
	
	//	Member variables
	Redist_options redistOptions; ///< Instantiate redistancing options.
	grid_in_type &r_grid_in; ///< Define reference to input grid.
	
	DistFromSol distFromSol; ///< Instantiate distance from solution in terms of change, residual, numb. point in NB.
	int final_iter = 0; ///< Will be set to the final iteration when redistancing ends.
	
	/// Transform the half-bandwidth in no_of_grid_points into physical half-bandwidth kappa.
	double kappa = ceil(redistOptions.width_NB_in_grid_points / 2.0) * get_biggest_spacing(g_temp);
	/**@brief Artificial timestep for the redistancing iterations.
	 * @see get_time_step_CFL(g_temp_type &grid), get_time_step(), set_user_time_step()
	 */
	double time_step;
	int order_upwind_gradient;
	//	Member functions
#ifdef SE_CLASS1
	/** @brief Checks if narrow band thickness >= 4 grid points. Else, sets it to 4 grid points.
		*
		* @details Makes sure, that the narrow band within which the convergence criteria are checked during the
		* redistancing, is thick enough.
   */
	void assure_minimal_thickness_of_NB()
	{
		if (redistOptions.width_NB_in_grid_points < 8)
		{
			std::cout << "The narrow band thickness that you chose for the convergence check is very small. Consider "
			             "setting redist_options.width_NB_in_grid_points to at least 8" << std::endl;
		} // check narrow band width of convergence check if defined SE_CLASS1
	}
#endif // SE_CLASS1
	/** @brief Copies values from input grid to internal temporary grid and initializes ghost layer with minimum
	 * value of input grid.
	 */
	template<size_t Phi_0_in>
	void init_temp_grid()
	{
		double min_value = get_min_val<Phi_0_in>(r_grid_in); // get minimum Phi_0 value on the input grid
		init_grid_and_ghost<Phi_n_temp>(g_temp, min_value); // init. Phi_n_temp (incl. ghost) with min. Phi_0
		copy_gridTogrid<Phi_0_in, Phi_n_temp>(r_grid_in, g_temp); // Copy Phi_0 from the input grid to Phi_n_temp
	}

	/** @brief Run one timestep of re-distancing and compute Phi_n+1.
	 *
	 * @param phi_n Phi value on current node and current time.
	 * @param phi_n_magnOfGrad Gradient magnitude of current Phi from upwinding FD.
	 * @param dt Time step.
	 * @param sgn_phi_n Sign of the current Phi, should be the smooth sign.
	 *
	 * @return Phi_n+1 which is the Phi of the next time step on current node.
	 *
	 */
	double get_phi_nplus1(double phi_n, double phi_n_magnOfGrad, double dt, double sgn_phi_n)
	{
		return phi_n + dt * sgn_phi_n * (1 - phi_n_magnOfGrad);
	}
	
	/** @brief Go one re-distancing time-step on the whole grid.
    *
    * @param grid Internal temporary grid.
    */
	void go_one_redistancing_step_whole_grid(g_temp_type &grid)
	{
		get_upwind_gradient<Phi_n_temp, Phi_0_sign_temp, Phi_grad_temp>(grid, order_upwind_gradient, true);
		grid.template ghost_get<Phi_n_temp, Phi_grad_temp>();
		auto dom = grid.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			const double phi_n = grid.template get<Phi_n_temp>(key);
			const double phi_n_magnOfGrad = grid.template get<Phi_grad_temp>(key).norm();
			double epsilon = phi_n_magnOfGrad * grid.getSpacing()[0];
			grid.template get<Phi_n_temp>(key) = get_phi_nplus1(phi_n, phi_n_magnOfGrad, time_step,
			                                                         smooth_S(phi_n, epsilon));
			++dom;
		}
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
	
	/** @brief Re-computes the member variables distFromSol.change, distFromSol.residual, distFromSol.count for the
	 * Phi of the current iteration. Needed to check how far current solution is from fulfilling the user-defined convergence criteria.
	 *
	 * @param grid Internal temporary grid.
	 */
	void update_distFromSol(g_temp_type &grid)
	{
		double max_residual = 0;
		double max_change = 0;
		int count = 0;
		auto dom = grid.getDomainIterator();
		while (dom.isNext())
		{
			auto key = dom.get();
			if (lays_inside_NB(grid.template get<Phi_n_temp>(key)))
			{
				count++;
				double phi_n_magnOfGrad = grid.template get<Phi_grad_temp>(key).norm();
				double epsilon = phi_n_magnOfGrad * grid.getSpacing()[0];
				double phi_nplus1 = get_phi_nplus1(grid.template get<Phi_n_temp>(key), phi_n_magnOfGrad, time_step,
				                                   smooth_S(grid.template get<Phi_n_temp>(key), epsilon));
				
				if (abs(phi_nplus1 - grid.template get<Phi_n_temp>(key)) > max_change)
				{
					max_change = abs(phi_nplus1 - grid.template get<Phi_n_temp>(key));
				}
				
				if (abs(phi_n_magnOfGrad - 1) > max_residual) { max_residual = abs(phi_n_magnOfGrad - 1); }
			}
			++dom;
		}
		auto &v_cl = create_vcluster();
		v_cl.max(max_change);
		v_cl.max(max_residual);
		v_cl.sum(count);
		v_cl.execute();
		
		// Update member variable distFromSol
		distFromSol.change   = max_change;
		distFromSol.residual = max_residual;
		distFromSol.count    = count;
	}
	
	/** @brief Prints out the iteration number, max. change, max. residual and number of points in the narrow band of
	 * the current re-distancing iteration.
	 *
	 * @param grid Internal temporary grid.
	 * @param iter Current re-distancing iteration.
	 */
	void print_out_iteration_change_residual(g_temp_type &grid, size_t iter)
	{
		update_distFromSol(grid);
		auto &v_cl = create_vcluster();
		if (v_cl.rank() == 0)
		{
			if (iter == 0)
			{
				std::cout << "Iteration,MaxChange,MaxResidual,NumberOfNarrowBandPoints" << std::endl;
			}
			std::cout << iter
			<< "," << to_string_with_precision(distFromSol.change, 15)
			<< "," << to_string_with_precision(distFromSol.residual, 15)
			<< "," << distFromSol.count << std::endl;
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
		update_distFromSol(grid);
		if (redistOptions.convTolChange.check && redistOptions.convTolResidual.check)
		{
			steady_state = (
					distFromSol.change <= redistOptions.convTolChange.value &&
					distFromSol.residual <= redistOptions.convTolResidual.value &&
					distFromSol.count > 0
					);
		}
		else
		{
			if (redistOptions.convTolChange.check)
			{
				steady_state = (distFromSol.change <= redistOptions.convTolChange.value && distFromSol.count > 0);
			}       // Use the normalized total change between two iterations in the narrow bands steady-state criterion
			if (redistOptions.convTolResidual.check)
			{
				steady_state = (distFromSol.residual <= redistOptions.convTolResidual.value && distFromSol.count > 0);
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
		int i = 0;
		while (i < redistOptions.max_iter)
		{
			for (int j = 0; j < redistOptions.interval_check_convergence; j++)
			{
				go_one_redistancing_step_whole_grid(grid);
				++i;
			}
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
							std::cout << "Steady state criterion reached at iteration: " << i << std::endl;
							std::cout << "FinalIteration,MaxChange,MaxResidual" << std::endl;
						}
						print_out_iteration_change_residual(grid, i);
					}
					break;
				}
			}
		}
		update_distFromSol(grid);
		final_iter = i;
		// If save_temp_grid set true, save the temporary grid as hdf5 that can be reloaded onto a grid
		if (redistOptions.save_temp_grid)
		{
			get_upwind_gradient<Phi_n_temp, Phi_0_sign_temp, Phi_grad_temp>(g_temp, order_upwind_gradient, true);
			g_temp.setPropNames({"Phi_n", "Phi_grad_temp", "Phi_0_sign_temp"});
			g_temp.save("g_temp_redistancing.hdf5"); // HDF5 file}
		}
	}
};

#endif //REDISTANCING_SUSSMAN_REDISTANCINGSUSSMAN_HPP
