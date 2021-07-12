//
// Created by jstark on 27.05.21.
//
#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN  // in only one cpp file
#include <boost/test/unit_test.hpp>

// Include redistancing files
#include "util/PathsAndFiles.hpp"
#include "level_set/redistancing_Sussman/RedistancingSussman.hpp"
#include "level_set/redistancing_Sussman/NarrowBand.hpp"
// Include header files for testing
#include "Draw/DrawSphere.hpp"
#include "Draw/DrawDisk.hpp"
#include "l_norms/LNorms.hpp"
#include "analytical_SDF/AnalyticalSDF.hpp"


BOOST_AUTO_TEST_SUITE(ConvergenceTestSuite)
	
	BOOST_AUTO_TEST_CASE(RedistancingSussman_unit_sphere_convergence_test_space)
	{
		const double EPSILON = std::numeric_limits<double>::epsilon();
		
		const size_t grid_dim = 3;
		// some indices
		const size_t x = 0;
		const size_t y = 1;
		const size_t z = 2;
		
		const size_t Phi_0_grid = 0;
		const size_t SDF_sussman_grid = 1;
		const size_t SDF_exact_grid = 2;
		const size_t Error_grid = 3;
		
		for (size_t N = 32; N <= 128; N *= 2)
		{
			const double dt = 0.000165334;
			const size_t sz[grid_dim] = {N, N, N};
			const double radius = 1.0;
			const double box_lower = 0.0;
			const double box_upper = 4.0 * radius;
			Box<grid_dim, double> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(0);
			typedef aggregate<double, double, double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
			
			const double center[grid_dim] = {0.5 * (box_upper - box_lower), 0.5 * (box_upper - box_lower),
			                                 0.5 * (box_upper - box_lower)};
			init_grid_with_sphere<Phi_0_grid>(g_dist, radius, center[x], center[y],
			                                  center[z]); // Initialize sphere onto grid


//			for (int order=1; order<=5; order+=2)
			{
				Redist_options redist_options;
				redist_options.min_iter = 1e4;
				redist_options.max_iter = 1e4;

//				redist_options.order_space_op                       = order;
				redist_options.order_space_op = 1;
				redist_options.order_timestepper = 1;
				
				// set both convergence criteria to false s.t. termination only when max_iterations reached
				redist_options.convTolChange.check = false;    // define here which of the convergence criteria above should be used. If both are true, termination only occurs when both are fulfilled or when iter > max_iter
				redist_options.convTolResidual.check = false;    // (default: false)
				
				redist_options.interval_check_convergence = 100;        // interval of #iterations at which
				// convergence is checked (default: 100)
				redist_options.width_NB_in_grid_points = 8;        // width of narrow band in number of grid points. Must be at least 4, in order to have at least 2 grid points on each side of the interface. (default: 4)
				redist_options.print_current_iterChangeResidual = true;     // if true, prints out every current iteration + corresponding change from the previous iteration + residual from SDF (default: false)
				redist_options.print_steadyState_iter = true;     // if true, prints out the final iteration number when steady state was reached + final change + residual (default: true)
				
				RedistancingSussman<grid_in_type> redist_obj(g_dist,
				                                             redist_options);   // Instantiation of Sussman-redistancing class
				redist_obj.set_user_time_step(dt);
				std::cout << "dt set to = " << dt << std::endl;
				// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
				redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
				
				// Compute exact signed distance function at each grid point
				init_analytic_sdf_sphere<SDF_exact_grid>(g_dist, radius, center[x], center[y], center[z]);
				
				// Compute the absolute error between analytical and numerical solution at each grid point
				get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
				
				
				/////////////////////////////////////////////////////////////////////////////////////////////
				//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
				//	interface)
				size_t bc[grid_dim] = {PERIODIC, PERIODIC, PERIODIC};
				typedef aggregate<double> props_nb;
				typedef vector_dist<grid_dim, double, props_nb> vd_type;
				Ghost<grid_dim, double> ghost_vd(0);
				vd_type vd_narrow_band(0, box, bc, ghost_vd);
				vd_narrow_band.setPropNames({"error"});
				const size_t Error_vd = 0;
				// Compute the L_2- and L_infinity-norm and save to file
				size_t narrow_band_width = 8;
				NarrowBand<grid_in_type> narrowBand_points(g_dist,
				                                           narrow_band_width); // Instantiation of NarrowBand class
				narrowBand_points.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
				                                                                                                 vd_narrow_band);
				vd_narrow_band.write("vd_nb8p_error_N" + std::to_string(N) + "_RKorder" +
						                     std::to_string(redist_options.order_timestepper), FORMAT_BINARY);
//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
				// Compute the L_2- and L_infinity-norm and save to file
				L_norms lNorms_vd;
				lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
				write_lnorms_to_file(N, lNorms_vd,
				                     "l_norms_vd_absError_8p_order" + std::to_string(redist_options.order_timestepper),
				                     "./");

//				switch(order)
//				{
//					case 1:
//						BOOST_CHECK(lNorms_vd.l2   < 0.03369 + EPSILON);
//						BOOST_CHECK(lNorms_vd.linf < 0.06307 + EPSILON);
//						break;
//					case 3:
//						BOOST_CHECK(lNorms_vd.l2   < 0.02794 + EPSILON);
//						BOOST_CHECK(lNorms_vd.linf < 0.0586704 + EPSILON);
//						break;
//					case 5:
//						BOOST_CHECK(lNorms_vd.l2   < 0.0187199 + EPSILON);
//						BOOST_CHECK(lNorms_vd.linf < 0.0367638 + EPSILON);
//						break;
//				}
			
			}
		}
	}
	BOOST_AUTO_TEST_CASE(RedistancingSussman_unit_disk_convergence_test_time)
	{
		// CFL dt for N=64 and c=1 is 0.00100781
		const double EPSILON = std::numeric_limits<double>::epsilon();
		
		const size_t grid_dim = 2;
		// some indices
		const size_t x = 0;
		const size_t y = 1;
		
		const size_t Phi_0_grid = 0;
		const size_t SDF_sussman_grid = 1;
		const size_t SDF_exact_grid = 2;
		const size_t Error_grid = 3;
		
		size_t N = 64;
		double dt = 0.00100781;
		double tf = 1e4 * dt;
		dt *= 4.0;
		for (int i = 0; i < 5; i++)
		{
			dt /= 2.0;
			size_t max_iter = (size_t) std::round(tf / dt);
			const size_t sz[grid_dim] = {N, N};
			const double radius = 1.0;
			const double box_lower = 0.0;
			const double box_upper = 4.0 * radius;
			Box<grid_dim, double> box({box_lower, box_lower}, {box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(0);
			typedef aggregate<double, double, double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
			
			const double center[grid_dim] = {0.5 * (box_upper - box_lower), 0.5 * (box_upper - box_lower)};
			init_grid_with_disk<Phi_0_grid>(g_dist, radius, center[x], center[y]); // Initialize sphere onto grid


//			for (int order=1; order<=5; order+=2)
			{
				Redist_options redist_options;
				redist_options.min_iter = max_iter;
				redist_options.max_iter = max_iter;
				
				redist_options.order_space_op = 1;
				redist_options.order_timestepper = 1;
				
				// set both convergence criteria to false s.t. termination only when max_iterations reached
				redist_options.convTolChange.check = false;    // define here which of the convergence criteria above should be used. If both are true, termination only occurs when both are fulfilled or when iter > max_iter
				redist_options.convTolResidual.check = false;    // (default: false)
				
				redist_options.interval_check_convergence = 100;        // interval of #iterations at which
				// convergence is checked (default: 100)
				redist_options.width_NB_in_grid_points = 8;        // width of narrow band in number of grid points. Must be at least 4, in order to have at least 2 grid points on each side of the interface. (default: 4)
				redist_options.print_current_iterChangeResidual = true;     // if true, prints out every current iteration + corresponding change from the previous iteration + residual from SDF (default: false)
				redist_options.print_steadyState_iter = true;     // if true, prints out the final iteration number when steady state was reached + final change + residual (default: true)
				
				RedistancingSussman<grid_in_type> redist_obj(g_dist,
				                                             redist_options);   // Instantiation of Sussman-redistancing class
				redist_obj.set_user_time_step(dt);
				std::cout << "dt set to = " << dt << std::endl;
//				std::cout << "dt for N = " << N << " is " << redist_obj.get_time_step() << std::endl;
				// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
				redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
				
				// Compute exact signed distance function at each grid point
				init_analytic_sdf_circle<SDF_exact_grid>(g_dist, radius, center[x], center[y]);
				
				// Compute the absolute error between analytical and numerical solution at each grid point
				get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
				
				g_dist.write(
						"grid_RKorder" + std::to_string(redist_options.order_timestepper) + "_i_" + std::to_string(i),
						FORMAT_BINARY);
				/////////////////////////////////////////////////////////////////////////////////////////////
				//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
				//	interface)
				size_t bc[grid_dim] = {PERIODIC, PERIODIC};
				typedef aggregate<double> props_nb;
				typedef vector_dist<grid_dim, double, props_nb> vd_type;
				Ghost<grid_dim, double> ghost_vd(0);
				vd_type vd_narrow_band(0, box, bc, ghost_vd);
				vd_narrow_band.setPropNames({"error"});
				const size_t Error_vd = 0;
				// Compute the L_2- and L_infinity-norm and save to file
				size_t narrow_band_width = 8;
				NarrowBand<grid_in_type> narrowBand_points(g_dist,
				                                           narrow_band_width); // Instantiation of NarrowBand class
				narrowBand_points.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
				                                                                                                 vd_narrow_band);
				vd_narrow_band.write("vd_error_N" + std::to_string(N) + "_RKorder" +
						                     std::to_string(redist_options.order_timestepper), FORMAT_BINARY);
//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
				// Compute the L_2- and L_infinity-norm and save to file
				L_norms lNorms_vd;
				lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
				write_lnorms_to_file(dt, lNorms_vd,
				                     "l_norms_RKorder" + std::to_string(redist_options.order_timestepper), "./");

//				switch(order)
//				{
//					case 1:
//						BOOST_CHECK(lNorms_vd.l2   < 0.03369 + EPSILON);
//						BOOST_CHECK(lNorms_vd.linf < 0.06307 + EPSILON);
//						break;
//					case 3:
//						BOOST_CHECK(lNorms_vd.l2   < 0.02794 + EPSILON);
//						BOOST_CHECK(lNorms_vd.linf < 0.0586704 + EPSILON);
//						break;
//					case 5:
//						BOOST_CHECK(lNorms_vd.l2   < 0.0187199 + EPSILON);
//						BOOST_CHECK(lNorms_vd.linf < 0.0367638 + EPSILON);
//						break;
//				}
			
			}
		}
	}

	BOOST_AUTO_TEST_CASE(RedistancingSussman_1D_convergence_test_time)
	{
		// CFL dt in 1D for N=64 and c=1 is 0.000503905
		const double EPSILON = std::numeric_limits<double>::epsilon();
		
		std::cout << "epsilon = " << EPSILON << std::endl;
		const size_t grid_dim = 1;
		// some indices
		const size_t x = 0;
		
		const size_t Phi_0_grid = 0;
		const size_t SDF_sussman_grid = 1;
		const size_t SDF_exact_grid = 2;
		const size_t Error_grid = 3;
		
		size_t N = 64;
		double dt = 0.000503905;
		double tf = 1e4 * dt;
		dt *= 2.0;
//		for (int i = 0; i < 3; i++)
		int i = 0;
		{
			dt /= 2.0;
			size_t max_iter = (size_t) std::round(tf / dt);
			const size_t sz[grid_dim] = {N};
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower}, {box_upper});
			Ghost<grid_dim, long int> ghost(0);
			typedef aggregate<double, double, double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
			
			// Now we initialize the grid with the indicator function and the analytical solution
			auto dom = g_dist.getDomainIterator();
			while(dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> coords = g_dist.getPos(key);
				g_dist.template get<Phi_0_grid>(key) = sgn(coords.get(x));
				g_dist.template get<SDF_exact_grid>(key) = coords.get(x);
				++dom;
			}
			
			g_dist.write("grid_1D_preRedistancing", FORMAT_BINARY); // Save the grid
			
			Redist_options redist_options;
//			redist_options.min_iter = max_iter;
//			redist_options.max_iter = max_iter;
			redist_options.order_space_op    = 1;
			redist_options.order_timestepper = 1;
			// set both convergence criteria to false s.t. termination only when max_iterations reached
			redist_options.convTolChange.check = true;    // define here which of the convergence criteria above should
			// be used. If both are true, termination only occurs when both are fulfilled or when iter > max_iter
			redist_options.convTolChange.value = EPSILON;
			
			redist_options.convTolResidual.check = false;    // (default: false)
			redist_options.interval_check_convergence = 100;        // interval of #iterations at which
			// convergence is checked (default: 100)
			redist_options.width_NB_in_grid_points = 32;        // width of narrow band in number of grid points. Must
			// be at least 4, in order to have at least 2 grid points on each side of the interface. (default: 4)
			redist_options.print_current_iterChangeResidual = true;     // if true, prints out every current iteration + corresponding change from the previous iteration + residual from SDF (default: false)
			redist_options.print_steadyState_iter = true;     // if true, prints out the final iteration number when steady state was reached + final change + residual (default: true)
			
			RedistancingSussman<grid_in_type> redist_obj(g_dist,
			                                             redist_options);   // Instantiation of Sussman-redistancing class
			redist_obj.set_user_time_step(dt);
//			std::cout << "CFL dt = " << redist_obj.get_time_step() << std::endl;
			std::cout << "dt set to = " << dt << std::endl;
			// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
			redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
			
			// Compute the absolute error between analytical and numerical solution at each grid point
			get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
			g_dist.write("grid_RKorder" + std::to_string(redist_options.order_timestepper) + "_i_" + std::to_string(i),
			             FORMAT_BINARY);
			/////////////////////////////////////////////////////////////////////////////////////////////
			//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
			//	interface)
			size_t bc[grid_dim] = {NON_PERIODIC};
			typedef aggregate<double> props_nb;
			typedef vector_dist<grid_dim, double, props_nb> vd_type;
			Ghost<grid_dim, double> ghost_vd(0);
			vd_type vd_narrow_band(0, box, bc, ghost_vd);
			vd_narrow_band.setPropNames({"error"});
			const size_t Error_vd = 0;
			// Compute the L_2- and L_infinity-norm and save to file
			size_t narrow_band_width = 64;
			NarrowBand<grid_in_type> narrowBand_points(g_dist, narrow_band_width); // Instantiation of NarrowBand class
			narrowBand_points.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
			                                                                                                 vd_narrow_band);
			vd_narrow_band.write(
					"vd_error_N" + std::to_string(N) + "_RKorder" + std::to_string(redist_options.order_timestepper),
					FORMAT_BINARY);
//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
			// Compute the L_2- and L_infinity-norm and save to file
			L_norms lNorms_vd;
			lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
			write_lnorms_to_file(dt, lNorms_vd, "l_norms_RKorder" + std::to_string(redist_options.order_timestepper),
			                     "./");
			
		}
	}

BOOST_AUTO_TEST_SUITE_END()


