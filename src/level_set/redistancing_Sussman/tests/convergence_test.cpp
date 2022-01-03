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
#include "l_norms/LNorms.hpp"
#include "analytical_SDF/AnalyticalSDF.hpp"



BOOST_AUTO_TEST_SUITE(ConvergenceTestSuite)
	typedef double phi_type;
	typedef double space_type;
	const phi_type EPSILON = std::numeric_limits<phi_type>::epsilon();
	
	space_type get_min_required_time_step(size_t Nmax, space_type b_low, space_type b_up, int dims)
	{
		space_type dx = (b_up - b_low) / (space_type) Nmax;
		space_type sum = 0;
		for (int d = 0; d < dims; d++)
		{
			sum += 1 / (dx * dx);
		}
		return 0.5 / sum;
	}

	BOOST_AUTO_TEST_CASE(RedistancingSussman_convergence_test_1D)
	{
		// CFL dt in 1D for N=64 and c=1 is 0.000503905
		
		const size_t grid_dim = 1;
		// some indices
		const size_t x = 0;
		
		const size_t Phi_0_grid = 0;
		const size_t SDF_sussman_grid = 1;
		const size_t SDF_exact_grid = 2;
		const size_t Error_grid = 3;
		
		const space_type box_lower = -1.0;
		const space_type box_upper = 1.0;
		Box<grid_dim, space_type> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(0);
		typedef aggregate<phi_type, phi_type, phi_type, phi_type> props;
		typedef grid_dist_id<grid_dim, space_type, props> grid_in_type;
		
//		space_type dt = 0.000503905;
		size_t Nmin = 32;
		size_t Nmax = 4096;
		space_type dt = get_min_required_time_step(Nmax, box_lower, box_upper, grid_dim);
		for (size_t N = Nmin; N <= Nmax; N*=2)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
			
			// Now we initialize the grid with the indicator function and the analytical solution
			auto dom = g_dist.getDomainIterator();
			while(dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, space_type> coords = g_dist.getPos(key);
				g_dist.template get<Phi_0_grid>(key) = sgn(coords.get(x));
				g_dist.template get<SDF_exact_grid>(key) = coords.get(x);
				++dom;
			}
			
			{
				Redist_options<phi_type> redist_options;
				redist_options.min_iter = 1e3;
				redist_options.max_iter = 1e16;
				
				// set both convergence criteria to false s.t. termination only when max_iterations reached
				redist_options.convTolChange.check = true;    // define here which of the convergence criteria above should
				// be used. If both are true, termination only occurs when both are fulfilled or when iter > max_iter
				redist_options.convTolChange.value = EPSILON;
				
				redist_options.convTolResidual.check = false;    // (default: false)
				redist_options.interval_check_convergence = 1000;        // interval of #iterations at which
				// convergence is checked (default: 100)
				redist_options.width_NB_in_grid_points = 10;        // width of narrow band in number of grid points.
				// Must
				// be at least 4, in order to have at least 2 grid points on each side of the interface. (default: 4)
				redist_options.print_current_iterChangeResidual = true;     // if true, prints out every current iteration + corresponding change from the previous iteration + residual from SDF (default: false)
				redist_options.print_steadyState_iter = true;     // if true, prints out the final iteration number when steady state was reached + final change + residual (default: true)
				
				RedistancingSussman<grid_in_type, phi_type> redist_obj(g_dist,
				                                             redist_options);   // Instantiation of Sussman-redistancing class
				std::cout << "CFL dt = " << redist_obj.get_time_step() << std::endl;
				redist_obj.set_user_time_step(dt);
				std::cout << "dt set to = " << dt << std::endl;
				// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
				redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
				
				// Compute the absolute error between analytical and numerical solution at each grid point
				get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
//				g_dist.write("grid_postRedistancing", FORMAT_BINARY);
				/////////////////////////////////////////////////////////////////////////////////////////////
				//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
				//	interface)
				size_t bc[grid_dim] = {NON_PERIODIC};
				typedef aggregate<phi_type> props_nb;
				typedef vector_dist<grid_dim, space_type, props_nb> vd_type;
				Ghost<grid_dim, space_type> ghost_vd(0);
				vd_type vd_narrow_band(0, box, bc, ghost_vd);
				vd_narrow_band.setPropNames({"error"});
				const size_t Error_vd = 0;
				// Compute the L_2- and L_infinity-norm and save to file
				size_t narrow_band_width = redist_options.width_NB_in_grid_points - 2;
				NarrowBand<grid_in_type> narrowBand_points(g_dist, narrow_band_width); // Instantiation of NarrowBand class
				narrowBand_points.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
				                                                                                                 vd_narrow_band);
//			    vd_narrow_band.write("vd_error_N" + std::to_string(N), FORMAT_BINARY);
//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
				// Compute the L_2- and L_infinity-norm and save to file
				L_norms lNorms_vd;
				lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
				write_lnorms_to_file(N, lNorms_vd, "l_norms_", "./");
			}
			
		}
	}
	
	BOOST_AUTO_TEST_CASE(RedistancingSussman_convergence_test_3D)
	{
		const size_t grid_dim = 3;
		// some indices
		const size_t x                      = 0;
		const size_t y                      = 1;
		const size_t z                      = 2;
		
		const size_t Phi_0_grid             = 0;
		const size_t SDF_sussman_grid       = 1;
		const size_t SDF_exact_grid         = 2;
		const size_t Error_grid             = 3;
		
		for (size_t N=32; N <=128; N*=2)
		{
			const space_type dt = 0.000165334;
			const size_t sz[grid_dim] = {N, N, N};
			const space_type radius = 1.0;
			const space_type box_lower = 0.0;
			const space_type box_upper = 4.0 * radius;
			Box<grid_dim, space_type> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(0);
			typedef aggregate<phi_type, phi_type, phi_type, phi_type> props;
			typedef grid_dist_id<grid_dim, space_type, props > grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
			
			const space_type center[grid_dim] = {0.5*(box_upper+box_lower), 0.5*(box_upper+box_lower), 0.5*(box_upper+box_lower)};
			init_grid_with_sphere<Phi_0_grid>(g_dist, radius, center[x], center[y], center[z]); // Initialize sphere onto grid
			
			int order = 1;
			{
				Redist_options<phi_type> redist_options;
				redist_options.min_iter                             = 1e4;
				redist_options.max_iter                             = 1e4;
				
				// set both convergence criteria to false s.t. termination only when max_iterations reached
				redist_options.convTolChange.check                  = false;    // define here which of the convergence criteria above should be used. If both are true, termination only occurs when both are fulfilled or when iter > max_iter
				redist_options.convTolResidual.check                = false;    // (default: false)
				
				redist_options.interval_check_convergence           = 100;        // interval of #iterations at which
				// convergence is checked (default: 100)
				redist_options.width_NB_in_grid_points              = 8;        // width of narrow band in number of grid points. Must be at least 4, in order to have at least 2 grid points on each side of the interface. (default: 4)
				redist_options.print_current_iterChangeResidual     = true;     // if true, prints out every current iteration + corresponding change from the previous iteration + residual from SDF (default: false)
				redist_options.print_steadyState_iter               = true;     // if true, prints out the final iteration number when steady state was reached + final change + residual (default: true)
				
				RedistancingSussman<grid_in_type, phi_type> redist_obj(g_dist, redist_options);   // Instantiation of
				// Sussman-redistancing class
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
				typedef aggregate<phi_type> props_nb;
				typedef vector_dist<grid_dim, space_type, props_nb> vd_type;
				Ghost<grid_dim, space_type> ghost_vd(0);
				vd_type vd_narrow_band(0, box, bc, ghost_vd);
				vd_narrow_band.setPropNames({"error"});
				const size_t Error_vd = 0;
				// Compute the L_2- and L_infinity-norm and save to file
				size_t narrow_band_width = 8;
				NarrowBand<grid_in_type> narrowBand_points(g_dist, narrow_band_width); // Instantiation of NarrowBand class
				narrowBand_points.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
				                                                                                                 vd_narrow_band);
//				vd_narrow_band.write("vd_nb8p_error_N" + std::to_string(N) + "_order" +
//				std::to_string(order), FORMAT_BINARY);
//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
				// Compute the L_2- and L_infinity-norm and save to file
				L_norms lNorms_vd;
				lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
				write_lnorms_to_file(N, lNorms_vd, "l_norms_vd_absError_8p_order" + std::to_string(order), "./");
				
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
	
BOOST_AUTO_TEST_SUITE_END()


