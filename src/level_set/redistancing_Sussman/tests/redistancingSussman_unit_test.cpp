//
// Created by jstark on 12.05.21.
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

BOOST_AUTO_TEST_SUITE(RedistancingSussmanTestSuite)
	
	BOOST_AUTO_TEST_CASE(RedistancingSussman_3D_test)
	{
		const double EPSILON = std::numeric_limits<double>::epsilon();
		const size_t grid_dim = 3;
		// some indices
		const size_t x                      = 0;
		const size_t y                      = 1;
		const size_t z                      = 2;
		
		const size_t Phi_0_grid             = 0;
		const size_t SDF_sussman_grid       = 1;
		const size_t SDF_exact_grid         = 2;
		const size_t Error_grid             = 3;
		
		size_t N = 128;
		const double dt = 0.000165334; // CFL-condition for N=128
		const size_t sz[grid_dim] = {N, N, N};
		const double radius = 1.0;
		const double box_lower = -2.0;
		const double box_upper = 2.0;
		Box<grid_dim, double> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
		Ghost<grid_dim, long int> ghost(0);
		typedef aggregate<double, double, double, double> props;
		typedef grid_dist_id<grid_dim, double, props > grid_in_type;
		grid_in_type g_dist(sz, box, ghost);
		g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
		
		const double center[grid_dim] = {0.5*(box_upper+box_lower),
										 0.5*(box_upper+box_lower),
										 0.5*(box_upper+box_lower)};
		init_grid_with_sphere<Phi_0_grid>(g_dist, radius, center[x], center[y], center[z]); // Initialize sphere onto grid
		
		
		Redist_options redist_options;
		redist_options.min_iter                             = 1e4;
		redist_options.max_iter                             = 1e4;
		
		redist_options.convTolChange.check                  = false;
		redist_options.convTolResidual.check                = false;
		
		redist_options.interval_check_convergence           = 1e3;
		redist_options.width_NB_in_grid_points              = 8;
		redist_options.print_current_iterChangeResidual     = true;
		redist_options.print_steadyState_iter               = true;
		
		RedistancingSussman<grid_in_type> redist_obj(g_dist, redist_options);   // Instantiation of Sussman-redistancing class
		std::cout << "New CFL timestep = " << redist_obj.get_time_step() << std::endl;
//		redist_obj.set_user_time_step(dt);
//		std::cout << "dt set to = " << dt << std::endl;
		// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
		redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
		
		double finalIter     = redist_obj.get_finalIteration();
		double finalChange   = redist_obj.get_finalChange();
		double finalResidual = redist_obj.get_finalResidual();
		double finalNumberPointsInNarrowBand = redist_obj.get_finalNumberNbPoints();
		
		
		// Compute exact signed distance function at each grid point
		init_analytic_sdf_sphere<SDF_exact_grid>(g_dist, radius, center[x], center[y], center[z]);
		
		// Compute the absolute error between analytical and numerical solution at each grid point
		get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
		//	interface)
		size_t bc[grid_dim] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
		typedef aggregate<double> props_nb;
		typedef vector_dist<grid_dim, double, props_nb> vd_type;
		Ghost<grid_dim, double> ghost_vd(0);
		vd_type vd_narrow_band(0, box, bc, ghost_vd);
		vd_narrow_band.setPropNames({"error"});
		size_t narrow_band_width = 8;
		NarrowBand<grid_in_type> narrowBand(g_dist, narrow_band_width); // Instantiation of NarrowBand class
		const size_t Error_vd = 0;
		narrowBand.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
		                                                                                          vd_narrow_band);
		// Compute the L_2- and L_infinity-norm and save to file
		L_norms lNorms_vd;
		lNorms_vd = get_l_norms_vector<0>(vd_narrow_band);
		std::cout << lNorms_vd.l2 << ", " << lNorms_vd.linf << std::endl;
		
		// l-norms original implementation of 1st order upwinding:
		// 32,0.033643,0.063065; dt = 0.000165334; 1e4 iterations
		// 0.0336846, 0.0630651
		// new impl
		// order 1: 0.0336846, 0.0630651
		// order 3: 0.02794, 0.0586704
		// order 5: 0.0187199, 0.0367638
		
		
		BOOST_CHECK(lNorms_vd.l2 <   0.03369 + EPSILON);
		BOOST_CHECK(lNorms_vd.linf < 0.06307 + EPSILON);
		
	}
	
	BOOST_AUTO_TEST_CASE(RedistancingSussman_2D_test)
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
		
		const double center[grid_dim] = {0.5 * (box_upper + box_lower), 0.5 * (box_upper + box_lower)};
		init_grid_with_disk<Phi_0_grid>(g_dist, radius, center[x], center[y]); // Initialize disk onto grid
		
		
		Redist_options redist_options;
		redist_options.min_iter = 1e3;
		redist_options.max_iter = 1e4;
		
		redist_options.convTolChange.check = false;
//		redist_options.convTolChange.value = 1e-3;
		redist_options.convTolResidual.check = false;
		
		redist_options.interval_check_convergence = 1e3;
		redist_options.width_NB_in_grid_points = 8;
		redist_options.print_current_iterChangeResidual = true;
		redist_options.print_steadyState_iter = true;
		redist_options.save_temp_grid = true;
		
		RedistancingSussman<grid_in_type> redist_obj(g_dist, redist_options);
		
		redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
		int finalIter        = redist_obj.get_finalIteration();
		double finalChange   = redist_obj.get_finalChange();
		double finalResidual = redist_obj.get_finalResidual();
		int finalNumberPointsInNarrowBand = redist_obj.get_finalNumberNbPoints();
		
		std::cout << "Final distance printed from variables:"
		<< finalIter << ", "
		<< to_string_with_precision(finalChange, 15) << ", "
		<< to_string_with_precision(finalResidual, 15) << ", "
		<< finalNumberPointsInNarrowBand << std::endl;
		
		
		// Compute exact signed distance function at each grid point
		init_analytic_sdf_circle<SDF_exact_grid>(g_dist, radius, center[x], center[y]);
		
		// Compute the absolute error between analytical and numerical solution at each grid point
		get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
		
//		g_dist.write("grid_original_sdf_N" + std::to_string(N), FORMAT_BINARY);
		/////////////////////////////////////////////////////////////////////////////////////////////
		//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
		//	interface)
		size_t bc[grid_dim] = {NON_PERIODIC, NON_PERIODIC};
		typedef aggregate<double, double, double> props_nb;
		typedef vector_dist<grid_dim, double, props_nb> vd_type;
		Ghost<grid_dim, double> ghost_vd(0);
		vd_type vd_narrow_band(0, box, bc, ghost_vd);
		vd_narrow_band.setPropNames({"SDF_sussman", "SDF_exact", "Error"});
		
		const size_t SDF_sussman_vd = 0;
		const size_t dPhi_vd        = 1;
//		const size_t Error_vd       = 2;
		// Compute the L_2- and L_infinity-norm and save to file
		size_t narrow_band_width = 8;
		NarrowBand<grid_in_type> narrowBand_points(g_dist, narrow_band_width); // Instantiation of NarrowBand class
		narrowBand_points.get_narrow_band_copy_three_scalar_properties<SDF_sussman_grid,
		SDF_sussman_grid, SDF_exact_grid, Error_grid,
		0, 1, 2>(g_dist, vd_narrow_band);
		
//		narrowBand_points.get_narrow_band<SDF_sussman_grid, SDF_sussman_vd, dPhi_vd>(g_dist, vd_narrow_band);
		
		vd_narrow_band.write("vd_corrected_withGradients_N" + std::to_string(N), FORMAT_BINARY);
//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
		// Compute the L_2- and L_infinity-norm and save to file
//		L_norms lNorms_vd;
//		lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
//		write_lnorms_to_file(N, lNorms_vd, "l_norms_testoriginal", "./");
	}
	
	BOOST_AUTO_TEST_CASE(RedistancingSussman_1D_test)
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
		
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(0);
		typedef aggregate<double, double, double, double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_in_type;
	
		size_t N = 64;
		{
			const size_t sz[grid_dim] = {N};
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
			
			{
				Redist_options redist_options;
				redist_options.min_iter = 1e3;
				redist_options.max_iter = 1e16;
				
				redist_options.convTolChange.check = true;
				redist_options.convTolChange.value = EPSILON;
				
				redist_options.convTolResidual.check = false;
				redist_options.interval_check_convergence = 1e3;
				redist_options.width_NB_in_grid_points = 10;
				redist_options.print_current_iterChangeResidual = true;
				redist_options.print_steadyState_iter = true;
				
				RedistancingSussman<grid_in_type> redist_obj(g_dist,
				                                             redist_options);
				redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
				
				// Compute the absolute error between analytical and numerical solution at each grid point
				get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
				g_dist.write("grid_1D_subcellfix", FORMAT_BINARY);
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
				size_t narrow_band_width = redist_options.width_NB_in_grid_points - 2;
				NarrowBand<grid_in_type> narrowBand_points(g_dist, narrow_band_width); // Instantiation of NarrowBand class
				narrowBand_points.get_narrow_band_copy_specific_property<SDF_sussman_grid, Error_grid, Error_vd>(g_dist,
				                                                                                                 vd_narrow_band);
				vd_narrow_band.write("vd_error_1D_subcellfix_N" + std::to_string(N), FORMAT_BINARY);
	//				vd_narrow_band.save("test_data/output/vd_nb8p_error" + std::to_string(N) + ".bin");
				// Compute the L_2- and L_infinity-norm and save to file
				L_norms lNorms_vd;
				lNorms_vd = get_l_norms_vector<Error_vd>(vd_narrow_band);
				write_lnorms_to_file(N, lNorms_vd, "l_norms_", "./");
			}
			
		}
	}
	
BOOST_AUTO_TEST_SUITE_END()

