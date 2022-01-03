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
#include "l_norms/LNorms.hpp"
#include "analytical_SDF/AnalyticalSDF.hpp"

BOOST_AUTO_TEST_SUITE(RedistancingSussmanTestSuite)
	
	BOOST_AUTO_TEST_CASE(RedistancingSussman_unit_sphere)
	{
		typedef double phi_type;
		typedef double space_type;
		const phi_type EPSILON = std::numeric_limits<phi_type>::epsilon();
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
		const space_type dt = 0.000165334; // CFL-condition for N=128
		const size_t sz[grid_dim] = {N, N, N};
		const space_type radius = 1.0;
		const space_type box_lower = -2.0;
		const space_type box_upper = 2.0;
		Box<grid_dim, space_type> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
		Ghost<grid_dim, long int> ghost(0);
		typedef aggregate<phi_type, phi_type, phi_type, phi_type> props;
		typedef grid_dist_id<grid_dim, space_type, props > grid_in_type;
		grid_in_type g_dist(sz, box, ghost);
		g_dist.setPropNames({"Phi_0", "SDF_sussman", "SDF_exact", "Relative error"});
		
		const space_type center[grid_dim] = {0.5*(box_upper+box_lower),
										 0.5*(box_upper+box_lower),
										 0.5*(box_upper+box_lower)};
		init_grid_with_sphere<Phi_0_grid>(g_dist, radius, center[x], center[y], center[z]); // Initialize sphere onto grid
		
		
		Redist_options<phi_type> redist_options;
		redist_options.min_iter                             = 1e4;
		redist_options.max_iter                             = 1e4;
		
		redist_options.convTolChange.check                  = false;
		redist_options.convTolResidual.check                = false;
		
		redist_options.interval_check_convergence           = 1e3;
		redist_options.width_NB_in_grid_points              = 8;
		redist_options.print_current_iterChangeResidual     = true;
		redist_options.print_steadyState_iter               = true;
		
		RedistancingSussman<grid_in_type, phi_type> redist_obj(g_dist, redist_options);   // Instantiation of
		// Sussman-redistancing class
		std::cout << "New CFL timestep = " << redist_obj.get_time_step() << std::endl;
//		redist_obj.set_user_time_step(dt);
//		std::cout << "dt set to = " << dt << std::endl;
		// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
		redist_obj.run_redistancing<Phi_0_grid, SDF_sussman_grid>();
		
		// Compute exact signed distance function at each grid point
		init_analytic_sdf_sphere<SDF_exact_grid>(g_dist, radius, center[x], center[y], center[z]);
		
		// Compute the absolute error between analytical and numerical solution at each grid point
		get_absolute_error<SDF_sussman_grid, SDF_exact_grid, Error_grid>(g_dist);
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		//	Get narrow band: Place particles on interface (narrow band width e.g. 4 grid points on each side of the
		//	interface)
		size_t bc[grid_dim] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
		typedef aggregate<phi_type> props_nb;
		typedef vector_dist<grid_dim, space_type, props_nb> vd_type;
		Ghost<grid_dim, space_type> ghost_vd(0);
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

BOOST_AUTO_TEST_SUITE_END()


