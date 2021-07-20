//
// Created by jstark on 10.06.21.
//
#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN  // in only one cpp file
#include <boost/test/unit_test.hpp>

// Include redistancing files
#include "level_set/redistancing_Sussman/NarrowBand.hpp"
// Include header files for testing
#include "Draw/DrawSphere.hpp"
#include "analytical_SDF/AnalyticalSDF.hpp"

BOOST_AUTO_TEST_SUITE(NarrowBandTestSuite)
	BOOST_AUTO_TEST_CASE(NarrowBand_unit_sphere)
	{
		const size_t dims = 3;
		// some indices
		const size_t x                      = 0;
		const size_t y                      = 1;
		const size_t z                      = 2;

		const size_t Phi_0_grid             = 0;
		const size_t SDF_exact_grid         = 1;

		const size_t SDF_vd                 = 0;
		const size_t Gradient_vd            = 1;
		const size_t magnOfGrad_vd          = 2;


		size_t N = 32;
		const double dt = 0.000165334;
		const size_t sz[dims] = {N, N, N};
		const double radius = 1.0;
		const double box_lower = 0.0;
		const double box_upper = 4.0 * radius;
		Box<dims, double> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
		Ghost<dims, long int> ghost(0);
		typedef aggregate<double, double> props;
		typedef grid_dist_id<dims, double, props > grid_in_type;
		grid_in_type g_dist(sz, box, ghost);
		g_dist.setPropNames({"Phi_0", "SDF_exact"});

		const double center[dims] = {0.5*(box_upper+box_lower), 0.5*(box_upper+box_lower), 0.5*(box_upper+box_lower)};
		init_grid_with_sphere<Phi_0_grid>(g_dist, radius, center[x], center[y], center[z]); // Initialize sphere onto grid
		// Compute exact signed distance function at each grid point
		init_analytic_sdf_sphere<SDF_exact_grid>(g_dist, radius, center[x], center[y], center[z]);

		/////////////////////////////////////////////////////////////////////////////////////////////
		//	Get narrow band: Place particles on interface
		size_t bc[dims] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
		typedef aggregate<double, Point<dims, double>, double> props_nb;
		typedef vector_dist<dims, double, props_nb> vd_type;
		Ghost<dims, double> ghost_vd(0);
		vd_type vd_narrow_band(0, box, bc, ghost_vd);
		vd_narrow_band.setPropNames({"SDF_exact", "Gradient of Phi", "Gradient magnitude of Phi"});
		size_t narrow_band_width = 8;
		NarrowBand<grid_in_type> narrowBand(g_dist, narrow_band_width); // Instantiation of NarrowBand class
		narrowBand.get_narrow_band<SDF_exact_grid, SDF_vd, Gradient_vd, magnOfGrad_vd>(g_dist, vd_narrow_band);

		BOOST_CHECK(vd_narrow_band.size_local() == 6568);
	}
BOOST_AUTO_TEST_SUITE_END()
