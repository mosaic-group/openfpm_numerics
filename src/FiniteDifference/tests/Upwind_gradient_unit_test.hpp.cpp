//
// Created by jstark on 16.06.21.
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
#include "level_set/redistancing_Sussman/tests/l_norms/LNorms.hpp"
#include "Gaussian.hpp"

BOOST_AUTO_TEST_SUITE(UpwindGradientTestSuite)
	const size_t f_gaussian   = 0;
	const size_t df_gaussian  = 1;
	const size_t ENO_plus     = 2;
	const size_t ENO_minus    = 3;
	const size_t Error_plus   = 4;
	const size_t Error_minus  = 5;
	
	
	BOOST_AUTO_TEST_CASE(Upwind_gradient_1D_test)
	{
		for (size_t N = 32; N <= 256; N *= 2)
		{
			const size_t grid_dim = 1;
			const size_t sz[grid_dim] = {N};
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower}, {box_upper});
			Ghost<grid_dim, long int> ghost(3);
			typedef aggregate<double, Point<grid_dim, double>, Point<grid_dim, double>, Point<grid_dim, double>, double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"f_gaussian", "df_gaussian", "ENO_plus", "ENO_minus", "Error_plus", "Error_minus"});
			
			double mu = 0.5 * (box_upper - abs(box_lower));
			double sigma = 0.1 * (box_upper - box_lower);
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<f_gaussian>(key) = gaussian(p, mu, sigma);
				++gdom;
			}
			
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				
				for (int d = 0; d < grid_dim; d++)
				{
					// Analytical solution
					g_dist.getProp<df_gaussian>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<f_gaussian>(key);
					// ENO_plus
					g_dist.getProp<ENO_plus>(key)[d]  = ENO_3_Plus<f_gaussian>(g_dist, key, d);
					// ENO_minus
					g_dist.getProp<ENO_minus>(key)[d]  = ENO_3_Minus<f_gaussian>(g_dist, key, d);
				}
				
				++dom;
			}
			// Get the error between analytical and numerical solution
			get_relative_error<df_gaussian, ENO_plus, Error_plus>(g_dist);
			get_relative_error<df_gaussian, ENO_minus, Error_minus>(g_dist);
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_plus>(g_dist);
				write_lnorms_to_file(N, lNorms, "l_norms_ENO_plus", "./");
			}
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_minus>(g_dist);
				write_lnorms_to_file(N, lNorms, "l_norms_ENO_minus", "./");
			}
			
			g_dist.write("grid_gaussian_ENO_1D_N" + std::to_string(N), FORMAT_BINARY);
			
		}
	}
BOOST_AUTO_TEST_SUITE_END()
