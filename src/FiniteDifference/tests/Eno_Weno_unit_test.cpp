//
// Created by jstark on 14.06.21.
//

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// Include header files for testing
#include "level_set/redistancing_Sussman/tests/l_norms/LNorms.hpp"
#include "Gaussian.hpp"
#include "FiniteDifference/Eno_Weno.hpp"

BOOST_AUTO_TEST_SUITE(EnoWenoTestSuite)
	const size_t f_gaussian   = 0;
	const size_t df_gaussian  = 1;
	const size_t ENO_plus     = 2;
	const size_t ENO_minus    = 3;
	const size_t Error_plus   = 4;
	const size_t Error_minus  = 5;
	const double EPSILON = std::numeric_limits<double>::epsilon();
	
	// 32,0.079885,0.186274
	// 64,0.014594,0.045086
	//128,0.002300,0.007982
	double l2_norms_ENO []    = {0.079885, 0.014594, 0.002300};
	double linf_norms_ENO []  = {0.186274, 0.045086, 0.007982};
	
	// 32,0.033198,0.097141
	// 64,0.000935,0.002459
	//128,0.000052,0.000143
	double l2_norms_WENO []   = {0.033198, 0.000935, 0.000052};
	double linf_norms_WENO [] = {0.097141, 0.002459, 0.000143};
	
	BOOST_AUTO_TEST_CASE(Eno_1D_1stDerivative_test)
	{
		int count = 0;
		for (size_t N = 32; N <= 128; N *= 2, ++count)
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
			get_relative_error<ENO_plus, df_gaussian, Error_plus>(g_dist);
			get_relative_error<ENO_minus, df_gaussian, Error_minus>(g_dist);
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_plus>(g_dist);
				BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms_ENO[count] + 0.000001 + EPSILON, "Checking L2-norm ENO");
				BOOST_CHECK_MESSAGE(lNorms.linf < linf_norms_ENO[count] + 0.000001 + EPSILON, "Checking Linf-norm "
				                                                                             "ENO");
//				write_lnorms_to_file(N, lNorms, "l_norms_ENO_plus", "./");
			}
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_minus>(g_dist);
				BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms_ENO[count] + 0.000001 + EPSILON, "Checking L2-norm ENO");
				BOOST_CHECK_MESSAGE(lNorms.linf < linf_norms_ENO[count] + 0.000001 + EPSILON, "Checking Linf-norm ENO");
//				write_lnorms_to_file(N, lNorms, "l_norms_ENO_minus", "./");
			}
//			g_dist.write("grid_gaussian_ENO_1D_N" + std::to_string(N), FORMAT_BINARY);
		}
	}
	const size_t WENO_plus     = 2;
	const size_t WENO_minus    = 3;
	BOOST_AUTO_TEST_CASE(Weno_1D_1stDerivative_test)
	{
		int count = 0;
		for (size_t N = 32; N <= 128; N *= 2, ++count)
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
			
			g_dist.setPropNames({"f_gaussian", "df_gaussian", "WENO_plus", "WENO_minus", "Error_plus", "Error_minus"});
			
			double mu = 0.5 * (box_upper - abs(box_lower));
			double sigma = 0.1 * (box_upper - box_lower);
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<f_gaussian>(key)  = gaussian(p, mu, sigma);
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
					// WENO_plus
					g_dist.getProp<WENO_plus>(key)[d]   = WENO_5_Plus<f_gaussian>(g_dist, key, d);
					// WENO_minus
					g_dist.getProp<WENO_minus>(key)[d]  = WENO_5_Minus<f_gaussian>(g_dist, key, d);
				}
				
				++dom;
			}
			// Get the error between analytical and numerical solution
			get_relative_error<WENO_plus, df_gaussian, Error_plus>(g_dist);
			get_relative_error<WENO_minus, df_gaussian, Error_minus>(g_dist);
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_plus>(g_dist);
				BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms_WENO[count] + 0.000001 + EPSILON, "Checking L2-norm WENO");
				BOOST_CHECK_MESSAGE(lNorms.linf < linf_norms_WENO[count] + 0.000001 + EPSILON, "Checking Linf-norm "
																							  "WENO");
//				write_lnorms_to_file(N, lNorms, "l_norms_WENO_plus", "./");
			}
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_minus>(g_dist);
				BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms_WENO[count] + 0.000001 + EPSILON, "Checking L2-norm WENO");
				BOOST_CHECK_MESSAGE(lNorms.linf < linf_norms_WENO[count] + 0.000001 + EPSILON, "Checking Linf-norm "
				                                                                              "WENO");
//				write_lnorms_to_file(N, lNorms, "l_norms_WENO_minus", "./");
			}
//			if(N==128) g_dist.write("grid_gaussian_WENO_1D_N" + std::to_string(N), FORMAT_BINARY);
		}
	}
	
	BOOST_AUTO_TEST_CASE(Eno_2D_1stDerivative_test)
	{
		int count = 0;
		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t grid_dim = 2;
			const size_t sz[grid_dim] = {N, N};
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower, box_lower}, {box_upper, box_upper});
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
				g_dist.getProp<f_gaussian>(key)  = gaussian(p, mu, sigma);
				++gdom;
			}
			
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				for (int d = 0; d < grid_dim; d++)
				{
					g_dist.getProp<df_gaussian>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<f_gaussian>(key);
					g_dist.getProp<ENO_plus>(key)[d]    = ENO_3_Plus<f_gaussian>(g_dist, key, d);
					g_dist.getProp<ENO_minus>(key)[d]   = ENO_3_Minus<f_gaussian>(g_dist, key, d);
				}
				
				++dom;
			}
			
			// Get the error between analytical and numerical solution
			get_relative_error<ENO_plus, df_gaussian, Error_plus>(g_dist);
			get_relative_error<ENO_minus, df_gaussian, Error_minus>(g_dist);
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_plus>(g_dist);
				BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms_ENO[count] + 0.000001 + EPSILON, "Checking L2-norm ENO");
				BOOST_CHECK_MESSAGE(lNorms.linf < linf_norms_ENO[count] + 0.000001 + EPSILON, "Checking Linf-norm "
				                                                                              "ENO");
//				write_lnorms_to_file(N, lNorms, "l_norms_2D_ENO_plus", "./");
			}
			
			{
				L_norms lNorms;
				lNorms = get_l_norms_grid<Error_minus>(g_dist);
				BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms_ENO[count] + 0.000001 + EPSILON, "Checking L2-norm ENO");
				BOOST_CHECK_MESSAGE(lNorms.linf < linf_norms_ENO[count] + 0.000001 + EPSILON, "Checking Linf-norm "
				                                                                             "ENO");
//				write_lnorms_to_file(N, lNorms, "l_norms_2D_ENO_minus", "./");
			}
//			g_dist.write("grid_gaussian_2D_N" + std::to_string(N), FORMAT_BINARY);
		}
	}
BOOST_AUTO_TEST_SUITE_END()