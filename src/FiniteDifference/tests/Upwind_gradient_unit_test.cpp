//
// Created by jstark on 16.06.21.
//

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// Include header files for testing
#include "level_set/redistancing_Sussman/tests/l_norms/LNorms.hpp"
#include "Gaussian.hpp"
#include "FiniteDifference/Upwind_gradient.hpp"
#include "level_set/redistancing_Sussman/HelpFunctions.hpp"

BOOST_AUTO_TEST_SUITE(UpwindGradientTestSuite)
	const size_t f_gaussian   = 0;
	const size_t Sign         = 1;
	const size_t df_gaussian  = 2;
	const size_t df_upwind    = 3;
	const size_t Error        = 4;
	
	BOOST_AUTO_TEST_CASE(Upwind_gradient_order1_1D_test)
	{
		const int convergence_order = 1;
		const size_t grid_dim  = 1;
		
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double, int, Point<grid_dim, double>, Point<grid_dim, double>, double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_in_type;
		
		double mu = 0.5 * (box_upper - abs(box_lower));
		double sigma = 0.1 * (box_upper - box_lower);
		int count = 0;
		size_t N = 32;
//		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"f_gaussian", "Sign", "df_gaussian", "df_upwind", "Error"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<f_gaussian>(key) = gaussian(p, mu, sigma);
				// Initialize the sign of the gaussian function
				g_dist.getProp<Sign>(key) = sgn(g_dist.getProp<f_gaussian>(key));
				++gdom;
			}
			
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				
				for (int d = 0; d < grid_dim; d++)
				{
					// Analytical gradient
					g_dist.getProp<df_gaussian>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<f_gaussian>(key);
				}
				++dom;
			}
			
			// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
			get_upwind_gradient<f_gaussian, Sign, df_upwind>(g_dist, convergence_order, false);
			
			// Get the error between analytical and numerical solution
			get_relative_error<df_upwind, df_gaussian, Error>(g_dist);
			
			L_norms lNorms;
			lNorms = get_l_norms_grid<Error>(g_dist);
			// std::cout << N << ", " << lNorms.l2 << ", " << lNorms.linf << std::endl;
			
			// write_lnorms_to_file(N, lNorms, "l_norms_upwindGrad_convOrder_" + std::to_string(convergence_order), ""
//																							"./");
			// if(N==128) g_dist.write("grid_gaussian_upwindGradient_N" + std::to_string(N), FORMAT_BINARY);
			
			// 32,0.389542,0.899459
			if (N==32) BOOST_CHECK_MESSAGE(lNorms.l2 < 0.389542 + 0.000001, "Checking L2-norm upwind gradient order "
			+ std::to_string(convergence_order));
			
			
		}
	}
	BOOST_AUTO_TEST_CASE(Upwind_gradient_order3_1D_test)
	{
		const int convergence_order = 3;
		const size_t grid_dim  = 1;
		
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double, int, Point<grid_dim, double>, Point<grid_dim, double>, double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_in_type;
		
		double mu = 0.5 * (box_upper - abs(box_lower));
		double sigma = 0.1 * (box_upper - box_lower);
		int count = 0;
		size_t N = 32;
//		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"f_gaussian", "Sign", "df_gaussian", "df_upwind", "Error"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<f_gaussian>(key) = gaussian(p, mu, sigma);
				// Initialize the sign of the gaussian function
				g_dist.getProp<Sign>(key) = sgn(g_dist.getProp<f_gaussian>(key));
				++gdom;
			}
			
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				
				for (int d = 0; d < grid_dim; d++)
				{
					// Analytical gradient
					g_dist.getProp<df_gaussian>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<f_gaussian>(key);
				}
				++dom;
			}
			
			// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
			get_upwind_gradient<f_gaussian, Sign, df_upwind>(g_dist, convergence_order, false);
			
			// Get the error between analytical and numerical solution
			get_relative_error<df_upwind, df_gaussian, Error>(g_dist);
			
			L_norms lNorms;
			lNorms = get_l_norms_grid<Error>(g_dist);
			std::cout << N << ", " << lNorms.l2 << ", " << lNorms.linf << std::endl;
			
//			write_lnorms_to_file(N, lNorms, "l_norms_upwindGrad_convOrder_" + std::to_string(convergence_order), ""
//			                                                                                                     "./");
//			if(N==128) g_dist.write("grid_gaussian_upwindGradient_N" + std::to_string(N), FORMAT_BINARY);
			
			if (N==32) BOOST_CHECK_MESSAGE(lNorms.l2 < 0.086679 + 0.000001, "Checking L2-norm upwind gradient order "
						+ std::to_string(convergence_order));
		}
	}
	BOOST_AUTO_TEST_CASE(Upwind_gradient_order5_1D_test)
	{
		const int convergence_order = 5;
		const size_t grid_dim  = 1;
		
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double, int, Point<grid_dim, double>, Point<grid_dim, double>, double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_in_type;
		
		double mu = 0.5 * (box_upper - abs(box_lower));
		double sigma = 0.1 * (box_upper - box_lower);
		int count = 0;
		size_t N = 32;
//		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"f_gaussian", "Sign", "df_gaussian", "df_upwind", "Error"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<f_gaussian>(key) = gaussian(p, mu, sigma);
				// Initialize the sign of the gaussian function
				g_dist.getProp<Sign>(key) = sgn(g_dist.getProp<f_gaussian>(key));
				++gdom;
			}
			
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				
				for (int d = 0; d < grid_dim; d++)
				{
					// Analytical gradient
					g_dist.getProp<df_gaussian>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<f_gaussian>(key);
				}
				++dom;
			}
			
			// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
			get_upwind_gradient<f_gaussian, Sign, df_upwind>(g_dist, convergence_order, false);
			
			// Get the error between analytical and numerical solution
			get_relative_error<df_upwind, df_gaussian, Error>(g_dist);
			
			L_norms lNorms;
			lNorms = get_l_norms_grid<Error>(g_dist);
			std::cout << N << ", " << lNorms.l2 << ", " << lNorms.linf << std::endl;
			
//			write_lnorms_to_file(N, lNorms, "l_norms_upwindGrad_convOrder_" + std::to_string(convergence_order), ""
//			                                                                                                     "./");
//			if(N==128) g_dist.write("grid_gaussian_upwindGradient_N" + std::to_string(N), FORMAT_BINARY);
			
			if (N==32) BOOST_CHECK_MESSAGE(lNorms.l2 < 0.032152 + 0.000001, "Checking L2-norm upwind gradient order "
						+ std::to_string(convergence_order));
		}
	}
BOOST_AUTO_TEST_SUITE_END()
