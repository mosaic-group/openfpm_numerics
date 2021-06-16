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

BOOST_AUTO_TEST_SUITE(FDOrder1TestSuite)
	const size_t Field              = 0;
	const size_t AnalyticalGradient = 1;
	const size_t NumericalGradient  = 2;
	const size_t Error              = 3;
	
	// 32, 0.573022, 1.33305
	// 64, 0.288525, 1
	//128, 0.171455, 1
	double l2_norms [] = {0.573022, 0.288525, 0.171455};
	
	const double EPSILON = std::numeric_limits<double>::epsilon();
	
	BOOST_AUTO_TEST_CASE(Forward_difference_1D_test)
	{
		const size_t grid_dim  = 1;
		
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double, Point<grid_dim, double>, Point<grid_dim, double>, double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_in_type;
		
		double mu = 0.5 * (box_upper - abs(box_lower));
		double sigma = 0.1 * (box_upper - box_lower);
		
		int count = 0;
		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"Field", "AnalyticalGradient", "NumericalGradient", "Error"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<Field>(key) = gaussian(p, mu, sigma);
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
					g_dist.getProp<AnalyticalGradient>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<Field>(key);
				    // 1st order Finite Difference gradient
					g_dist.getProp<NumericalGradient>(key)[d]  = FD_forward<Field>(g_dist, key, d);
				}
				++dom;
			}
			
			// Get the error between analytical and numerical solution
//			get_absolute_error<NumericalGradient, AnalyticalGradient, Error>(g_dist);
			get_relative_error<NumericalGradient, AnalyticalGradient, Error>(g_dist);
			
			L_norms lNorms;
			lNorms = get_l_norms_grid<Error>(g_dist);
			BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms[count] + 0.00001 + EPSILON, "Checking L2-norm ENO");
//			write_lnorms_to_file(N, lNorms, "l_norms_FDfwd", "./");
			std::cout << N << ", " << lNorms.l2 << ", " << lNorms.linf << std::endl;
//			if (N==128) g_dist.write("grid_gaussian_FDfwd_N" + std::to_string(N), FORMAT_BINARY);
		}
	}
	BOOST_AUTO_TEST_CASE(Backward_difference_1D_test)
	{
		const size_t grid_dim  = 1;
		
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double, Point<grid_dim, double>, Point<grid_dim, double>, double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_in_type;
		
		double mu = 0.5 * (box_upper - abs(box_lower));
		double sigma = 0.1 * (box_upper - box_lower);
		
		int count = 0;
		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"Field", "AnalyticalGradient", "NumericalGradient", "Error"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<Field>(key) = gaussian(p, mu, sigma);
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
					g_dist.getProp<AnalyticalGradient>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<Field>(key);
					// 1st order Finite Difference gradient
					g_dist.getProp<NumericalGradient>(key)[d]  = FD_backward<Field>(g_dist, key, d);
				}
				++dom;
			}
			
			// Get the error between analytical and numerical solution
//			get_absolute_error<NumericalGradient, AnalyticalGradient, Error>(g_dist);
			get_relative_error<NumericalGradient, AnalyticalGradient, Error>(g_dist);
			
			L_norms lNorms;
			lNorms = get_l_norms_grid<Error>(g_dist);
			BOOST_CHECK_MESSAGE(lNorms.l2   < l2_norms[count] + 0.00001 + EPSILON, "Checking L2-norm ENO");
//			write_lnorms_to_file(N, lNorms, "l_norms_FDbwd", "./");
			std::cout << N << ", " << lNorms.l2 << ", " << lNorms.linf << std::endl;
//			if (N==128) g_dist.write("grid_gaussian_FDbwd_N" + std::to_string(N), FORMAT_BINARY);
		}
	}
BOOST_AUTO_TEST_SUITE_END()
