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
	const double EPSILON = 1e-14;
	const size_t F_GAUSSIAN   = 0;
	const size_t VELOCITY         = 1;
	const size_t dF_GAUSSIAN  = 2;
	const size_t dF_UPWIND    = 3;
	const size_t ERROR        = 4;
	
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
		double sigma = 0.3 * (box_upper - box_lower);
		int count = 0;
		size_t N = 32;
//		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"F_GAUSSIAN", "VELOCITY", "dF_GAUSSIAN", "dF_UPWIND", "ERROR"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<F_GAUSSIAN>(key) = gaussian(p, mu, sigma);
				// Initialize the sign of the gaussian function
				g_dist.getProp<VELOCITY>(key) = sgn(g_dist.getProp<F_GAUSSIAN>(key));
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
					g_dist.getProp<dF_GAUSSIAN>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<F_GAUSSIAN>(key);
				}
				++dom;
			}
			
			// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
			get_upwind_gradient<F_GAUSSIAN, VELOCITY, dF_UPWIND>(g_dist, convergence_order, false);
			
			// Get the error between analytical and numerical solution
			get_relative_error<dF_UPWIND, dF_GAUSSIAN, ERROR>(g_dist);
			
			LNorms<double> lNorms;
			lNorms.get_l_norms_grid<ERROR>(g_dist);
			
			if (N==32) BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.38954184748195 + EPSILON, "Checking L2-norm upwind gradient "
																				"order " + std::to_string
																				(convergence_order));
			
			
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
		double sigma = 0.3 * (box_upper - box_lower);
		int count = 0;
		size_t N = 32;
//		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"F_GAUSSIAN", "VELOCITY", "dF_GAUSSIAN", "dF_UPWIND", "ERROR"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<F_GAUSSIAN>(key) = gaussian(p, mu, sigma);
				// Initialize the sign of the gaussian function
				g_dist.getProp<VELOCITY>(key) = sgn(g_dist.getProp<F_GAUSSIAN>(key));
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
					g_dist.getProp<dF_GAUSSIAN>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<F_GAUSSIAN>(key);
				}
				++dom;
			}
			
			// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
			get_upwind_gradient<F_GAUSSIAN, VELOCITY, dF_UPWIND>(g_dist, convergence_order, false);
			
			// Get the error between analytical and numerical solution
			get_relative_error<dF_UPWIND, dF_GAUSSIAN, ERROR>(g_dist);
			
			LNorms<double> lNorms;
			lNorms.get_l_norms_grid<ERROR>(g_dist);
			
			if (N==32) BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.08667855716144 + EPSILON, "Checking L2-norm upwind gradient "
																				"order "
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
		double sigma = 0.3 * (box_upper - box_lower);
		int count = 0;
		size_t N = 32;
//		for (size_t N = 32; N <= 128; N *= 2, ++count)
		{
			const size_t sz[grid_dim] = {N};
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"F_GAUSSIAN", "VELOCITY", "dF_GAUSSIAN", "dF_UPWIND", "ERROR"});
			
			auto gdom = g_dist.getDomainGhostIterator();
			while (gdom.isNext())
			{
				auto key = gdom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				// Initialize grid and ghost with gaussian function
				g_dist.getProp<F_GAUSSIAN>(key) = gaussian(p, mu, sigma);
				// Initialize the sign of the gaussian function
				g_dist.getProp<VELOCITY>(key) = sgn(g_dist.getProp<F_GAUSSIAN>(key));
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
					g_dist.getProp<dF_GAUSSIAN>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<F_GAUSSIAN>(key);
				}
				++dom;
			}
			
			// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
			get_upwind_gradient<F_GAUSSIAN, VELOCITY, dF_UPWIND>(g_dist, convergence_order, false);
			
			// Get the error between analytical and numerical solution
			get_relative_error<dF_UPWIND, dF_GAUSSIAN, ERROR>(g_dist);
			
			LNorms<double> lNorms;
			lNorms.get_l_norms_grid<ERROR>(g_dist);

			if (N==32) BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.03215172234342 + EPSILON, "Checking L2-norm upwind gradient "
																				"order "
						+ std::to_string(convergence_order));
		}
	}
	BOOST_AUTO_TEST_CASE(Upwind_gradient_3D_test)
	{
		{
			const size_t grid_dim  = 3;
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(3);
			typedef aggregate<double, int, Point<grid_dim, double>, Point<grid_dim, double>, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			
			double mu = 0.5 * (box_upper - abs(box_lower));
			double sigma = 0.3 * (box_upper - box_lower);
			
			size_t N = 32;
//			for(size_t N = 32; N <=256; N *=2)
			{
				const size_t sz[grid_dim] = {N, N, N};
				grid_in_type g_dist(sz, box, ghost);
				g_dist.setPropNames({"F_GAUSSIAN", "VELOCITY", "dF_GAUSSIAN", "dF_UPWIND", "ERROR"});
				
				auto gdom = g_dist.getDomainGhostIterator();
				while (gdom.isNext())
				{
					auto key = gdom.get();
					Point<grid_dim, double> p = g_dist.getPos(key);
					// Initialize grid and ghost with gaussian function
					g_dist.getProp<F_GAUSSIAN>(key) = gaussian(p, mu, sigma);
					// Initialize the sign of the gaussian function
					g_dist.getProp<VELOCITY>(key) = sgn(g_dist.getProp<F_GAUSSIAN>(key));
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
						g_dist.getProp<dF_GAUSSIAN>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<F_GAUSSIAN>(key);
					}
					++dom;
				}
				
				for (int convergence_order = 1; convergence_order <=5; convergence_order +=2 )
				{
					// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
					get_upwind_gradient<F_GAUSSIAN, VELOCITY, dF_UPWIND>(g_dist, convergence_order, false);
					// Get the error between analytical and numerical solution
					get_relative_error<dF_UPWIND, dF_GAUSSIAN, ERROR>(g_dist);
					
					LNorms<double> lNorms;
					lNorms.get_l_norms_grid<ERROR>(g_dist);

					if(N==32)
					{
						switch(convergence_order)
						{
							case 1:
								BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.38954184748194 + EPSILON, "Checking L2-norm upwind gradient order " + std::to_string(convergence_order));
								break;
							case 3:
								BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.08667855716144 + EPSILON, "Checking L2-norm upwind gradient order " + std::to_string(convergence_order));
								break;
							case 5:
								BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.02285996528578 + EPSILON, "Checking L2-norm upwind gradient order " + std::to_string(convergence_order));

								break;
							default:
								std::cout << "Checking only implemented for convergence order 1, 3 and 5." << std::endl;
						}
					}
				}
				
			}
		}
	}
	
	BOOST_AUTO_TEST_CASE(Upwind_gradient_3D_vector_velocity_test)
	{
		{
			const size_t grid_dim  = 3;
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower, box_lower, box_lower}, {box_upper, box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(3);
			typedef aggregate<double, double[grid_dim], Point<grid_dim, double>, Point<grid_dim, double>, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			
			double mu = 0.5 * (box_upper - abs(box_lower));
			double sigma = 0.3 * (box_upper - box_lower);
			
			size_t N = 32;
//			for(size_t N = 32; N <=256; N *=2)
			{
				const size_t sz[grid_dim] = {N, N, N};
				grid_in_type g_dist(sz, box, ghost);
				g_dist.setPropNames({"F_GAUSSIAN", "Velocity", "dF_GAUSSIAN", "dF_UPWIND", "ERROR"});
				
				auto gdom = g_dist.getDomainGhostIterator();
				while (gdom.isNext())
				{
					auto key = gdom.get();
					Point<grid_dim, double> p = g_dist.getPos(key);
					// Initialize grid and ghost with gaussian function
					g_dist.getProp<F_GAUSSIAN>(key) = gaussian(p, mu, sigma);
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
						g_dist.getProp<dF_GAUSSIAN>(key)[d] = hermite_polynomial(p.get(d), sigma, 1) * g_dist.getProp<F_GAUSSIAN>(key);
						// Initialize the velocity vector field
						g_dist.getProp<VELOCITY>(key)[d] = g_dist.getProp<dF_GAUSSIAN>(key)[d];
					}
					++dom;
				}
				
				for (int convergence_order = 1; convergence_order <=5; convergence_order +=2 )
				{
					// Upwind gradient with specific convergence order, not becoming one-sided at the boundary
					get_upwind_gradient<F_GAUSSIAN, VELOCITY, dF_UPWIND>(g_dist, convergence_order, false);
					// Get the error between analytical and numerical solution
					get_relative_error<dF_UPWIND, dF_GAUSSIAN, ERROR>(g_dist);
					
					LNorms<double> lNorms;
					lNorms.get_l_norms_grid<ERROR>(g_dist);
					
					if(N==32)
					{
						switch(convergence_order)
						{
							case 1:
								BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.38954184748194 + EPSILON, "Checking L2-norm upwind"
																							 " gradient with vector "
																							 "velocity order " +
																							 std::to_string(convergence_order));
								break;
							case 3:
								BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.08667855716144 + EPSILON, "Checking L2-norm upwind"
																							 " gradient with "
																							 "vector velocity order " +
																							 std::to_string(convergence_order));
								break;
							case 5:
								BOOST_CHECK_MESSAGE(lNorms.l2 <= 0.02285996528578 + EPSILON, "Checking L2-norm upwind"
																							 " gradient with "
																							 "vector velocity order " +
																							 std::to_string(convergence_order));
								
								break;
							default:
								std::cout << "Checking only implemented for convergence order 1, 3 and 5." << std::endl;
						}
					}
				}
				
			}
		}
	}
BOOST_AUTO_TEST_SUITE_END()
