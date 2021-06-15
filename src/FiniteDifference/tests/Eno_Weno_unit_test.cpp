//
// Created by jstark on 14.06.21.
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

BOOST_AUTO_TEST_SUITE(EnoWenoTestSuite)
	const size_t f_gaussian   = 0;
	const size_t df_gaussian  = 1;
	const size_t ddf_gaussian = 2;
	const size_t Error1       = 3;
	const size_t Error2       = 4;
	
	BOOST_AUTO_TEST_CASE(Eno_1D_1stDerivative_test)
	{
		for (size_t N = 64; N <= 64; N *= 2)
		{
			const size_t grid_dim = 1;
			const size_t sz[grid_dim] = {N};
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower}, {box_upper});
			Ghost<grid_dim, long int> ghost(0);
			typedef aggregate<double, Point<grid_dim, double>, Point<grid_dim, double>, double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			
			g_dist.setPropNames({"f_gaussian", "df_gaussian", "ddf_gaussian", "Error gradient", "Error Laplace"});

			double mu = 0.5 * (box_upper - abs(box_lower));
			double sigma = 0.1 * (box_upper - box_lower);
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				g_dist.getProp<f_gaussian>(key) = gaussian(p, mu, sigma);
				for (int d = 0; d < grid_dim; d++)
				{
					g_dist.getProp<df_gaussian>(key)[d] = gaussian_dx(p.get(d), mu, sigma);
//					g_dist.getProp<df_gaussian>(key)[d] = 1;

				}

				++dom;
			}
			
			g_dist.write("grid_gaussian_1D_N" + std::to_string(N), FORMAT_BINARY);
			
		}
	}
	BOOST_AUTO_TEST_CASE(Eno_2D_1stDerivative_test)
	{
		for(size_t N=32; N<=128; N*=2)
		{
			const size_t grid_dim = 2;
			const size_t sz[grid_dim] = {N, N};
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower, box_lower}, {box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(0);
			typedef aggregate<double, Point<grid_dim, double>, Point<grid_dim, double>, double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"f_gaussian", "df_gaussian", "ddf_gaussian", "Error gradient", "Error Laplace"});
			
			double mu = 0.5 * (box_upper - abs(box_lower));
			double sigma = 0.1 * (box_upper - box_lower);
			auto dom = g_dist.getDomainIterator();
			while (dom.isNext())
			{
				auto key = dom.get();
				Point<grid_dim, double> p = g_dist.getPos(key);
				g_dist.getProp<f_gaussian>(key) = gaussian(p, mu, sigma);
				for (int d = 0; d < grid_dim; d++)
				{
					g_dist.getProp<df_gaussian>(key)[d] = gaussian_dx(p.get(d), mu, sigma);
				}
				
				++dom;
			}
			
			g_dist.write("grid_gaussian_2D_N" + std::to_string(N), FORMAT_BINARY);
		}
	}
BOOST_AUTO_TEST_SUITE_END()