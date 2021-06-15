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
#include "FiniteDifference/Eno_Weno.hpp"

BOOST_AUTO_TEST_SUITE(EnoWenoTestSuite)
	const size_t f_gaussian   = 0;
	const size_t df_gaussian  = 1;
	const size_t ENO_plus     = 2;
	const size_t ENO_minus    = 3;
	const size_t Error_plus   = 4;
	const size_t Error_minus  = 5;
	
	BOOST_AUTO_TEST_CASE(Eno_1stDerivative_1D_test)
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
			/*
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
			*/
//			g_dist.write("grid_gaussian_1D_N" + std::to_string(N));
		
		}
	}
	BOOST_AUTO_TEST_CASE(Eno_1stDerivative_2D_test)
	{
		for(size_t N=32; N<=128; N*=2)
		{
			const size_t grid_dim = 2;
			const size_t sz[grid_dim] = {N, N};
			const double box_lower = -1.0;
			const double box_upper = 1.0;
			Box<grid_dim, double> box({box_lower, box_lower}, {box_upper, box_upper});
			Ghost<grid_dim, long int> ghost(6);
			typedef aggregate<double, Point<grid_dim, double>, Point<grid_dim, double>, Point<grid_dim, double>,
			        double, double> props;
			typedef grid_dist_id<grid_dim, double, props> grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"f_gaussian", "df_gaussian", "ENO_plus", "ENO_minus", "Error plus", "Error minus"});
			
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
//					template<typename GridDist, typename GridKey, unsigned int field>
//					double ENO_3_Plus(GridKey key, GridDist &gd, size_t x)
					g_dist.getProp<ENO_plus>(key)[d] = ENO_3_Plus<f_gaussian>(g_dist, key, d);

				}
				++dom;
			}
			
			
			
			
			
			g_dist.write("grid_gaussian_2D_N" + std::to_string(N));
		}
	}
BOOST_AUTO_TEST_SUITE_END()