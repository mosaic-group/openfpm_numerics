//
// Created by jstark on 03.01.22.
//
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// Include header files for testing
#include "level_set/redistancing_Sussman/HelpFunctions.hpp"
#include "level_set/redistancing_Sussman/HelpFunctionsForGrid.hpp"

BOOST_AUTO_TEST_SUITE(HelpFunctionsTestSuite)
	const size_t Field              = 0;

	BOOST_AUTO_TEST_CASE(get_min_val_test)
	{
		const size_t grid_dim  = 1;
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_type;
		const size_t sz[grid_dim] = {32};
		grid_type g_dist(sz, box, ghost);
		
		double i = 3.3;
		auto dom = g_dist.getDomainIterator();
		while(dom.isNext())
		{
			i -= 0.1;
			auto key = dom.get();
			g_dist.template get<Field>(key) = i;
			++dom;
		}
		auto correct_value = i;
		auto min_value = get_min_val<Field>(g_dist);
//		BOOST_CHECK_MESSAGE(min_value == correct_value, "Checking if minimum value stored in grid is returned.");
	}
	
	BOOST_AUTO_TEST_CASE(get_max_val_test)
	{
		const size_t grid_dim  = 1;
		const double box_lower = -1.0;
		const double box_upper = 1.0;
		Box<grid_dim, double> box({box_lower}, {box_upper});
		Ghost<grid_dim, long int> ghost(3);
		typedef aggregate<double> props;
		typedef grid_dist_id<grid_dim, double, props> grid_type;
		const size_t sz[grid_dim] = {32};
		grid_type g_dist(sz, box, ghost);
		
		double i = 0.0;
		auto dom = g_dist.getDomainIterator();
		while(dom.isNext())
		{
			i += 0.1;
			auto key = dom.get();
			g_dist.template get<Field>(key) = i;
			++dom;
		}
		auto correct_value = i;
		auto max_value = get_max_val<Field>(g_dist);
//		BOOST_CHECK_MESSAGE(max_value == correct_value, "Checking if maximum value stored in grid is returned.");
	}
BOOST_AUTO_TEST_SUITE_END()
