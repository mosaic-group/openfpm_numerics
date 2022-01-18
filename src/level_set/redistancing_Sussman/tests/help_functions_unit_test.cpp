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
		
		// Each processor assigns smaller_value to the first grid node in its domain
		double smaller_value = 0.1;
		auto dom = g_dist.getDomainIterator();
		auto key_0 = dom.get();
		g_dist.template get<Field>(key_0) = smaller_value;
		++dom;
		
		// Afterwards we loop over the other grid nodes and assign them another bigger number
		double bigger_value = 0.5;
		while(dom.isNext())
		{
			auto key = dom.get();
			g_dist.template get<Field>(key) = bigger_value;
			++dom;
		}
		
		// Now we check if get_min_value returns smaller_value
		auto min_value = get_min_val<Field>(g_dist);
		BOOST_CHECK_MESSAGE(min_value == smaller_value, "Checking if smallest value stored in grid is returned.");
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
		
		// Each processor assigns smaller_value to the first grid node in its domain
		double smaller_value = 0.1;
		auto dom = g_dist.getDomainIterator();
		auto key_0 = dom.get();
		g_dist.template get<Field>(key_0) = smaller_value;
		++dom;
		
		// Afterwards we loop over the other grid nodes and assign them another bigger number
		double bigger_value = 0.5;
		while(dom.isNext())
		{
			auto key = dom.get();
			g_dist.template get<Field>(key) = bigger_value;
			++dom;
		}
		
		// Now we check if get_max_value returns bigger_value
		auto max_value = get_max_val<Field>(g_dist);
		BOOST_CHECK_MESSAGE(max_value == bigger_value, "Checking if smallest value stored in grid is returned.");
	}
BOOST_AUTO_TEST_SUITE_END()
