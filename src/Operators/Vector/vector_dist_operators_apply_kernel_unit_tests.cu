/*
 * vector_dist_operators_apply_kernel_unit_tests.hpp
 *
 *  Created on: Jun 19, 2016
 *      Author: i-bird
 */

#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Vector/vector_dist.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Operators/Vector/tests/vector_dist_operators_tests_util.hpp"

BOOST_AUTO_TEST_CASE( vector_dist_operators_apply_host_gpu_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	vector_dist_gpu<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>,float>> vd(512,box,bc,ghost);

	check_all_apply_ker<comp_host>::check(vd);
}


BOOST_AUTO_TEST_CASE( vector_dist_operators_apply_kernel_gpu_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	vector_dist_gpu<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>,float>> vd(512,box,bc,ghost);

	check_all_apply_ker<comp_dev>::check(vd);
}

