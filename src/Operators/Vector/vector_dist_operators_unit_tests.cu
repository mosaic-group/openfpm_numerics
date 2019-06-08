#include "config.h"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "Vector/vector_dist.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Operators/Vector/tests/vector_dist_operators_tests_util.hpp"

BOOST_AUTO_TEST_SUITE( vector_dist_operators_test_gpu )

BOOST_AUTO_TEST_CASE(vector_dist_operators_list_ker_test)
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.2);

	vector_dist_gpu<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> vd(127*create_vcluster().size(),box,bc,ghost);

	typedef decltype(vd.toKernel()) vector_dist_kernel;

	vector_dist_ker_list<vector_dist_kernel> & vdkl = vd.private_get_vector_dist_ker_list();

	auto vdk = vd.toKernel();
	BOOST_REQUIRE_EQUAL(vdkl.n_entry(),0);

	{
		auto vA = getV<A>(vd,vdk);

		BOOST_REQUIRE_EQUAL(vdkl.n_entry(),1);

		{
			auto vB = getV<B>(vd,vdk);
			auto vC = getV<C>(vd,vdk);

			auto vVA = getV<VA>(vd,vdk);
			auto vVB = getV<VB>(vd,vdk);
			auto vVC = getV<VC>(vd,vdk);

			BOOST_REQUIRE_EQUAL(vdkl.n_entry(),6);

			// fill vd with some value
			fill_values<comp_dev>(vd);

			vd.map(RUN_ON_DEVICE);

			// check that the entry has been updated
			BOOST_REQUIRE_EQUAL(vdkl.check(vd.toKernel()),true);

			vd.template ghost_get<0,1,2,3,4,5>(RUN_ON_DEVICE);

			// check that the entry has been updated
			BOOST_REQUIRE_EQUAL(vdkl.check(vd.toKernel()),true);
		}

		BOOST_REQUIRE_EQUAL(vdkl.n_entry(),1);
	}

	BOOST_REQUIRE_EQUAL(vdkl.n_entry(),0);
}

BOOST_AUTO_TEST_CASE( vector_dist_operators_gpu_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	// vector type
	typedef vector_dist_gpu<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> vtype;

	vector_dist_gpu<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> vd(100,box,bc,ghost);

	check_all_expressions<comp_host>::check(vd);
	check_all_expressions<comp_dev>::check(vd);
}

BOOST_AUTO_TEST_SUITE_END()
