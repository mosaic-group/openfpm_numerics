/*
 * vector_dist_operators_unit_tests.hpp
 *
 *  Created on: Jun 11, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_UNIT_TESTS_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "tests/vector_dist_operators_tests_util.hpp"


////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( vector_dist_operators_test )

BOOST_AUTO_TEST_CASE( vector_dist_operators_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	// vector type
	typedef vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> vtype;

	vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>>> vd(100,box,bc,ghost);

	check_all_expressions<comp_host>::check(vd);
}


void reset(vector_dist<3,float,aggregate<float,float,float,float,float,float,float,float>> & v1)
{
	auto it = v1.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		v1.template getProp<0>(p) = -1.0;
		v1.template getProp<1>(p) = -1.0;
		v1.template getProp<2>(p) = -1.0;
		v1.template getProp<3>(p) = -1.0;
		v1.template getProp<4>(p) = -1.0;
		v1.template getProp<5>(p) = -1.0;
		v1.template getProp<6>(p) = -1.0;
		v1.template getProp<7>(p) = -1.0;

		++it;
	}
}

template<unsigned int i> void loop_check(vector_dist<3,float,aggregate<float,float,float,float,float,float,float,float>> & v1, float check)
{
	auto it = v1.getDomainIterator();
	bool ret = true;

	while (it.isNext())
	{
		auto p = it.get();

		ret &= v1.template getProp<i>(p) == check;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);
}

void check(vector_dist<3,float,aggregate<float,float,float,float,float,float,float,float>> & v1, size_t i)
{
	float check = -1.0;

	if (0 < i) check = 0.0;
	else check = -1.0;
	loop_check<0>(v1,check);

	if (1 < i) check = 1.0;
	else check = -1.0;
	loop_check<1>(v1,check);

	if (2 < i) check = 2.0;
	else check = -1.0;
	loop_check<2>(v1,check);

	if (3 < i) check = 3.0;
	else check = -1.0;
	loop_check<3>(v1,check);

	if (4 < i) check = 4.0;
	else check = -1.0;
	loop_check<4>(v1,check);

	if (5 < i) check = 5.0;
	else check = -1.0;
	loop_check<5>(v1,check);

	if (6 < i) check = 6.0;
	else check = -1.0;
	loop_check<6>(v1,check);

	if (7 < i) check = 7.0;
	else check = -1.0;
	loop_check<7>(v1,check);
}

BOOST_AUTO_TEST_CASE( vector_dist_operators_assign_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	vector_dist<3,float,aggregate<float,float,float,float,float,float,float,float>> vd(100,box,bc,ghost);

	auto v1 = getV<0>(vd);
	auto v2 = getV<1>(vd);
	auto v3 = getV<2>(vd);
	auto v4 = getV<3>(vd);
	auto v5 = getV<4>(vd);
	auto v6 = getV<5>(vd);
	auto v7 = getV<6>(vd);
	auto v8 = getV<7>(vd);

	auto v_pos = getV<POS_PROP>(vd);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0);

	check(vd,2);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0,
		   v3,2.0);

	check(vd,3);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0,
		   v3,2.0,
		   v4,3.0);

	check(vd,4);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0,
		   v3,2.0,
		   v4,3.0,
		   v5,4.0);

	check(vd,5);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0,
		   v3,2.0,
		   v4,3.0,
		   v5,4.0,
		   v6,5.0);

	check(vd,6);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0,
		   v3,2.0,
		   v4,3.0,
		   v5,4.0,
		   v6,5.0,
		   v7,6.0);

	check(vd,7);

	reset(vd);

	assign(v1,0.0,
		   v2,1.0,
		   v3,2.0,
		   v4,3.0,
		   v5,4.0,
		   v6,5.0,
		   v7,6.0,
		   v8,7.0);

	check(vd,8);

	// Self assign;
	assign(v_pos,v_pos,
		   v_pos,v_pos);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_UNIT_TESTS_HPP_ */
