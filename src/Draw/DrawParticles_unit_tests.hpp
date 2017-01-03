/*
 * DrawParticles_unit_tests.hpp
 *
 *  Created on: Jan 3, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_UNIT_TESTS_HPP_

#include "PointIterator.hpp"

BOOST_AUTO_TEST_SUITE( draw_particles )

BOOST_AUTO_TEST_CASE(point_iterator)
{
	size_t sz[] = {30,17,28};

	Box<3,double> domain({-1.2,0.5,-0.6},{1.0,3.1,1.312});
	Box<3,double> sub_domain({-0.2,0.8,0.2},{1.0,1.1,1.0});

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(0.01);

	vector_dist<3,double,aggregate<double>> vd(0,domain,bc,ghost);

	PointIterator<3,double,CartDecomposition<3,double>> p(vd.getDecomposition(),sz,domain,sub_domain);

	size_t cnt = 0;

	bool good = true;
	while (p.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = p.get().get(0);
		vd.getLastPos()[1] = p.get().get(1);
		vd.getLastPos()[2] = p.get().get(2);

		good &= sub_domain.isInside(p.get());

		cnt++;

		++p;
	}

	Vcluster & v_cl = create_vcluster();

	v_cl.sum(cnt);

	BOOST_REQUIRE_EQUAL(cnt,sz[0]*sz[1]*sz[2]);
	BOOST_REQUIRE_EQUAL(good,true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_UNIT_TESTS_HPP_ */
