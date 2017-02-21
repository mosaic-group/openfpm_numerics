/*
 * DrawParticles_unit_tests.hpp
 *
 *  Created on: Jan 3, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_UNIT_TESTS_HPP_

#include "PointIterator.hpp"
#include "PointIteratorSkin.hpp"
#include "DrawParticles.hpp"

BOOST_AUTO_TEST_SUITE( draw_particles )

BOOST_AUTO_TEST_CASE(point_iterator)
{
	//! [DrawBox_example]

	size_t sz[] = {23,27,20};

	Box<3,double> domain({-1.2,0.5,-0.6},{1.0,3.1,1.3});
	Box<3,double> sub_domain({-0.15,0.75,0.15},{1.05,1.15,1.05});
	Box<3,double> sub_real({-0.1,0.80,0.2},{1.0,1.1,1.0});

	size_t sz_sub[] = {12,4,9};

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(0.01);

	vector_dist<3,double,aggregate<double>> vd(0,domain,bc,ghost);

	auto p = DrawParticles::DrawBox(vd,sz,domain,sub_domain);

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

	//! [DrawBox_example]

	Vcluster & v_cl = create_vcluster();

	v_cl.sum(cnt);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(cnt,sz_sub[0]*sz_sub[1]*sz_sub[2]);
	BOOST_REQUIRE_EQUAL(good,true);

	Box<3,double> marg = p.getBoxMargins();

	BOOST_REQUIRE_CLOSE(marg.getLow(0), sub_real.getLow(0) ,0.0001);
	BOOST_REQUIRE_CLOSE(marg.getLow(1), sub_real.getLow(1) ,0.0001);
	BOOST_REQUIRE_CLOSE(marg.getLow(2), sub_real.getLow(2) ,0.0001);

	BOOST_REQUIRE_CLOSE(marg.getHigh(0), sub_real.getHigh(0) ,0.0001);
	BOOST_REQUIRE_CLOSE(marg.getHigh(1), sub_real.getHigh(1) ,0.0001);
	BOOST_REQUIRE_CLOSE(marg.getHigh(2), sub_real.getHigh(2) ,0.0001);
}

BOOST_AUTO_TEST_CASE(point_iterator_skin)
{
	size_t sz[] = {23,27,20};

	Box<3,double> domain({-1.2,0.5,-0.6},{1.0,3.1,1.3});
	Box<3,double> sub_domainA({-0.15,0.75,0.15},{1.05,1.15,1.05});
	Box<3,double> sub_domainB({-0.25,0.65,0.05},{0.95,1.05,1.05});

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(0.01);

	vector_dist<3,double,aggregate<double>> vd(0,domain,bc,ghost);

	auto p = DrawParticles::DrawSkin(vd,sz,domain,sub_domainA, sub_domainB);

	size_t cnt = 0;

	bool good = true;
	while (p.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = p.get().get(0);
		vd.getLastPos()[1] = p.get().get(1);
		vd.getLastPos()[2] = p.get().get(2);

		if (!((p.get().get(0) > -0.21 && p.get().get(0) < -0.19) || (p.get().get(1) > 0.69 && p.get().get(1) < 0.71) || (p.get().get(2) > 0.09 && p.get().get(2) < 0.11)) )
			good &= false;

		if (sub_domainA.isInside(p.get()))
			good &= false;

		cnt++;

		++p;
	}

	Vcluster & v_cl = create_vcluster();

	v_cl.sum(cnt);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(good,true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_UNIT_TESTS_HPP_ */
