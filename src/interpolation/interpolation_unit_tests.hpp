/*
 * interpolation_unit_tests.hpp
 *
 *  Created on: May 5, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_UNIT_TESTS_HPP_

#include "interpolation/mp4_kernel.hpp"
#include "interpolation/interpolation.hpp"
#include "interpolation/z_spline.hpp"

BOOST_AUTO_TEST_SUITE( interpolation_test )

template<typename grid, unsigned int mom_p> void momenta_grid(grid & gd,typename grid::stype (& mom_tot)[grid::dims])
{
	auto it = gd.getDomainGhostIterator();

	for (size_t i = 0 ; i < grid::dims ; i++)
		mom_tot[i] = 0.0;

	while (it.isNext())
	{
		auto key = it.get();
		auto key_g = gd.getGKey(key);

		for (size_t i = 0 ; i < grid::dims ; i++)
		{
			typename grid::stype coord = gd.spacing(i)*key_g.get(i);

			mom_tot[i] += boost::math::pow<mom_p>(coord)*gd.template getProp<0>(key);
		}

		++it;
	}
}

template<typename grid, unsigned int mom_p> void momenta_grid_domain(grid & gd,typename grid::stype (& mom_tot)[grid::dims])
{
	auto it = gd.getDomainIterator();

	for (size_t i = 0 ; i < grid::dims ; i++)
		mom_tot[i] = 0.0;

	while (it.isNext())
	{
		auto key = it.get();
		auto key_g = gd.getGKey(key);

		for (size_t i = 0 ; i < grid::dims ; i++)
		{
			typename grid::stype coord = gd.spacing(i)*key_g.get(i);

			mom_tot[i] += boost::math::pow<mom_p>(coord)*gd.template getProp<0>(key);
		}

		++it;
	}
}

template<typename vector, unsigned int mom_p> void momenta_vector(vector & vd,typename vector::stype (& mom_tot)[vector::dims])
{
	auto it = vd.getDomainIterator();

	for (size_t i = 0 ; i < vector::dims ; i++)
		mom_tot[i] = 0.0;

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0 ; i < vector::dims ; i++)
		{
			typename vector::stype coord = vd.getPos(key)[i];

			mom_tot[i] += boost::math::pow<mom_p>(coord)*vd.template getProp<0>(key);
		}

		++it;
	}
}


BOOST_AUTO_TEST_CASE( interpolation_full_test )
{
	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	size_t sz[2] = {64,64};

	Ghost<2,long int> gg(3);
	Ghost<2,float> gv(0.01);

	size_t bc_v[2] = {PERIODIC,PERIODIC};

	{
	vector_dist<2,float,aggregate<float>> vd(4096,domain,bc_v,gv);
	grid_dist_id<2,float,aggregate<float>> gd(vd.getDecomposition(),sz,gg);

	// set one particle on vd

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (double)rand()/RAND_MAX;
		vd.getPos(p)[1] = (double)rand()/RAND_MAX;
		vd.getPos(p)[2] = (double)rand()/RAND_MAX;

		vd.getProp<0>(p) = 5.0/*(double)rand()/RAND_MAX*/;

		++it;
	}

	vd.map();

	// Reset the grid

	auto it2 = gd.getDomainGhostIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		gd.template get<0>(key) = 0.0;

		++it2;
	}

	interpolate<decltype(vd),decltype(gd),mp4_kernel<float>> inte(vd,gd);

	inte.p2m<0,0>(vd,gd);

	float mg[2];
	float mv[2];

	momenta_grid<decltype(gd),0>(gd,mg);
	momenta_vector<decltype(vd),0>(vd,mv);

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);

	momenta_grid<decltype(gd),1>(gd,mg);
	momenta_vector<decltype(vd),1>(vd,mv);

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);

	momenta_grid<decltype(gd),2>(gd,mg);
	momenta_vector<decltype(vd),2>(vd,mv);

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);

	auto & v_cl = create_vcluster();

	// We have to do a ghost get before interpolating m2p
	// Before doing mesh to particle particle must be arranged
	// into a grid like

	vd.clear();

	auto it4 = vd.getGridIterator(sz);

	while(it4.isNext())
	{
		auto key = it4.get();

		vd.add();

		vd.getLastPos()[0] = key.get(0) * it4.getSpacing(0) + domain.getLow(0) + 0.1*it4.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it4.getSpacing(1) + domain.getLow(1) + 0.1*it4.getSpacing(1);

		vd.getLastProp<0>() = 0.0;

		++it4;
	}

	// Reset also the grid

	auto it5 = gd.getDomainGhostIterator();

	while(it5.isNext())
	{
		auto key = it5.get();

		gd.get<0>(key) = 0.0;

		++it5;
	}
	gd.ghost_get<0>();

	grid_key_dx<2> start({3,3});
	grid_key_dx<2> stop({(long int)gd.size(0) - 4,(long int)gd.size(1) - 4});

	auto it6 = gd.getSubDomainIterator(start,stop);
	while(it6.isNext())
	{
		auto key = it6.get();

		gd.get<0>(key) = 5.0/*(double)rand()/RAND_MAX;*/;

		++it6;
	}
	gd.ghost_get<0>();

	vd.map();
	gd.ghost_get<0>();
	inte.m2p<0,0>(gd,vd);

	momenta_grid_domain<decltype(gd),0>(gd,mg);
	momenta_vector<decltype(vd),0>(vd,mv);

	v_cl.sum(mg[0]);
	v_cl.sum(mg[1]);
	v_cl.sum(mv[0]);
	v_cl.sum(mv[1]);
	v_cl.execute();

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);

	momenta_grid_domain<decltype(gd),1>(gd,mg);
	momenta_vector<decltype(vd),1>(vd,mv);

	v_cl.sum(mg[0]);
	v_cl.sum(mg[1]);
	v_cl.sum(mv[0]);
	v_cl.sum(mv[1]);
	v_cl.execute();

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);

	momenta_grid_domain<decltype(gd),2>(gd,mg);
	momenta_vector<decltype(vd),2>(vd,mv);

	v_cl.sum(mg[0]);
	v_cl.sum(mg[1]);
	v_cl.sum(mv[0]);
	v_cl.sum(mv[1]);
	v_cl.execute();

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);

	}
}

BOOST_AUTO_TEST_CASE( int_kernel_test )
{
		mp4_kernel<float> mp4;

		float tot = 0.0;

		// Check momenta 0

		tot += mp4.value(-1.3f,0);
		tot += mp4.value(-0.3f,1);
		tot += mp4.value(0.7f,2);
		tot += mp4.value(1.7f,3);

		BOOST_REQUIRE_CLOSE(tot,1.0f,0.001);

		// Check momenta 1

		tot = 0.0;

		tot += -1.3f*mp4.value(-1.3f,0);
		tot += -0.3f*mp4.value(-0.3f,1);
		tot += 0.7f*mp4.value(0.7f,2);
		tot += 1.7f*mp4.value(1.7f,3);

		BOOST_REQUIRE_SMALL(tot,0.001f);

		// Check momenta 2

		tot = 0.0;

		tot += (1.3f)*(1.3f)*mp4.value(-1.3f,0);
		tot += (0.3f)*(0.3f)*mp4.value(-0.3f,1);
		tot += (0.7f)*(0.7f)*mp4.value(0.7f,2);
		tot += (1.7f)*(1.7f)*mp4.value(1.7f,3);

		BOOST_REQUIRE_SMALL(tot,0.001f);


		//////// Check zeta 1

		tot = 0.0;

		z_kernel<float,1> zk1;

		tot += zk1.value(-0.3f,0);
		tot += zk1.value(0.7f,1);

		BOOST_REQUIRE_CLOSE(tot,1.0f,0.001);

		//////// zeta 2 is equivalent to mp4 we do not test

		//////// zeta 3

		z_kernel<float,3> zk3;

		tot = 0.0;

		// Check momenta 0

		tot += zk3.value(-2.3f,0);
		tot += zk3.value(-1.3f,1);
		tot += zk3.value(-0.3f,2);
		tot += zk3.value(0.7f,3);
		tot += zk3.value(1.7f,4);
		tot += zk3.value(2.7f,5);

		BOOST_REQUIRE_CLOSE(tot,1.0f,0.001);

		// Check momenta 1

		tot = 0.0;

		tot += -2.3*zk3.value(-2.3f,0);
		tot += -1.3*zk3.value(-1.3f,1);
		tot += -0.3*zk3.value(-0.3f,2);
		tot += 0.7*zk3.value(0.7f,3);
		tot += 1.7*zk3.value(1.7f,4);
		tot += 2.7*zk3.value(2.7f,5);

		BOOST_REQUIRE_SMALL(tot,0.001f);

		// Check momenta 2

		tot = 0.0;

		tot += 2.3*2.3*zk3.value(-2.3f,0);
		tot += 1.3*1.3*zk3.value(-1.3f,1);
		tot += 0.3*0.3*zk3.value(-0.3f,2);
		tot += 0.7*0.7*zk3.value(0.7f,3);
		tot += 1.7*1.7*zk3.value(1.7f,4);
		tot += 2.7*2.7*zk3.value(2.7f,5);

		BOOST_REQUIRE_SMALL(tot,0.001f);

		// Check momenta 3

		tot = 0.0;

		tot += -2.3*-2.3*-2.3*zk3.value(-2.3f,0);
		tot += -1.3*-1.3*-1.3*zk3.value(-1.3f,1);
		tot += -0.3*-0.3*-0.3*zk3.value(-0.3f,2);
		tot += 0.7*0.7*0.7*zk3.value(0.7f,3);
		tot += 1.7*1.7*1.7*zk3.value(1.7f,4);
		tot += 2.7*2.7*2.7*zk3.value(2.7f,5);

		BOOST_REQUIRE_SMALL(tot,0.001f);


		// z4

		z_kernel<float,4> zk4;

		// Check momenta 0

		tot = 0.0;

		tot += zk4.value(-3.3f,0);
		tot += zk4.value(-2.3f,1);
		tot += zk4.value(-1.3f,2);
		tot += zk4.value(-0.3f,3);
		tot += zk4.value(0.7f,4);
		tot += zk4.value(1.7f,5);
		tot += zk4.value(2.7f,6);
		tot += zk4.value(3.7f,7);

		BOOST_REQUIRE_CLOSE(tot,1.0f,0.001);

		// Check momenta 1

		tot = 0.0;

		tot += -3.3*zk4.value(-3.3f,0);
		tot += -2.3*zk4.value(-2.3f,1);
		tot += -1.3*zk4.value(-1.3f,2);
		tot += -0.3*zk4.value(-0.3f,3);
		tot += 0.7*zk4.value(0.7f,4);
		tot += 1.7*zk4.value(1.7f,5);
		tot += 2.7*zk4.value(2.7f,6);
		tot += 3.7*zk4.value(3.7f,7);

		BOOST_REQUIRE_SMALL(tot,0.001f);

		// Check momenta 2

		tot = 0.0;

		tot += 3.3*3.3*zk4.value(-3.3f,0);
		tot += 2.3*2.3*zk4.value(-2.3f,1);
		tot += 1.3*1.3*zk4.value(-1.3f,2);
		tot += 0.3*0.3*zk4.value(-0.3f,3);
		tot += 0.7*0.7*zk4.value(0.7f,4);
		tot += 1.7*1.7*zk4.value(1.7f,5);
		tot += 2.7*2.7*zk4.value(2.7f,6);
		tot += 3.7*3.7*zk4.value(3.7f,7);

		BOOST_REQUIRE_SMALL(tot,0.001f);

		// Check momenta 3

		tot = 0.0;

		tot += -3.3*-3.3*-3.3*zk4.value(-3.3f,0);
		tot += -2.3*-2.3*-2.3*zk4.value(-2.3f,1);
		tot += -1.3*-1.3*-1.3*zk4.value(-1.3f,2);
		tot += -0.3*-0.3*-0.3*zk4.value(-0.3f,3);
		tot += 0.7*0.7*0.7*zk4.value(0.7f,4);
		tot += 1.7*1.7*1.7*zk4.value(1.7f,5);
		tot += 2.7*2.7*2.7*zk4.value(2.7f,6);
		tot += 3.7*3.7*3.7*zk4.value(3.7f,7);

		BOOST_REQUIRE_SMALL(tot,0.001f);

		// Check momenta 4

		tot = 0.0;

		tot += -3.3*-3.3*-3.3*-3.3*zk4.value(-3.3f,0);
		tot += -2.3*-2.3*-2.3*-2.3*zk4.value(-2.3f,1);
		tot += -1.3*-1.3*-1.3*-1.3*zk4.value(-1.3f,2);
		tot += -0.3*-0.3*-0.3*-0.3*zk4.value(-0.3f,3);
		tot += 0.7*0.7*0.7*0.7*zk4.value(0.7f,4);
		tot += 1.7*1.7*1.7*1.7*zk4.value(1.7f,5);
		tot += 2.7*2.7*2.7*2.7*zk4.value(2.7f,6);
		tot += 3.7*3.7*3.7*3.7*zk4.value(3.7f,7);

		BOOST_REQUIRE_SMALL(tot,0.001f);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_UNIT_TESTS_HPP_ */
