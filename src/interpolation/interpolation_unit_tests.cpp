/*
 * interpolation_unit_tests.hpp
 *
 *  Created on: May 5, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_UNIT_TESTS_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interpolation/mp4_kernel.hpp"
#include "interpolation/lambda_kernel.hpp"
#include "interpolation/z_spline.hpp"
#include "interpolation.hpp"
#include <boost/math/special_functions/pow.hpp>
#include <Vector/vector_dist.hpp>
#include <Operators/Vector/vector_dist_operators.hpp>
#include <FiniteDifference/FD_op.hpp>
#include <Grid/grid_dist_id.hpp>

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

template<unsigned int dim, typename T,typename Kernel, typename grid, typename vector>
void interp_test(grid & gd, vector & vd, bool single_particle,unsigned int numberOfMomenta)
{
	// Reset the grid

	auto it2 = gd.getDomainGhostIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		gd.template get<0>(key) = 0.0;

		++it2;
	}

	interpolate<vector,grid,Kernel> inte(vd,gd);

	if (single_particle == false)
	{inte.template p2m<0,0>(vd,gd);}
	else
	{
		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			inte.template p2m<0,0>(vd,gd,p);

			++it;
		}
	}

	T mg[dim];
	T mv[dim];
    vd.write("Particles");
    gd.write("Grid");

    if(numberOfMomenta>=0){
        momenta_grid<grid,0>(gd,mg);
        momenta_vector<vector,0>(vd,mv);
        for (size_t i = 0 ; i < dim ; i++)
        {BOOST_REQUIRE_CLOSE(mg[i],mv[i],0.001);}
	}

	if(numberOfMomenta>=1){
        momenta_grid<grid,1>(gd,mg);
        momenta_vector<vector,1>(vd,mv);
	    for (size_t i = 0 ; i < dim ; i++)
	    {BOOST_REQUIRE_CLOSE(mg[i],mv[i],0.001);}
    }

	if(numberOfMomenta>=2){
        momenta_grid<grid,2>(gd,mg);
        momenta_vector<vector,2>(vd,mv);
        for (size_t i = 0 ; i < dim ; i++)
        {BOOST_REQUIRE_CLOSE(mg[i],mv[i],0.001);}
	}

    if(numberOfMomenta>=3){
        momenta_grid<grid,3>(gd,mg);
        momenta_vector<vector,3>(vd,mv);
        for (size_t i = 0 ; i < dim ; i++)
        {BOOST_REQUIRE_CLOSE(mg[i],mv[i],0.001);}
        }

    if(numberOfMomenta>=4){
        momenta_grid<grid,4>(gd,mg);
        momenta_vector<vector,4>(vd,mv);
        for (size_t i = 0 ; i < dim ; i++)
        {BOOST_REQUIRE_CLOSE(mg[i],mv[i],0.001);}
    }
}



BOOST_AUTO_TEST_CASE( interpolation_full_single_test_2D )
{
	Box<2,float> domain({0.0,0.0},{1.0,1.0});
	size_t sz[2] = {64,64};

	Ghost<2,long int> gg(2);
	Ghost<2,float> gv(0.01);

	size_t bc_v[2] = {PERIODIC,PERIODIC};

	vector_dist<2,float,aggregate<float>> vd(65536,domain,bc_v,gv);
	grid_dist_id<2,float,aggregate<float>> gd(vd.getDecomposition(),sz,gg);

	// set one particle on vd

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (double)rand()/RAND_MAX;
		vd.getPos(p)[1] = (double)rand()/RAND_MAX;

		vd.getProp<0>(p) = 5.0;

		++it;
	}

	vd.map();

	interp_test<2,float,mp4_kernel<float>>(gd,vd,true,2);
}
    BOOST_AUTO_TEST_CASE( interpolation_full_single_test_2D_double )
    {
        Box<2,double> domain({0.0,0.0},{1.0,1.0});
        size_t sz[2] = {64,64};

        Ghost<2,long int> gg(3);
        Ghost<2,double> gv(0.01);

        size_t bc_v[2] = {PERIODIC,PERIODIC};

        vector_dist<2,double,aggregate<double>> vd(65536,domain,bc_v,gv);
        grid_dist_id<2,double,aggregate<double>> gd(vd.getDecomposition(),sz,gg);

        // set one particle on vd

        auto it = vd.getDomainIterator();

        while (it.isNext())
        {
            auto p = it.get();

            vd.getPos(p)[0] = (double)rand()/RAND_MAX;
            vd.getPos(p)[1] = (double)rand()/RAND_MAX;

            vd.getProp<0>(p) = 5.0;

            ++it;
        }

        vd.map();

        interp_test<2,double,lambda4_4kernel<double>>(gd,vd,true,2);
    }


BOOST_AUTO_TEST_CASE( interpolation_full_test_2D )
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

		vd.getProp<0>(p) = 5.0;

		++it;
	}

	vd.map();

	// Reset the grid
	interp_test<2,float,mp4_kernel<float>>(gd,vd,false,2);

	float mg[2];
	float mv[2];

	auto & v_cl = create_vcluster();

	interpolate<decltype(vd),decltype(gd),mp4_kernel<float>> inte(vd,gd);

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

BOOST_AUTO_TEST_CASE( interpolation_full_single_test_3D )
{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t sz[3] = {64,64,64};

	Ghost<3,long int> gg(2);
	Ghost<3,double> gv(0.01);

	size_t bc_v[3] = {PERIODIC,PERIODIC,PERIODIC};

	vector_dist<3,double,aggregate<double>> vd(65536,domain,bc_v,gv);
	grid_dist_id<3,double,aggregate<double>> gd(vd.getDecomposition(),sz,gg);

	// set one particle on vd

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// coverty[dont_call]
		vd.getPos(p)[0] = (double)rand()/RAND_MAX;
		// coverty[dont_call]
		vd.getPos(p)[1] = (double)rand()/RAND_MAX;
		// coverty[dont_call]
		vd.getPos(p)[2] = (double)rand()/RAND_MAX;

		vd.getProp<0>(p) = 5.0;

		++it;
	}

	vd.map();

	// Reset the grid
	interp_test<3,double,mp4_kernel<float>>(gd,vd,true,2);
}

BOOST_AUTO_TEST_CASE( interpolation_getSubCheck )
{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t sz[3] = {64,64,64};

	Ghost<3,long int> gg(2);
	Ghost<3,double> gv(0.01);

	size_t bc_v[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	typedef vector_dist<3,double,aggregate<double>> vector;
	typedef grid_dist_id<3,double,aggregate<double>> grid;

	vector vd(0,domain,bc_v,gv);
	grid_dist_id<3,double,aggregate<double>> gd(vd.getDecomposition(),sz,gg);

	// interpolate
	interpolate<vector,grid,mp4_kernel<double>> inte(vd,gd);

	// decomposition
	auto & dec = vd.getDecomposition();

	int nl = dec.getNLocalSub();

	for (int i = 0 ; i < nl ; i++)
	{
		int  nll = dec.getLocalNIGhost(i);
		for (int j = 0 ; j < nll ; j++)
		{
			auto ibx = dec.getLocalIGhostBox(i,j);
			int x = dec.getLocalIGhostSub(i,j);
			auto bx = dec.getSubDomain(x);

			Point<3,double> p;
			for (int s = 0; s < 3 ; s++)
			{
				Point<3,double> p;
				for (int s1 = 0; s1 < 3 ; s1++)
				{
					p.get(s1) = (ibx.getHigh(s1) - ibx.getLow(s1)) / 2.0 + ibx.getLow(s1);
				}

				if (ibx.getLow(s) == bx.getHigh(s))
				{
					p.get(s) = ibx.getLow(s);
					int sub = inte.getSub(p);
					BOOST_REQUIRE_EQUAL(sub,i);
				}
				else if (ibx.getHigh(s) == bx.getLow(s))
				{
					p.get(s) = ibx.getHigh(s);
					int sub = inte.getSub(p);
					BOOST_REQUIRE_EQUAL(sub,x);
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( interpolation_full_test_3D )
{
	Box<3,double> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t sz[3] = {64,64,64};

	Ghost<3,long int> gg(2);
	Ghost<3,double> gv(0.01);

	size_t bc_v[3] = {PERIODIC,PERIODIC,PERIODIC};

	{
	vector_dist<3,double,aggregate<double>> vd(65536,domain,bc_v,gv);
	grid_dist_id<3,double,aggregate<double>> gd(vd.getDecomposition(),sz,gg);

	// set one particle on vd

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (double)rand()/RAND_MAX;
		vd.getPos(p)[1] = (double)rand()/RAND_MAX;
		vd.getPos(p)[2] = (double)rand()/RAND_MAX;

		vd.getProp<0>(p) = 5.0;

		++it;
	}

	vd.map();

	// Reset the grid

	// Reset the grid
	interp_test<3,double,mp4_kernel<float>>(gd,vd,false,2);

	auto & v_cl = create_vcluster();
	double mg[3];
	double mv[3];
	interpolate<decltype(vd),decltype(gd),mp4_kernel<double>> inte(vd,gd);

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
		vd.getLastPos()[2] = key.get(2) * it4.getSpacing(2) + domain.getLow(2) + 0.1*it4.getSpacing(2);

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

	grid_key_dx<3> start({3,3,3});
	grid_key_dx<3> stop({(long int)gd.size(0) - 4,(long int)gd.size(1) - 4,(long int)gd.size(2) - 4});

	auto it6 = gd.getSubDomainIterator(start,stop);
	while(it6.isNext())
	{
		auto key = it6.get();

		gd.get<0>(key) = 5.0;

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
	v_cl.sum(mg[2]);
	v_cl.sum(mv[0]);
	v_cl.sum(mv[1]);
	v_cl.sum(mv[2]);
	v_cl.execute();

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);
	BOOST_REQUIRE_CLOSE(mg[2],mv[2],0.001);

	momenta_grid_domain<decltype(gd),1>(gd,mg);
	momenta_vector<decltype(vd),1>(vd,mv);

	v_cl.sum(mg[0]);
	v_cl.sum(mg[1]);
	v_cl.sum(mg[2]);
	v_cl.sum(mv[0]);
	v_cl.sum(mv[1]);
	v_cl.sum(mv[2]);
	v_cl.execute();

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);
	BOOST_REQUIRE_CLOSE(mg[2],mv[2],0.001);

	momenta_grid_domain<decltype(gd),2>(gd,mg);
	momenta_vector<decltype(vd),2>(vd,mv);

	v_cl.sum(mg[0]);
	v_cl.sum(mg[1]);
	v_cl.sum(mg[2]);
	v_cl.sum(mv[0]);
	v_cl.sum(mv[1]);
	v_cl.sum(mv[2]);
	v_cl.execute();

	BOOST_REQUIRE_CLOSE(mg[0],mv[0],0.001);
	BOOST_REQUIRE_CLOSE(mg[1],mv[1],0.001);
	BOOST_REQUIRE_CLOSE(mg[2],mv[2],0.001);

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

BOOST_AUTO_TEST_CASE( int_kernel_test_double)
{
		mp4_kernel<double> mp4;

		double tot = 0.0;

		// Check momenta 0

		tot += mp4.value(-1.3,0);
		tot += mp4.value(-0.3,1);
		tot += mp4.value(0.7,2);
		tot += mp4.value(1.7,3);

		BOOST_REQUIRE_CLOSE(tot,1.0,0.001);

		// Check momenta 1

		tot = 0.0;

		tot += -1.3*mp4.value(-1.3,0);
		tot += -0.3*mp4.value(-0.3,1);
		tot += 0.7*mp4.value(0.7,2);
		tot += 1.7*mp4.value(1.7,3);

		BOOST_REQUIRE_SMALL(tot,0.001);

		// Check momenta 2

		tot = 0.0;

		tot += (1.3)*(1.3)*mp4.value(-1.3,0);
		tot += (0.3)*(0.3)*mp4.value(-0.3,1);
		tot += (0.7)*(0.7)*mp4.value(0.7,2);
		tot += (1.7)*(1.7)*mp4.value(1.7,3);

		BOOST_REQUIRE_SMALL(tot,0.001);


		//////// Check zeta 1

		tot = 0.0;

		z_kernel<double,1> zk1;

		tot += zk1.value(-0.3,0);
		tot += zk1.value(0.7,1);

		BOOST_REQUIRE_CLOSE(tot,1.0,0.001);

		//////// zeta 2 is equivalent to mp4 we do not test

		//////// zeta 3

		z_kernel<double,3> zk3;

		tot = 0.0;

		// Check momenta 0

		tot += zk3.value(-2.3,0);
		tot += zk3.value(-1.3,1);
		tot += zk3.value(-0.3,2);
		tot += zk3.value(0.7,3);
		tot += zk3.value(1.7,4);
		tot += zk3.value(2.7,5);

		BOOST_REQUIRE_CLOSE(tot,1.0,0.001);

		// Check momenta 1

		tot = 0.0;

		tot += -2.3*zk3.value(-2.3,0);
		tot += -1.3*zk3.value(-1.3,1);
		tot += -0.3*zk3.value(-0.3,2);
		tot += 0.7*zk3.value(0.7,3);
		tot += 1.7*zk3.value(1.7,4);
		tot += 2.7*zk3.value(2.7,5);

		BOOST_REQUIRE_SMALL(tot,0.001);

		// Check momenta 2

		tot = 0.0;

		tot += 2.3*2.3*zk3.value(-2.3,0);
		tot += 1.3*1.3*zk3.value(-1.3,1);
		tot += 0.3*0.3*zk3.value(-0.3,2);
		tot += 0.7*0.7*zk3.value(0.7,3);
		tot += 1.7*1.7*zk3.value(1.7,4);
		tot += 2.7*2.7*zk3.value(2.7,5);

		BOOST_REQUIRE_SMALL(tot,0.001);

		// Check momenta 3

		tot = 0.0;

		tot += -2.3*-2.3*-2.3*zk3.value(-2.3,0);
		tot += -1.3*-1.3*-1.3*zk3.value(-1.3,1);
		tot += -0.3*-0.3*-0.3*zk3.value(-0.3,2);
		tot += 0.7*0.7*0.7*zk3.value(0.7,3);
		tot += 1.7*1.7*1.7*zk3.value(1.7,4);
		tot += 2.7*2.7*2.7*zk3.value(2.7,5);

		BOOST_REQUIRE_SMALL(tot,0.001);


		// z4

		z_kernel<double,4> zk4;

		// Check momenta 0

		tot = 0.0;

		tot += zk4.value(-3.3,0);
		tot += zk4.value(-2.3,1);
		tot += zk4.value(-1.3,2);
		tot += zk4.value(-0.3,3);
		tot += zk4.value(0.7,4);
		tot += zk4.value(1.7,5);
		tot += zk4.value(2.7,6);
		tot += zk4.value(3.7,7);

		BOOST_REQUIRE_CLOSE(tot,1.0,0.001);

		// Check momenta 1

		tot = 0.0;

		tot += -3.3*zk4.value(-3.3,0);
		tot += -2.3*zk4.value(-2.3,1);
		tot += -1.3*zk4.value(-1.3,2);
		tot += -0.3*zk4.value(-0.3,3);
		tot += 0.7*zk4.value(0.7,4);
		tot += 1.7*zk4.value(1.7,5);
		tot += 2.7*zk4.value(2.7,6);
		tot += 3.7*zk4.value(3.7,7);

		BOOST_REQUIRE_SMALL(tot,0.001);

		// Check momenta 2

		tot = 0.0;

		tot += 3.3*3.3*zk4.value(-3.3,0);
		tot += 2.3*2.3*zk4.value(-2.3,1);
		tot += 1.3*1.3*zk4.value(-1.3,2);
		tot += 0.3*0.3*zk4.value(-0.3,3);
		tot += 0.7*0.7*zk4.value(0.7,4);
		tot += 1.7*1.7*zk4.value(1.7,5);
		tot += 2.7*2.7*zk4.value(2.7,6);
		tot += 3.7*3.7*zk4.value(3.7,7);

		BOOST_REQUIRE_SMALL(tot,0.001);

		// Check momenta 3

		tot = 0.0;

		tot += -3.3*-3.3*-3.3*zk4.value(-3.3,0);
		tot += -2.3*-2.3*-2.3*zk4.value(-2.3,1);
		tot += -1.3*-1.3*-1.3*zk4.value(-1.3,2);
		tot += -0.3*-0.3*-0.3*zk4.value(-0.3,3);
		tot += 0.7*0.7*0.7*zk4.value(0.7,4);
		tot += 1.7*1.7*1.7*zk4.value(1.7,5);
		tot += 2.7*2.7*2.7*zk4.value(2.7,6);
		tot += 3.7*3.7*3.7*zk4.value(3.7,7);

		BOOST_REQUIRE_SMALL(tot,0.001);

		// Check momenta 4

		tot = 0.0;

		tot += -3.3*-3.3*-3.3*-3.3*zk4.value(-3.3,0);
		tot += -2.3*-2.3*-2.3*-2.3*zk4.value(-2.3,1);
		tot += -1.3*-1.3*-1.3*-1.3*zk4.value(-1.3,2);
		tot += -0.3*-0.3*-0.3*-0.3*zk4.value(-0.3,3);
		tot += 0.7*0.7*0.7*0.7*zk4.value(0.7,4);
		tot += 1.7*1.7*1.7*1.7*zk4.value(1.7,5);
		tot += 2.7*2.7*2.7*2.7*zk4.value(2.7,6);
		tot += 3.7*3.7*3.7*3.7*zk4.value(3.7,7);

		BOOST_REQUIRE_SMALL(tot,0.001);
}

/*
BOOST_AUTO_TEST_CASE(InterpolationConvergenceP2M)
{
        size_t res;
        std::cout<<"Enter Res:";
        std::cin>>res;
        const size_t sz[2] = {res,res};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {PERIODIC, PERIODIC};
        double spacing[2];
        spacing[0] = 2 *M_PI / (sz[0]);
        spacing[1] = 2 *M_PI / (sz[1]);
        Ghost<2,long int> gg(3);
        double rCut = 3.0 * spacing[0];
        Ghost<2, double> ghost(rCut);

        vector_dist<2, double, aggregate<double, double>> particles(0, box,bc,ghost);
        grid_dist_id<2, double, aggregate<double, double>> gd(particles.getDecomposition(),sz,gg);
        double sigma2 = spacing[0] / (40.0);
        std::normal_distribution<> gaussian{0, sigma2};
        std::mt19937 rng{6666666};
        auto it = particles.getGridIterator(sz);
        while (it.isNext()) {
            particles.add();
            auto key = it.get();
            double x=key.get(0) * spacing[0] + gaussian(rng);
            double y=key.get(1) * spacing[1] + gaussian(rng);
            particles.getLastPos()[0] = x;
            particles.getLastPos()[1] = y;
            // Here fill the function value
            particles.template getLastProp<0>() = sin(particles.getLastPos()[0]) + sin(particles.getLastPos()[0]);
            ++it;
        }
        particles.map();
        particles.ghost_get<0>();

        auto itG=gd.getDomainIterator();
        while(itG.isNext())
        {
            auto key=itG.get();
            gd.template getProp<1>(key) = sin(gd.getPos(key)[0]) + sin(gd.getPos(key)[0]);
            ++itG;
        }

        particles.write("InitP");
        gd.write("Grid");

        auto Pu=getV<0>(particles);
        auto Gu=FD::getV<0>(gd);
        typedef vector_dist<2, double, aggregate<double, double>> particle_type;
        typedef grid_dist_id<2, double, aggregate<double, double>> gd_type;
        typedef z_kernel<double,4> kerneltype; //mp4_kernel<double>
        typedef lambda4_4kernel<double> kerneltype2;
        interpolate<particle_type,gd_type,kerneltype2> inte2m(particles,gd);
        Gu=0;
        gd.ghost_get<0>();
        inte2m.template p2m<0,0>(particles,gd);
        gd.template ghost_put<add_,0>();
        gd.ghost_get<0>();
        particles.write("InitPAfter");
        gd.write("GridAfter");


        auto it2 = gd.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(gd.template getProp<1>(p) - gd.template  getProp<0>(p)) > worst) {
                worst = fabs(gd. template  getProp<1>(p) - gd.template  getProp<0>(p));
            }
            ++it2;
        }
        std::cout<<worst<<std::endl;
        //BOOST_REQUIRE(worst < 0.03);
}

BOOST_AUTO_TEST_CASE(InterpolationConvergenceM2P)
{
        size_t res;
        std::cout<<"Enter Res:";
        std::cin>>res;
        const size_t sz[2] = {res,res};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {PERIODIC, PERIODIC};
        double spacing[2];
        spacing[0] = 2 *M_PI / (sz[0]);
        spacing[1] = 2 *M_PI / (sz[1]);
        Ghost<2,long int> gg(3);
        double rCut = 3.0 * spacing[0];
        Ghost<2, double> ghost(rCut);

        vector_dist<2, double, aggregate<double, double>> particles(0, box,bc,ghost);
        grid_dist_id<2, double, aggregate<double, double>> gd(particles.getDecomposition(),sz,gg);
        double sigma2 = spacing[0] * spacing[1];
        std::normal_distribution<> gaussian{0, sigma2};
        std::mt19937 rng{6666666};
        auto it = particles.getGridIterator(sz);
        while (it.isNext()) {
            particles.add();
            auto key = it.get();
            double x=key.get(0) * spacing[0] + gaussian(rng);
            double y=key.get(1) * spacing[1] + gaussian(rng);
            particles.getLastPos()[0] = x;
            particles.getLastPos()[1] = y;
            // Here fill the function value
            particles.template getLastProp<0>() = 0;
            particles.template getLastProp<1>() = sin(particles.getLastPos()[0]) + sin(particles.getLastPos()[0]);
            ++it;
        }
        particles.map();
        particles.ghost_get<0>();

        auto itG=gd.getDomainIterator();
        while(itG.isNext())
        {
            auto key=itG.get();
            gd.template getProp<1>(key) = sin(gd.getPos(key)[0]) + sin(gd.getPos(key)[0]);
            ++itG;
        }

        particles.write("InitP");
        gd.write("Grid");

        auto Pu=getV<0>(particles);
        auto Gu=FD::getV<0>(gd);
        typedef vector_dist<2, double, aggregate<double, double>> particle_type;
        typedef grid_dist_id<2, double, aggregate<double, double>> gd_type;
        typedef z_kernel<double,4> kerneltype; //mp4_kernel<double>
        typedef lambda4_4kernel<double> kerneltype2;
        interpolate<particle_type,gd_type,kerneltype2> inte2m(particles,gd);
        Gu=0;
        gd.ghost_get<0>();
        inte2m.template m2p<1,0>(gd,particles);
        particles.ghost_get<0>();
        particles.write("InitPAfter");
        gd.write("GridAfter");

        auto it2 = particles.getDomainIterator();

        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(particles.template getProp<1>(p) - particles.template  getProp<0>(p)) > worst) {
                worst = fabs(particles. template  getProp<1>(p) - particles.template  getProp<0>(p));
            }
            ++it2;
        }
        std::cout<<worst<<std::endl;
        //BOOST_REQUIRE(worst < 0.03);

}



BOOST_AUTO_TEST_CASE(InterpolationMoving)
{

        size_t res;
        std::cin>>res;
        const size_t sz[2] = {res,res};
        Box<2, double> box({0, 0}, {2 * M_PI, 2 * M_PI});
        size_t bc[2] = {PERIODIC, PERIODIC};
        double spacing[2];
        spacing[0] = 2 *M_PI / (sz[0]);
        spacing[1] = 2 *M_PI / (sz[1]);
        Ghost<2,long int> gg(3);
        double rCut = 3.0 * spacing[0];
        Ghost<2, double> ghost(rCut);

        vector_dist<2, double, aggregate<double, VectorS<2, double>>> particles(0, box,bc,ghost),particlesMoved(0, box,bc,ghost);
        grid_dist_id<2, double, aggregate<double, VectorS<2, double>>> gd(particles.getDecomposition(),sz,gg);

        auto it = particles.getGridIterator(sz);
        while (it.isNext()) {
            particles.add();
            auto key = it.get();
            double x=key.get(0) * spacing[0];
            double y=key.get(1) * spacing[1];
            particles.getLastPos()[0] = x;
            particles.getLastPos()[1] = y;
            // Here fill the function value
            particles.template getLastProp<1>() = 1.0;
            particles.template getLastProp<1>() = 0;
            particles.template getLastProp<0>() = 0;
            if((x-3.14)*(x-3.14)+(y-3.14)*(y-3.14)<1)
            {
                particles.template getLastProp<0>() = 1;
            }
            ++it;
        }
        particles.map();
        particles.ghost_get<0>();

        particles.write("InitP");
        gd.write("Grid");

        auto Pu=getV<0>(particles);
        auto Pmu=getV<0>(particlesMoved);
        auto Gu=FD::getV<0>(gd);
        typedef vector_dist<2, double, aggregate<double, VectorS<2, double>>> vd;
        typedef grid_dist_id<2, double, aggregate<double, VectorS<2, double>>> gd_type;
        interpolate<vd,gd_type,mp4_kernel<double>> inte2m(particlesMoved,gd);
        interpolate<vd,gd_type,mp4_kernel<double>> inte2p(particles,gd);
        double t=0,dt=0.5;
        int ctr=0;
        while(t<10)
        {
            particlesMoved.clear();
            auto it=particles.getDomainIterator();
            while(it.isNext())
            {
                auto p=it.get();
                double xp=particles.getPos(p)[0],yp=particles.getPos(p)[1];
                particlesMoved.add();
                particlesMoved.getLastPos()[0] = xp+dt*particles.getProp<1>(p)[0];
                particlesMoved.getLastPos()[1] = yp+dt*particles.getProp<1>(p)[1];
                particlesMoved.getLastProp<0>() = particles.getProp<0>(p);
                ++it;
            }
            particlesMoved.map();
            particlesMoved.ghost_get<0>();
            Gu=0;
            gd.ghost_get<0>();
            inte2m.template p2m<0,0>(particlesMoved,gd);
            gd.template ghost_put<add_,0>();
            gd.ghost_get<0>();
            Pu=0;
            inte2p.template m2p<0,0>(gd,particles);
            particles.write_frame("InitP",ctr);
            gd.write_frame("Grid",ctr);
            ctr++;
            t+=dt;
        }*/

/*

        auto it2 = domain.getDomainIterator();

        double worst = 0.0;

        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) > worst) {
                worst = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
            }
            ++it2;
        }
*/

 //       domain.deleteGhost();
 //       BOOST_REQUIRE(worst < 0.03);

//}
BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_INTERPOLATION_INTERPOLATION_UNIT_TESTS_HPP_ */
