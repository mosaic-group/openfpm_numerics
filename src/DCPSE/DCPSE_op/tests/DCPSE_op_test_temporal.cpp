/*
 * DCPSE_op_test_temporal.cpp
 *
 *  Created on: Sep 15, 2020
 *      Author: i-bird
 */

#ifndef DCPSE_OP_TEST_TEMPORAL_CPP_
#define DCPSE_OP_TEST_TEMPORAL_CPP_

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "../DCPSE_op.hpp"

BOOST_AUTO_TEST_SUITE(temporal_test_suite)

    BOOST_AUTO_TEST_CASE(temporal_test)
    {
		size_t grd_sz = 17;
		double boxsize = 10;
		const size_t sz[3] = {grd_sz, grd_sz, grd_sz};
		Box<3, double> box({0, 0, 0}, {boxsize, boxsize, boxsize});
		size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
		double Lx = box.getHigh(0);
		double Ly = box.getHigh(1);
		double Lz = box.getHigh(2);
		double spacing = box.getHigh(0) / (sz[0] - 1);
		double rCut = 3.9 * spacing;
		double rCut2 = 3.9 * spacing;
		int ord = 2;
		int ord2 = 2;
		double sampling_factor = 4.0;
		double sampling_factor2 = 2.4;
		Ghost<3, double> ghost(rCut);
		auto &v_cl = create_vcluster();

		vector_dist<3, double, aggregate<double,VectorS<3, double>,double[3][3],double,VectorS<3,double>,double[3][3]> > Particles(0, box, bc, ghost);

		auto it = Particles.getGridIterator(sz);
		while (it.isNext())
		{
			Particles.add();
			auto key = it.get();
			double x = key.get(0) * it.getSpacing(0);
			Particles.getLastPos()[0] = x;
			double y = key.get(1) * it.getSpacing(1);
			Particles.getLastPos()[1] = y;
			double z = key.get(2) * it.getSpacing(1);
			Particles.getLastPos()[2] = z;
			++it;
		}

		Particles.map();
		Particles.ghost_get<0>();

		constexpr int x = 0;
		constexpr int y = 1;
		constexpr int z = 2;


		constexpr int sScalar = 0;
		constexpr int sVector = 1;
		constexpr int sTensor = 2;
		constexpr int dScalar = 0;
		constexpr int dVector = 1;
		constexpr int dTensor = 2;
		auto Pos = getV<PROP_POS>(Particles);
		auto sS = getV<sScalar>(Particles);
		auto sV = getV<sVector>(Particles);
		auto sT = getV<sTensor>(Particles);
		auto dS = getV<dScalar>(Particles);
		auto dV = getV<dVector>(Particles);
		auto dT = getV<dTensor>(Particles);

		//Particles_subset.write("Pars");
		Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dx(Particles, ord,
																								 rCut, sampling_factor,
																								 support_options::RADIUS);
        texp_v<double> TVx,TdxVx;
		texp_v<VectorS<3, double>> TV;
        texp_v<double[3][3]> TT;


		auto it3 = Particles.getDomainIterator();

		sS = 5;
		sV[0] = 1;
		sV[1] = 2;
		sV[2] = 3;
		sT[0][0] = 1;
		sT[0][1] = 2;
		sT[0][2] = 3;
		sT[1][0] = 4;
		sT[1][1] = 5;
		sT[1][2] = 6;
		sT[2][0] = 7;
		sT[2][1] = 8;
		sT[2][2] = 9;

        TVx=sS;
        dS = TVx;

        {
		auto it3 = Particles.getDomainIterator();

		bool match = true;
		while (it3.isNext())
		{
			auto key = it3.get();

			match &= Particles.template getProp<sScalar>(key) == Particles.template getProp<dScalar>(key);

			++it3;
		}

		BOOST_REQUIRE_EQUAL(match,true);
        }

		TVx=sV[0];
		dS = TVx;

        {
		auto it3 = Particles.getDomainIterator();

		bool match = true;
		while (it3.isNext())
		{
			auto key = it3.get();

			match &= Particles.template getProp<sVector>(key)[0] == Particles.template getProp<dScalar>(key);

			++it3;
		}
		BOOST_REQUIRE_EQUAL(match,true);
        }

        TVx=sT[0][0];
        dS = TVx;
        {
		auto it3 = Particles.getDomainIterator();

		bool match = true;
		while (it3.isNext())
		{
			auto key = it3.get();

			match &= Particles.template getProp<sVector>(key)[0] == Particles.template getProp<dScalar>(key);

			++it3;
		}
		BOOST_REQUIRE_EQUAL(match,true);
        }

        TV=sV;
        //TT=sT;
        TdxVx=Dx(sV[x]);
		TT[0][0] = Dx(sV[x]);
    }

BOOST_AUTO_TEST_SUITE_END()

#endif /* DCPSE_OP_TEST_TEMPORAL_CPP_ */
