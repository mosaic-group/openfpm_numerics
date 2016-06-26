/*
 * vector_dist_operators_apply_kernel_unit_tests.hpp
 *
 *  Created on: Jun 19, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_APPLY_KERNEL_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_APPLY_KERNEL_UNIT_TESTS_HPP_

template <typename vector,typename Kernel, typename NN_type> bool check_values_apply_kernel(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	// ### WIKI 13 ###
	//
	// Check that apply kernel work
	//
	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		float base2 = 0.0;
		auto p = it2.get();

		Point<3,float> xp = vd.getPos(p);

		float base1 = vd.template getProp<A>(p);

		float prp_x = vd.template getProp<VC>(p) * vd.template getProp<VB>(p) + norm(vd.template getProp<VB>(p));

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// position q
			Point<3,float> xq = vd.getPos(q);

			float prp_y = vd.template getProp<VC>(q) * vd.template getProp<VB>(q) + norm(vd.template getProp<VB>(q));

			base2 += ker.value(xp,xq,prp_x,prp_y);

			++Np;
		}

		base2 += vd.template getProp<C>(p);

		ret &= base1 == base2;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector,typename Kernel, typename NN_type> bool check_values_apply_kernel_reduce(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	float base1 = 0.0;
	float base2 = 0.0;
	float base3 = 0.0;

	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,float> xp = vd.getPos(p);

		float ker_accu = 0.0;
		float prp_x = vd.template getProp<VC>(p) * vd.template getProp<VB>(p) + norm(vd.template getProp<VB>(p));

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// position q
			Point<3,float> xq = vd.getPos(q);

			float prp_y = vd.template getProp<VC>(q) * vd.template getProp<VB>(q) + norm(vd.template getProp<VB>(q));

			ker_accu += ker.value(xp,xq,prp_x,prp_y);

			++Np;
		}

		base2 += ker_accu;

		++it2;
	}

	auto it3 = vd.getDomainIterator();
	while (it3.isNext())
	{
		auto p = it3.get();

		base1 = vd.template getProp<A>(p);
		base3 = vd.template getProp<C>(p) + base2;

		ret &= base1 == base3;

		++it3;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename vector,typename Kernel, typename NN_type> bool check_values_apply_kernel2(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	// ### WIKI 13 ###
	//
	// Check that apply kernel work
	//
	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		Point<3,float> base2 = 0.0;
		auto p = it2.get();

		Point<3,float> xp = vd.getPos(p);

		Point<3,float> base1 = vd.template getProp<VA>(p);

		Point<3,float> prp_x = 2.0 * vd.template getProp<VC>(p) + vd.template getProp<VB>(p);

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// position q
			Point<3,float> xq = vd.getPos(q);

			Point<3,float> prp_y = 2.0 * vd.template getProp<VC>(q) + vd.template getProp<VB>(q);

			base2 += ker.value(xp,xq,prp_x,prp_y);

			++Np;
		}

		base2 += vd.template getProp<VC>(p);

		ret &= base1 == base2;

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename vector,typename Kernel, typename NN_type> bool check_values_apply_kernel2_reduce(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;
	auto it = vd.getDomainIterator();

	Point<3,float> base1 = 0.0;
	Point<3,float> base2 = 0.0;
	Point<3,float> base3 = 0.0;

	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,float> xp = vd.getPos(p);

		Point<3,float> ker_accu = 0.0;
		Point<3,float> prp_x = 2.0f*vd.template getProp<VC>(p) + vd.template getProp<VB>(p);

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// position q
			Point<3,float> xq = vd.getPos(q);
			Point<3,float> prp_y = 2.0f*vd.template getProp<VC>(q) + vd.template getProp<VB>(q);

			ker_accu += ker.value(xp,xq,prp_x,prp_y);

			++Np;
		}

		base2 += ker_accu;

		++it2;
	}

	auto it3 = vd.getDomainIterator();
	while (it3.isNext())
	{
		auto p = it3.get();

		base1 = vd.template getProp<VA>(p);
		base3 = vd.template getProp<VC>(p) + base2;

		ret &= base1 == base3;

		++it3;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

BOOST_AUTO_TEST_CASE( vector_dist_operators_apply_kernel_test )
{
	if (create_vcluster().getProcessingUnits() > 3)
		return;

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost
	Ghost<3,float> ghost(0.05);

	// vector type
	typedef vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>,float>> vtype;

	vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>,float>> vd(512,box,bc,ghost);

	auto vA = getV<A>(vd);
	auto vB = getV<B>(vd);
	auto vC = getV<C>(vd);

	auto vVA = getV<VA>(vd);
	auto vVB = getV<VB>(vd);
	auto vVC = getV<VC>(vd);

	// fill vd with some value
	fill_values(vd);

	vd.map();
	vd.ghost_get<0,1,2,3,4,5,6>();

	// we apply an exponential kernel to calculate something

	auto cl = vd.getCellList(0.05);
	exp_kernel ker(0.2);

	vA = applyKernel_in(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	check_values_apply_kernel(vd,ker,cl);

	vVA = applyKernel_in(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	check_values_apply_kernel2(vd,ker,cl);

	vA = applyKernel_reduce(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	check_values_apply_kernel_reduce(vd,ker,cl);

	vVA = applyKernel_reduce(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	check_values_apply_kernel2_reduce(vd,ker,cl);

/*	vA = applyKernel_in_gen(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	check_values_apply_kernel(vd,ker,cl);

	vVA = applyKernel_in_gen(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	check_values_apply_kernel2(vd,ker,cl);

	vA = applyKernel_reduce_gen(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	check_values_apply_kernel_reduce(vd,ker,cl);

	vVA = applyKernel_reduce_gen(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	check_values_apply_kernel2_reduce(vd,ker,cl);*/
}


#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATORS_APPLY_KERNEL_UNIT_TESTS_HPP_ */
