/*
 * vector_dist_operators_tests_util.hpp
 *
 *  Created on: May 31, 2019
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_OPERATORS_TESTS_UTIL_HPP_
#define VECTOR_DIST_OPERATORS_TESTS_UTIL_HPP_

#include "Operators/Vector/vector_dist_operators.hpp"
#include "Space/Shape/Point.hpp"

constexpr int A = 0;
constexpr int B = 1;
constexpr int C = 2;

constexpr int VA = 3;
constexpr int VB = 4;
constexpr int VC = 5;

constexpr int TA = 6;

//////////////////// Here we define all the function to checl the operators


template <unsigned int prp, unsigned int impl, typename vector>
bool check_values(vector & v,float a)
{
	if (impl == comp_dev)
	{v.template deviceToHostProp<prp>();}

	bool ret = true;
	auto it = v.getDomainIterator();

	while (it.isNext())
	{
		ret &= v.template getProp<prp>(it.get()) == a;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector>
bool check_values_complex_expr(vector & vd)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<B>(key) + 2.0 + vd.template getProp<B>(key) - 2.0*vd.template getProp<C>(key) / 5.0;
		float base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_pos_sum(vector & vd, const rtype & p)
{
	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<VA>();
		vd.deviceToHostPos();
	}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		Point<vector::dims,typename vector::stype> xp = vd.getPos(key);

		rtype base1 = rtype(xp) + p;
		rtype base2 = vd.template getProp<VA>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_pos_integration(vector & vd, double dt)
{
	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<VA,VB>();
		vd.deviceToHostPos();
	}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		Point<vector::dims,typename vector::stype> xp = vd.template getProp<VB>(key);

		rtype base1 = rtype(xp) + rtype(vd.template getProp<VA>(key)) * dt;
		rtype base2 = vd.getPos(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_pos_sub(vector & vd, const rtype & p)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		Point<vector::dims,typename vector::stype> xp = vd.getPos(key);

		rtype base1 = rtype(xp) - p;
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_pos_sub_minus(vector & vd, const rtype & p)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		Point<vector::dims,typename vector::stype> xp = vd.getPos(key);

		rtype base1 = -(rtype(xp) - p);
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_point_sub(vector & vd, const rtype & p)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = -vd.template getProp<B>(key);
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sum(vector & vd, double d)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd.template getProp<B>(key) + d;
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sum(vector & vd1, vector & vd2)
{
	if (impl == comp_dev)
	{
		vd1.template deviceToHostProp<A,B>();
		vd2.template deviceToHostProp<C>();
	}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<B>(key) + vd2.template getProp<C>(key);
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sum_3(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<B>(key) + vd1.template getProp<C>(key) + vd1.template getProp<B>(key);
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sum_4(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = (vd1.template getProp<B>(key) + vd1.template getProp<C>(key)) + (vd1.template getProp<B>(key) + vd1.template getProp<C>(key));
		rtype base2 = vd1.template getProp<A>(key);


		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sub(vector & vd, double d)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd.template getProp<B>(key) - d;
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sub(double d, vector & vd)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = d - vd.template getProp<B>(key);
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sub(vector & vd1, vector & vd2)
{
	if (impl == comp_dev)
	{
		vd1.template deviceToHostProp<A,C>();
		vd2.template deviceToHostProp<B>();
	}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<C>(key) - vd2.template getProp<B>(key);
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sub_31(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<B>(key) - (vd1.template getProp<C>(key) + vd1.template getProp<B>(key));
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sub_32(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = (vd1.template getProp<C>(key) + vd1.template getProp<B>(key)) - vd1.template getProp<B>(key);
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_sub_4(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = (vd1.template getProp<C>(key) + vd1.template getProp<B>(key)) - (vd1.template getProp<C>(key) + vd1.template getProp<B>(key));
		rtype base2 = vd1.template getProp<A>(key);

		if (impl == comp_dev)
		{
			ret &= (double)norm(base1 - base2) < 0.00001;
		}
		else
		{
			ret &= base1 == base2;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_mul(vector & vd, double d)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd.template getProp<B>(key) * d;
		rtype base2 = vd.template getProp<A>(key);

			ret &= base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_mul(vector & vd1, vector & vd2)
{
	if (impl == comp_dev)
	{
		vd1.template deviceToHostProp<A,C>();
		vd2.template deviceToHostProp<B>();
	}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<C>(key) * vd2.template getProp<B>(key);
		rtype base2 = vd1.template getProp<A>(key);

		if (impl == comp_dev)
		{
			ret &= (double)norm(base1 - base2) < 0.00001;
		}
		else
		{
			ret &= base1 == base2;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_mul_3(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<B>(key) * (vd1.template getProp<B>(key) + vd1.template getProp<C>(key));
		rtype base2 = vd1.template getProp<A>(key);

		if (impl == comp_dev)
		{
			ret &= (double)norm(base1 - base2) < 0.00001;
		}
		else
		{
			ret &= base1 == base2;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_mul_4(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = (vd1.template getProp<B>(key) + vd1.template getProp<C>(key)) * (vd1.template getProp<B>(key) + vd1.template getProp<C>(key));
		rtype base2 = vd1.template getProp<A>(key);

		if (impl == comp_dev)
		{
			ret &= (double)norm(base1 - base2) < 0.00001;
		}
		else
		{
			ret &= base1 == base2;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}



template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_div(vector & vd, double d)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd.template getProp<B>(key) / d;
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_div(double d, vector & vd)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,B>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = d / vd.template getProp<B>(key);
		rtype base2 = vd.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_div(vector & vd1, vector & vd2)
{
	if (impl == comp_dev)
	{
		vd1.template deviceToHostProp<A,C>();
		vd2.template deviceToHostProp<B>();
	}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<C>(key) / vd2.template getProp<B>(key);
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template<unsigned int impl, unsigned int prp, typename vector>
bool check_values_pos_exp_slicer(vector & v)
{
	if (impl == comp_dev)
	{
		v.template deviceToHostProp<prp>();
	}

	bool ret = true;
	auto it = v.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		typename vector::stype base1 = -v.getPos(key)[1]*exp(-10.0*(v.getPos(key)[0]*v.getPos(key)[0]+v.getPos(key)[1]*v.getPos(key)[1]));
		typename vector::stype base2 = v.template getProp<prp>(key)[0];

		ret &= fabs(base1 - base2) < 1e-5;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_div_31(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = vd1.template getProp<B>(key) / (vd1.template getProp<B>(key) + vd1.template getProp<C>(key));
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_div_32(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = (vd1.template getProp<C>(key) + vd1.template getProp<B>(key)) / vd1.template getProp<B>(key);
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}


template <typename rtype, typename vector, unsigned int A, unsigned int B, unsigned int C, unsigned int impl>
bool check_values_div_4(vector & vd1)
{
	if (impl == comp_dev)
	{vd1.template deviceToHostProp<A,B,C>();}

	bool ret = true;
	auto it = vd1.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		rtype base1 = (vd1.template getProp<B>(key) + vd1.template getProp<C>(key)) / (vd1.template getProp<B>(key) + vd1.template getProp<C>(key));
		rtype base2 = vd1.template getProp<A>(key);

		ret &=  base1 == base2;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector>
bool check_values_scal_norm_dist(vector & vd)
{
	if (impl == comp_dev)
	{vd.template deviceToHostProp<A,VB,VC>();}

	bool ret = true;
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		float base1 = vd.template getProp<VB>(key) * vd.template getProp<VC>(key) + norm(vd.template getProp<VC>(key) + vd.template getProp<VB>(key)) + distance(vd.template getProp<VC>(key),vd.template getProp<VB>(key));
		float base2 = vd.template getProp<A>(key);

		if (impl == comp_dev)
		{
			ret &= (double)norm(base1 - base2) < 0.00001;
		}
		else
		{
			ret &= base1 == base2;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl,typename vector>
void fill_values(vector & v)
{
	auto it = v.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		v.getPos(p)[0] = (float)rand() / (float)RAND_MAX;
		v.getPos(p)[1] = (float)rand() / (float)RAND_MAX;
		v.getPos(p)[2] = (float)rand() / (float)RAND_MAX;

		v.template getProp<A>(p) = fabs(sin(p.getKey()+1.0));
		v.template getProp<B>(p) = fabs(sin(2.0*p.getKey()+3.0));
		v.template getProp<C>(p) = fabs(sin(3.0*p.getKey()+18.0));

		for (size_t k = 0 ; k < 3 ; k++)
		{
			v.template getProp<VA>(p)[k] = fabs(sin(p.getKey()+1.0+k));
			v.template getProp<VB>(p)[k] = fabs(sin(2.0*p.getKey()+1.0+3.0));
			v.template getProp<VC>(p)[k] = fabs(sin(3.0*p.getKey()+1.0+k));
		}

		++it;
	}

	if (impl == comp_dev)
	{
		v.template hostToDeviceProp<A,B,C,VA,VB,VC>();
		v.hostToDevicePos();
	}
}

template<unsigned int impl,typename vector_type,
		 typename vA_type,
		 typename vB_type,
		 typename vC_type,
		 typename vVA_type,
		 typename vVB_type,
		 typename vVC_type,
		 typename vPOS_type>
void check_all_expressions_imp(vector_type & vd,
						   vA_type vA,
						   vB_type vB,
						   vC_type vC,
						   vVA_type vVA,
						   vVB_type vVB,
						   vVC_type vVC,
						   vPOS_type vPOS)
{
	// vector type
	typedef vector_type vtype;

	vA = 1.0;
	vB = 2.0f;
	vC = 3.0;

	check_values<A,impl>(vd,1.0);
	check_values<B,impl>(vd,2.0);
	check_values<C,impl>(vd,3.0);

	vA = vB;
	check_values<A,impl>(vd,2.0);

	fill_values<impl>(vd);

	vA = vB + 2.0 + vB - 2.0*vC / 5.0;
	check_values_complex_expr<impl>(vd);

	// Various combination of 2 operator

	vA = vB + 2.0;
	check_values_sum<float,vtype,A,B,C,impl>(vd,2.0);
	vA = 2.0 + vB;
	check_values_sum<float,vtype,A,B,C,impl>(vd,2.0);
	vA = vC + vB;
	check_values_sum<float,vtype,A,B,C,impl>(vd,vd);

	vA = vB - 2.0;
	check_values_sub<float,vtype,A,B,C,impl>(vd,2.0);
	vA = 2.0 - vB;
	check_values_sub<float,vtype,A,B,C,impl>(2.0,vd);
	vA = vC - vB;
	check_values_sub<float,vtype,A,B,C,impl>(vd,vd);

	vA = vB * 2.0;
	check_values_mul<float,vtype,A,B,C,impl>(vd,2.0);
	vA = 2.0 * vB;
	check_values_mul<float,vtype,A,B,C,impl>(vd,2.0);
	vA = vC * vB;
	check_values_mul<float,vtype,A,B,C,impl>(vd,vd);

	vA = vB / 2.0;
	check_values_div<float,vtype,A,B,C,impl>(vd,2.0);
	vA = 2.0 / vB;
	check_values_div<float,vtype,A,B,C,impl>(2.0,vd);
	vA = vC / vB;
	check_values_div<float,vtype,A,B,C,impl>(vd,vd);

	// Variuos combination 3 operator

	vA = vB + (vC + vB);
	check_values_sum_3<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) + vB;
	check_values_sum_3<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) + (vC + vB);
	check_values_sum_4<float,vtype,A,B,C,impl>(vd);

	vA = vB - (vC + vB);
	check_values_sub_31<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) - vB;
	check_values_sub_32<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) - (vC + vB);
	check_values_sub_4<float,vtype,A,B,C,impl>(vd);

	vA = vB * (vC + vB);
	check_values_mul_3<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) * vB;
	check_values_mul_3<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) * (vC + vB);
	check_values_mul_4<float,vtype,A,B,C,impl>(vd);

	vA = vB / (vC + vB);
	check_values_div_31<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) / vB;
	check_values_div_32<float,vtype,A,B,C,impl>(vd);
	vA = (vC + vB) / (vC + vB);
	check_values_div_4<float,vtype,A,B,C,impl>(vd);

	if (impl == comp_host)
	{
		auto test = vC + vB;
		auto & v = test.getVector();
		BOOST_REQUIRE_EQUAL((void *)&v,(void *)&vd);
	}

	// We try with vectors

	// Various combination of 2 operator

	vVA = vVB + 2.0;
	check_values_sum<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,2.0f);
	vVA = 2.0 + vVB;
	check_values_sum<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,2.0f);
	vVA = vVC + vVB;
	check_values_sum<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,vd);

	vVA = vVB - 2.0;
	check_values_sub<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,2.0f);
	vVA = 2.0 - vVB;
	check_values_sub<VectorS<3,float>,vtype,VA,VB,VC,impl>(2.0f,vd);
	vVA = vVC - vVB;
	check_values_sub<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,vd);

	vVA = vVB * 2.0;
	check_values_mul<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,2.0f);
	vVA = 2.0 * vVB;
	check_values_mul<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,2.0f);
	vVA = vVC * vVB;
	check_values_mul<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,vd);

	vVA = vVB / 2.0;
	check_values_div<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,2.0f);
	vVA = 2.0 / vVB;
	check_values_div<VectorS<3,float>,vtype,VA,VB,VC,impl>(2.0f,vd);
	vVA = vVC / vVB;
	check_values_div<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,vd);

	if (impl == comp_host)
	{
		auto test = vVB / 2.0;
		auto & v = test.getVector();
		BOOST_REQUIRE_EQUAL((void *)&v,(void *)&vd);
	}

	// Variuos combination 3 operator

	vVA = vVB + (vVC + vVB);
	check_values_sum_3<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) + vVB;
	check_values_sum_3<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) + (vVC + vVB);
	check_values_sum_4<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);

	vVA = vVB - (vVC + vVB);
	check_values_sub_31<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) - vVB;
	check_values_sub_32<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) - (vVC + vVB);
	check_values_sub_4<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);

	vVA = vVB * (vVC + vVB);
	check_values_mul_3<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) * vVB;
	check_values_mul_3<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) * (vVC + vVB);
	check_values_mul_4<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vA = vVB * (vVC + vVB);
	check_values_mul_3<float,vtype,A,VB,VC,impl>(vd);
	vA = (vVC + vVB) * vVB;
	check_values_mul_3<float,vtype,A,VB,VC,impl>(vd);
	vA = (vVC + vVB) * (vVC + vVB);
	check_values_mul_4<float,vtype,A,VB,VC,impl>(vd);

	if (impl == comp_host)
	{
		auto test = (vVC + vVB) * (vVC + vVB);
		auto & v = test.getVector();
		BOOST_REQUIRE_EQUAL((void *)&v,(void *)&vd);
	}

	vVA = vVB / (vVC + vVB);
	check_values_div_31<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) / vVB;
	check_values_div_32<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);
	vVA = (vVC + vVB) / (vVC + vVB);
	check_values_div_4<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd);

	// normalization function

	vA = vVB * vVC + norm(vVC + vVB) + openfpm::distance(vVC,vVB);
	check_values_scal_norm_dist<impl>(vd);

	Point<3,float> p0({2.0,2.0,2.0});
	auto p0_e = getVExpr(p0);

	vVA = vPOS + p0_e;

	check_values_pos_sum<VectorS<3,float>,vtype,VA,VB,VC,impl>(vd,p0);

	vVA = vPOS - p0_e;
	check_values_pos_sub<Point<3,float>,vtype,VA,VB,VC,impl>(vd,p0);

	vVA = -(vPOS - p0_e);
	check_values_pos_sub_minus<Point<3,float>,vtype,VA,VB,VC,impl>(vd,p0);

	vVA = -vVB;
	check_values_point_sub<Point<3,float>,vtype,VA,VB,VC,impl>(vd,p0);

	if (impl == comp_host)
	{
		auto test = vPOS + p0_e;
		auto & v = test.getVector();
		BOOST_REQUIRE_EQUAL((void *)&v,(void *)&vd);
	}

	// Just check it compile testing it will test the same code
	// as the previuous one
	vVC = exp(vVB);
	vA = norm(vPOS);
	vVA = vPOS + 2.0;
	vVA = 2.0 + vPOS;
	vVA = vPOS + vPOS;
	vVA = vPOS - 2.0f;
	vVA = 2.0 - vPOS;
	vVA = vPOS - vPOS;

	if (impl == comp_host)
	{
		auto test = exp(vVB);
		auto & v = test.getVector();
		BOOST_REQUIRE_EQUAL((void *)&v,(void *)&vd);
	}

	vVA = vPOS * 2.0;
	vVA = 2.0 * vPOS;
	vVA = vPOS * vPOS;

	vVA = vPOS / 2.0f;
	vVA = 2.0f / vPOS;
	vVA = vPOS / vPOS;

	// Variuos combination 3 operator

	vVA = vPOS + (vPOS + vPOS);
	vVA = (vPOS + vPOS) + vPOS;
	vVA = (vPOS + vPOS) + (vPOS + vPOS);

	vVA = vPOS - (vPOS + vPOS);
	vVA = (vPOS + vPOS) - vPOS;
	vVA = (vVC + vPOS) - (vPOS + vPOS);

	vVA = vPOS * (vPOS + vPOS);
	vVA = (vPOS + vPOS) * vPOS;
	vVA = (vPOS + vPOS) * (vPOS + vPOS);
	vA = vPOS * (vPOS + vPOS);
	vA = (vPOS + vPOS) * vPOS;
	vA = (vPOS + vPOS) * (vPOS + vPOS);

	vVA = vPOS / (vPOS + vPOS);
	vVA = (vPOS + vPOS) / vPOS;
	vVA = (vPOS + vPOS) / (vPOS + vPOS);

	vVB = vPOS;
	double dt = 0.1;
	vPOS = vPOS + vVA * dt;

	check_values_pos_integration<Point<3,float>,vtype,VA,VB,VC,impl>(vd,dt);

	// Position with slicer (not tested on GPU)

	vVA[0]=-vPOS[1]*exp(-10.0*(vPOS[0]*vPOS[0]+vPOS[1]*vPOS[1]));
	check_values_pos_exp_slicer<impl,VA>(vd);
}

template<unsigned int impl>
struct check_all_expressions
{
	template<typename vector_type> static void check(vector_type & vd)
	{
		auto vA = getV<A>(vd);
		auto vB = getV<B>(vd);
		auto vC = getV<C>(vd);

		auto vVA = getV<VA>(vd);
		auto vVB = getV<VB>(vd);
		auto vVC = getV<VC>(vd);

		auto vPOS = getV<PROP_POS>(vd);

		check_all_expressions_imp<impl>(vd,vA,vB,vC,vVA,vVB,vVC,vPOS);
	}
};

template<>
struct check_all_expressions<comp_dev>
{
	template<typename vector_type> static void check(vector_type & vd)
	{
		auto vdk = vd.toKernel();

		auto vA = getV<A>(vdk);
		auto vB = getV<B>(vdk);
		auto vC = getV<C>(vdk);

		auto vVA = getV<VA>(vdk);
		auto vVB = getV<VB>(vdk);
		auto vVC = getV<VC>(vdk);

		auto vPOS = getV<PROP_POS>(vdk);

		check_all_expressions_imp<comp_dev>(vd,vA,vB,vC,vVA,vVB,vVC,vPOS);
	}
};


template <unsigned impl, typename vector,typename Kernel, typename NN_type>
bool check_values_apply_kernel(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;

	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<A,C,VB,VC>();
		vd.deviceToHostPos();
	}

	auto it = vd.getDomainIterator();

	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		float base2 = 0.0;
		auto p = it2.get();

		Point<3,float> xp = vd.getPos(p);

		float base1 = vd.template getProp<A>(p);
		float prp_x = vd.template getProp<VC>(p) * vd.template getProp<VB>(p) + norm(vd.template getProp<VB>(p));

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

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

		if (impl == comp_host)
		{ret &= base1 == base2;}
		else
		{ret &= fabs(base1 - base2) < 0.0001;}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector,typename Kernel, typename NN_type>
bool check_values_apply_kernel_reduce(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;

	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<A,C,VB,VC>();
		vd.deviceToHostPos();
	}

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
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

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

		if (impl == comp_host)
		{ret &= base1 == base3;}
		else
		{ret &= fabs(base1 - base3) < 0.001;}

		++it3;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector,typename Kernel, typename NN_type>
bool check_values_apply_kernel2(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;

	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<VA,VB,VC>();
		vd.deviceToHostPos();
	}

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
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

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

		if (impl == comp_host)
		{ret &= base1 == base2;}
		else
		{
			for (size_t i = 0 ; i < 3 ; i++)
			{ret &= fabs(base1.get(i) - base2.get(i)) < 0.0001;}
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector,typename Kernel, typename NN_type>
bool check_values_apply_kernel3(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;

	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<VA,VC>();
		vd.deviceToHostPos();
	}

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

		Point<3,float> prp_x = vd.template getProp<VC>(p);

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// position q
			Point<3,float> xq = vd.getPos(q);

			Point<3,float> prp_y = vd.template getProp<VC>(q);

			base2 += ker.value(xp,xq,prp_x,prp_y);

			++Np;
		}

		base2 += vd.template getProp<VC>(p);

		if (impl == comp_host)
		{ret &= base1 == base2;}
		else
		{
			for (size_t i = 0 ; i < 3 ; i++)
			{ret &= fabs(base1.get(i) - base2.get(i)) < 0.0001;}
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector,typename Kernel, typename NN_type>
bool check_values_apply_kernel2_reduce(vector & vd, Kernel & ker, NN_type & NN)
{
	bool ret = true;

	if (impl == comp_dev)
	{
		vd.template deviceToHostProp<VA,VB,VC>();
		vd.deviceToHostPos();
	}

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
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

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

		if (impl == comp_host)
		{ret &= base1 == base3;}
		else
		{
			for (size_t i = 0 ; i < 3 ; i++)
			{ret &= fabs(base1.get(i) - base3.get(i)) < 0.002;}

			if (ret == false)
			{
				int debug = 0;
				debug++;
			}
		}

		++it3;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	return ret;
}

template <unsigned int impl, typename vector,typename Kernel, typename NN_type>
bool check_values_apply_kernel3_reduce(vector & vd, Kernel & ker, NN_type & NN, const Point<2,float> & p)
{
	bool ret = true;

	if (impl == comp_dev)
	{
		vd.deviceToHostPos();
	}

	auto it = vd.getDomainIterator();

	Point<2,float> base2 = 0.0;

	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,float> xp = vd.getPos(p);

		Point<2,float> ker_accu = 0.0;

		// For each neighborhood particle
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(xp));

		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			if (q == p.getKey())	{++Np; continue;};

			// position q
			Point<3,float> xq = vd.getPos(q);

			ker_accu += ker.value(xp,xq);

			++Np;
		}

		base2 += ker_accu;

		++it2;
	}

	if (impl == comp_host)
	{
		BOOST_REQUIRE_EQUAL(p.get(0),base2.get(0));
		BOOST_REQUIRE_EQUAL(p.get(1),base2.get(1));
	}
	else
	{
		BOOST_REQUIRE(fabs(p.get(0) - base2.get(0)) < 0.001);
		BOOST_REQUIRE(fabs(p.get(1) - base2.get(1)) < 0.001);
	}

	return ret;
}

typedef vector_dist<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>,float>> vector_type;

#ifdef CUDA_GPU
typedef vector_dist_ker<3,float,aggregate<float,float,float,VectorS<3,float>,VectorS<3,float>,VectorS<3,float>,float>> vector_type_ker;
#endif

//! Exponential kernel
struct exp_kernel
{
	//! variance of the exponential kernel
	float var;

	//! Exponential kernel giving variance
	exp_kernel(float var)
	:var(var)
	{}

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 *
	 * \return the result
	 *
	 */
	__device__ __host__ inline float value(const Point<3,float> & p, const Point<3,float> & q,float pA,float pB)
	{
		float dist = norm(p-q);

		return (pA + pB) * exp(dist * dist / var);
	}

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 *
	 * \return the result
	 *
	 */
	__device__ __host__ inline Point<3,float> value(const Point<3,float> & p, const Point<3,float> & q,const Point<3,float> & pA, const Point<3,float> & pB)
	{
		float dist = norm(p-q);

		return (pA + pB) * exp(dist * dist / var);
	}

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 * \param vd1 original vector
	 *
	 * \return the result
	 *
	 */
	template<typename vector_t>
	__host__ __device__ inline float value(size_t p, size_t q, float pA, float pB, const vector_t & vd1)
	{
		Point<3,float> pp = vd1.getPos(p);
		Point<3,float> pq = vd1.getPos(q);

		float dist = norm(pp-pq);

		return (pA + pB) * exp(dist * dist / var);
	}

#ifdef CUDA_GPU

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 * \param vd1 original vector
	 *
	 * \return the result
	 *
	 */
	__device__ inline float value(size_t p, size_t q, float pA, float pB, const vector_type_ker & vd1)
	{
		Point<3,float> pp = vd1.getPos(p);
		Point<3,float> pq = vd1.getPos(q);

		float dist = norm(pp-pq);

		return (pA + pB) * exp(dist * dist / var);
	}

#endif

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 * \param vd1 original vector
	 *
	 * \return the result
	 *
	 */
	__host__ inline Point<3,float> value(size_t p, size_t q, const Point<3,float> & pA, const Point<3,float> & pB , const vector_type & vd1)
	{
		Point<3,float> pp = vd1.getPos(p);
		Point<3,float> pq = vd1.getPos(q);

		float dist = norm(pp-pq);

		return (pA + pB) * exp(dist * dist / var);
	}

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 * \param vd1 original vector
	 *
	 * \return the result
	 *
	 */
	template<typename vector_t>
	__host__ inline Point<3,float> value(size_t p, size_t q, const Point<3,float> & pA, const Point<3,float> & pB , const vector_t & vd1)
	{
		Point<3,float> pp = vd1.getPos(p);
		Point<3,float> pq = vd1.getPos(q);

		float dist = norm(pp-pq);

		return (pA + pB) * exp(dist * dist / var);
	}

#ifdef CUDA_GPU

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 * \param pA property value at p
	 * \param pB property value at q
	 * \param vd1 original vector
	 *
	 * \return the result
	 *
	 */
	__device__ inline Point<3,float> value(size_t p, size_t q, const Point<3,float> & pA, const Point<3,float> & pB , const vector_type_ker & vd1)
	{
		Point<3,float> pp = vd1.getPos(p);
		Point<3,float> pq = vd1.getPos(q);

		float dist = norm(pp-pq);

		return (pA + pB) * exp(dist * dist / var);
	}

#endif

	/*! \brief Result of the exponential kernel
	 *
	 * \param p position of the particle p
	 * \param q position of the particle q
	 *
	 * \return the result
	 *
	 */
	__device__ __host__ inline Point<2,float> value(const Point<3,float> & p, const Point<3,float> & q)
	{
		float dist = norm(p-q);

		return exp(dist * dist / var);
	}
};

template<unsigned int impl,typename vector,
		 typename vA_type,
		 typename vC_type,
		 typename vVA_type,
		 typename vVB_type,
		 typename vVC_type>
void vector_dist_op_ap_ker_impl(vector & vd, vA_type & vA,
											 vC_type & vC,
											 vVA_type & vVA,
											 vVB_type & vVB,
											 vVC_type & vVC,
											 unsigned int opt)
{
	// we apply an exponential kernel to calculate something

	auto cl = vd.template getCellListDev<impl>(0.05);
	auto cl_host = vd.template getCellListDev<comp_host>(0.05);
	exp_kernel ker(0.2);

	vA = applyKernel_in(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	check_values_apply_kernel<impl>(vd,ker,cl_host);

	vVA = applyKernel_in(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	check_values_apply_kernel2<impl>(vd,ker,cl_host);

	vA = rsum(applyKernel_in(vVC * vVB + norm(vVB),vd,cl,ker)) + vC;
	check_values_apply_kernel_reduce<impl>(vd,ker,cl_host);

	vVA = rsum(applyKernel_in(2.0*vVC + vVB ,vd,cl,ker)) + vVC;
	check_values_apply_kernel2_reduce<impl>(vd,ker,cl_host);

	vA = applyKernel_in_gen(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	check_values_apply_kernel<impl>(vd,ker,cl_host);

	vVA = applyKernel_in_gen(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	check_values_apply_kernel2<impl>(vd,ker,cl_host);

	vA = rsum(applyKernel_in_gen(vVC * vVB + norm(vVB),vd,cl,ker)) + vC;
	check_values_apply_kernel_reduce<impl>(vd,ker,cl_host);

	vVA = rsum(applyKernel_in_gen(2.0*vVC + vVB ,vd,cl,ker)) + vVC;
	check_values_apply_kernel2_reduce<impl>(vd,ker,cl_host);

	// Check it compile the code is the same
	vVA = applyKernel_in_gen(vVC,vd,cl,ker) + vVC;
	check_values_apply_kernel3<impl>(vd,ker,cl_host);

	vVA = applyKernel_in    (vVC,vd,cl,ker) + vVC;
	check_values_apply_kernel3<impl>(vd,ker,cl_host);

	Point<2,float> p = rsum(applyKernel_in_sim(vd,cl,ker)).get();
	check_values_apply_kernel3_reduce<impl>(vd,ker,cl_host,p);
}

template<typename vector,
		 typename vA_type,
		 typename vC_type,
		 typename vVA_type,
		 typename vVB_type,
		 typename vVC_type>
void vector_dist_op_ap_ker_impl_sort(vector & vd, vA_type & vA,
											 vC_type & vC,
											 vVA_type & vVA,
											 vVB_type & vVB,
											 vVC_type & vVC,
											 unsigned int opt)
{
	// we apply an exponential kernel to calculate something

	auto cl_gpu = vd.getCellListGPU(0.05);
	auto cl = cl_gpu.toKernel();
	auto cl_host = vd.template getCellListDev<comp_host>(0.05);
	exp_kernel ker(0.2);

	vA = applyKernel_in_sort(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	vd.template merge_sort<A>(cl_gpu);
	check_values_apply_kernel<comp_dev>(vd,ker,cl_host);

	vVA = applyKernel_in_sort(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	vd.template merge_sort<VA>(cl_gpu);
	check_values_apply_kernel2<comp_dev>(vd,ker,cl_host);

/*	vA = rsum(applyKernel_in_sort(vVC * vVB + norm(vVB),vd,cl,ker)) + vC;
	vd.template merge_sort<A>(cl_gpu);
	check_values_apply_kernel_reduce<comp_dev>(vd,ker,cl_host);

	vVA = rsum(applyKernel_in_sort(2.0*vVC + vVB ,vd,cl,ker)) + vVC;
	vd.template merge_sort<VA>(cl_gpu);
	check_values_apply_kernel2_reduce<comp_dev>(vd,ker,cl_host);*/

	vA = applyKernel_in_gen_sort(vVC * vVB + norm(vVB),vd,cl,ker) + vC;
	vd.template merge_sort<A>(cl_gpu);
	check_values_apply_kernel<comp_dev>(vd,ker,cl_host);

	vVA = applyKernel_in_gen_sort(2.0*vVC + vVB ,vd,cl,ker) + vVC;
	vd.template merge_sort<VA>(cl_gpu);
	check_values_apply_kernel2<comp_dev>(vd,ker,cl_host);

/*	vA = rsum(applyKernel_in_gen_sort(vVC * vVB + norm(vVB),vd,cl,ker)) + vC;
	vd.template merge_sort<A>(cl_gpu);
	check_values_apply_kernel_reduce<comp_dev>(vd,ker,cl_host);

	vVA = rsum(applyKernel_in_gen_sort(2.0*vVC + vVB ,vd,cl,ker)) + vVC;
	vd.template merge_sort<VA>(cl_gpu);
	check_values_apply_kernel2_reduce<comp_dev>(vd,ker,cl_host);*/

	// Check it compile the code is the same
	vVA = applyKernel_in_gen_sort(vVC,vd,cl,ker) + vVC;
	vd.template merge_sort<VA>(cl_gpu);
	check_values_apply_kernel3<comp_dev>(vd,ker,cl_host);

	vVA = applyKernel_in_sort(vVC,vd,cl,ker) + vVC;
	vd.template merge_sort<VA>(cl_gpu);
	check_values_apply_kernel3<comp_dev>(vd,ker,cl_host);

/*	Point<2,float> p = rsum(applyKernel_in_sim_sort(vd,cl,ker)).get();
	check_values_apply_kernel3_reduce<comp_dev>(vd,ker,cl_host,p);*/
}

template<unsigned int impl>
struct check_all_apply_ker
{
	template<typename vector_type> static void check(vector_type & vd)
	{
		// fill vd with some value
		fill_values<impl>(vd);

		vd.map();
		vd.template ghost_get<0,1,2,3,4,5,6>();

		auto vA = getV<A>(vd);
		auto vC = getV<C>(vd);

		auto vVA = getV<VA>(vd);
		auto vVB = getV<VB>(vd);
		auto vVC = getV<VC>(vd);

		vector_dist_op_ap_ker_impl<impl>(vd,vA,vC,vVA,vVB,vVC,NONE);
	}
};


template<>
struct check_all_apply_ker<comp_dev>
{
	template<typename vector_type> static void check(vector_type & vd)
	{
		auto vA = getV<A,comp_dev>(vd);
		auto vC = getV<C,comp_dev>(vd);

		auto vVA = getV<VA,comp_dev>(vd);
		auto vVB = getV<VB,comp_dev>(vd);
		auto vVC = getV<VC,comp_dev>(vd);

		// fill vd with some value
		fill_values<comp_dev>(vd);

		vd.map(RUN_ON_DEVICE);
		vd.template ghost_get<0,1,2,3,4,5,6>(RUN_ON_DEVICE);
		vd.template deviceToHostProp<0,1,2,3,4,5,6>();
		vd.deviceToHostPos();

		vector_dist_op_ap_ker_impl<comp_dev>(vd,vA,vC,vVA,vVB,vVC,RUN_ON_DEVICE);
	}
};


struct check_all_apply_ker_sort
{
	template<typename vector_type> static void check(vector_type & vd)
	{
		auto vA = getV_sort<A>(vd);
		auto vC = getV_sort<C>(vd);

		auto vVA = getV_sort<VA>(vd);
		auto vVB = getV_sort<VB>(vd);
		auto vVC = getV_sort<VC>(vd);

		// fill vd with some value
		fill_values<comp_dev>(vd);

		vd.map(RUN_ON_DEVICE);
		vd.template ghost_get<0,1,2,3,4,5,6>(RUN_ON_DEVICE);
		vd.template deviceToHostProp<0,1,2,3,4,5,6>();
		vd.deviceToHostPos();

		vector_dist_op_ap_ker_impl_sort(vd,vA,vC,vVA,vVB,vVC,RUN_ON_DEVICE);
	}
};

#endif /* VECTOR_DIST_OPERATORS_TESTS_UTIL_HPP_ */
