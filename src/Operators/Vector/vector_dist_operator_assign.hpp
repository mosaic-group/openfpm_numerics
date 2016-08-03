/*
 * vector_dist_operator_assign.ipp
 *
 *  Created on: Jul 21, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATOR_ASSIGN_HPP_
#define OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATOR_ASSIGN_HPP_

template<typename T>
struct construct_expression
{
	static inline const T & construct(const T & e)
	{
		return e;
	}
};

template<>
struct construct_expression<double>
{
	static inline vector_dist_expression<0,double> construct(double e)
	{
		return vector_dist_expression<0,double>(e);
	}
};

template<>
struct construct_expression<float>
{
	static inline vector_dist_expression<0,float> construct(const float e)
	{
		return vector_dist_expression<0,float>(e);
	}
};

template<typename prp1, typename expr1, typename prp2, typename expr2> void assign(prp1 & p1, const expr1 & v_e1, prp2 & p2, const expr2 & v_e2)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);

	v_exp1.init();
	v_exp2.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);

		++it;
	}
}

template<typename prp1, typename expr1, typename prp2, typename expr2, typename prp3, typename expr3> void assign(prp1 & p1, const expr1 & v_e1, prp2 & p2, const expr2 & v_e2, prp3 & p3, const expr3 & v_e3)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);
	auto v_exp3 = construct_expression<expr3>::construct(v_e3);

	v_exp1.init();
	v_exp2.init();
	v_exp3.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);
		pos_or_prop<typename prp3::vtype,prp3::prop>::value(p3.getVector(),key) = v_exp3.value(key);

		++it;
	}
}

template<typename prp1, typename expr1,
		 typename prp2, typename expr2,
		 typename prp3, typename expr3,
		 typename prp4, typename expr4>
void assign(prp1 & p1, const expr1 & v_e1,
			prp2 & p2, const expr2 & v_e2,
			prp3 & p3, const expr3 & v_e3,
			prp4 & p4, const expr4 & v_e4)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);
	auto v_exp3 = construct_expression<expr3>::construct(v_e3);
	auto v_exp4 = construct_expression<expr4>::construct(v_e4);

	v_exp1.init();
	v_exp2.init();
	v_exp3.init();
	v_exp4.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);
		pos_or_prop<typename prp3::vtype,prp3::prop>::value(p3.getVector(),key) = v_exp3.value(key);
		pos_or_prop<typename prp4::vtype,prp4::prop>::value(p4.getVector(),key) = v_exp4.value(key);

		++it;
	}
}


template<typename prp1, typename expr1,
		 typename prp2, typename expr2,
		 typename prp3, typename expr3,
		 typename prp4, typename expr4,
		 typename prp5, typename expr5>
void assign(prp1 & p1, const expr1 & v_e1,
			prp2 & p2, const expr2 & v_e2,
			prp3 & p3, const expr3 & v_e3,
			prp4 & p4, const expr4 & v_e4,
			prp5 & p5, const expr5 & v_e5)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);
	auto v_exp3 = construct_expression<expr3>::construct(v_e3);
	auto v_exp4 = construct_expression<expr4>::construct(v_e4);
	auto v_exp5 = construct_expression<expr5>::construct(v_e5);

	v_exp1.init();
	v_exp2.init();
	v_exp3.init();
	v_exp4.init();
	v_exp5.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);
		pos_or_prop<typename prp3::vtype,prp3::prop>::value(p3.getVector(),key) = v_exp3.value(key);
		pos_or_prop<typename prp4::vtype,prp4::prop>::value(p4.getVector(),key) = v_exp4.value(key);
		pos_or_prop<typename prp5::vtype,prp5::prop>::value(p5.getVector(),key) = v_exp5.value(key);

		++it;
	}
}

template<typename prp1, typename expr1,
		 typename prp2, typename expr2,
		 typename prp3, typename expr3,
		 typename prp4, typename expr4,
		 typename prp5, typename expr5,
		 typename prp6, typename expr6>
void assign(prp1 & p1, const expr1 & v_e1,
			prp2 & p2, const expr2 & v_e2,
			prp3 & p3, const expr3 & v_e3,
			prp4 & p4, const expr4 & v_e4,
			prp5 & p5, const expr5 & v_e5,
			prp6 & p6, const expr6 & v_e6)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);
	auto v_exp3 = construct_expression<expr3>::construct(v_e3);
	auto v_exp4 = construct_expression<expr4>::construct(v_e4);
	auto v_exp5 = construct_expression<expr5>::construct(v_e5);
	auto v_exp6 = construct_expression<expr6>::construct(v_e6);

	v_exp1.init();
	v_exp2.init();
	v_exp3.init();
	v_exp4.init();
	v_exp5.init();
	v_exp6.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);
		pos_or_prop<typename prp3::vtype,prp3::prop>::value(p3.getVector(),key) = v_exp3.value(key);
		pos_or_prop<typename prp4::vtype,prp4::prop>::value(p4.getVector(),key) = v_exp4.value(key);
		pos_or_prop<typename prp5::vtype,prp5::prop>::value(p5.getVector(),key) = v_exp5.value(key);
		pos_or_prop<typename prp6::vtype,prp6::prop>::value(p6.getVector(),key) = v_exp6.value(key);

		++it;
	}
}



template<typename prp1, typename expr1,
		 typename prp2, typename expr2,
		 typename prp3, typename expr3,
		 typename prp4, typename expr4,
		 typename prp5, typename expr5,
		 typename prp6, typename expr6,
		 typename prp7, typename expr7>
void assign(prp1 & p1, const expr1 & v_e1,
			prp2 & p2, const expr2 & v_e2,
			prp3 & p3, const expr3 & v_e3,
			prp4 & p4, const expr4 & v_e4,
			prp5 & p5, const expr5 & v_e5,
			prp6 & p6, const expr6 & v_e6,
			prp7 & p7, const expr7 & v_e7)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);
	auto v_exp3 = construct_expression<expr3>::construct(v_e3);
	auto v_exp4 = construct_expression<expr4>::construct(v_e4);
	auto v_exp5 = construct_expression<expr5>::construct(v_e5);
	auto v_exp6 = construct_expression<expr6>::construct(v_e6);
	auto v_exp7 = construct_expression<expr7>::construct(v_e7);

	v_exp1.init();
	v_exp2.init();
	v_exp3.init();
	v_exp4.init();
	v_exp5.init();
	v_exp6.init();
	v_exp7.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);
		pos_or_prop<typename prp3::vtype,prp3::prop>::value(p3.getVector(),key) = v_exp3.value(key);
		pos_or_prop<typename prp4::vtype,prp4::prop>::value(p4.getVector(),key) = v_exp4.value(key);
		pos_or_prop<typename prp5::vtype,prp5::prop>::value(p5.getVector(),key) = v_exp5.value(key);
		pos_or_prop<typename prp6::vtype,prp6::prop>::value(p6.getVector(),key) = v_exp6.value(key);
		pos_or_prop<typename prp7::vtype,prp7::prop>::value(p7.getVector(),key) = v_exp7.value(key);

		++it;
	}
}


template<typename prp1, typename expr1,
		 typename prp2, typename expr2,
		 typename prp3, typename expr3,
		 typename prp4, typename expr4,
		 typename prp5, typename expr5,
		 typename prp6, typename expr6,
		 typename prp7, typename expr7,
		 typename prp8, typename expr8>
void assign(prp1 & p1, const expr1 & v_e1,
			prp2 & p2, const expr2 & v_e2,
			prp3 & p3, const expr3 & v_e3,
			prp4 & p4, const expr4 & v_e4,
			prp5 & p5, const expr5 & v_e5,
			prp6 & p6, const expr6 & v_e6,
			prp7 & p7, const expr7 & v_e7,
			prp8 & p8, const expr8 & v_e8)
{
	auto v_exp1 = construct_expression<expr1>::construct(v_e1);
	auto v_exp2 = construct_expression<expr2>::construct(v_e2);
	auto v_exp3 = construct_expression<expr3>::construct(v_e3);
	auto v_exp4 = construct_expression<expr4>::construct(v_e4);
	auto v_exp5 = construct_expression<expr5>::construct(v_e5);
	auto v_exp6 = construct_expression<expr6>::construct(v_e6);
	auto v_exp7 = construct_expression<expr7>::construct(v_e7);
	auto v_exp8 = construct_expression<expr8>::construct(v_e8);

	v_exp1.init();
	v_exp2.init();
	v_exp3.init();
	v_exp4.init();
	v_exp5.init();
	v_exp6.init();
	v_exp7.init();
	v_exp8.init();

	auto it = p1.getVector().getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		pos_or_prop<typename prp1::vtype,prp1::prop>::value(p1.getVector(),key) = v_exp1.value(key);
		pos_or_prop<typename prp2::vtype,prp2::prop>::value(p2.getVector(),key) = v_exp2.value(key);
		pos_or_prop<typename prp3::vtype,prp3::prop>::value(p3.getVector(),key) = v_exp3.value(key);
		pos_or_prop<typename prp4::vtype,prp4::prop>::value(p4.getVector(),key) = v_exp4.value(key);
		pos_or_prop<typename prp5::vtype,prp5::prop>::value(p5.getVector(),key) = v_exp5.value(key);
		pos_or_prop<typename prp6::vtype,prp6::prop>::value(p6.getVector(),key) = v_exp6.value(key);
		pos_or_prop<typename prp7::vtype,prp7::prop>::value(p7.getVector(),key) = v_exp7.value(key);
		pos_or_prop<typename prp8::vtype,prp8::prop>::value(p8.getVector(),key) = v_exp8.value(key);

		++it;
	}
}

#endif /* OPENFPM_NUMERICS_SRC_OPERATORS_VECTOR_VECTOR_DIST_OPERATOR_ASSIGN_HPP_ */
