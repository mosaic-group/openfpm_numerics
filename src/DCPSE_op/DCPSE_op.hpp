/*
 * DCPSE_op.hpp
 *
 *  Created on: Jan 7, 2020
 *      Author: i-bird
 */

#ifndef DCPSE_OP_HPP_
#define DCPSE_OP_HPP_

template<typename operand_type>
class Lap_node
{
	operand_type arg;

public:

	typedef int it_is_a_node;

	Lap_node(operand_type & arg)
	:arg(arg)
	{}
};

class Lap
{
public:

	template<typename operand_type>
	Lap_node<operand_type> operator()(operand_type arg)
	{
		return Lap_node<operand_type>(arg);
	}
};



template<typename operand_type1, typename operand_type2>
class plus
{
	operand_type1 op1;
	operand_type2 op2;

public:

	typedef int it_is_a_node;

	plus(const operand_type1 & op1, const operand_type2 & op2)
	:op1(op1),op2(op2)
	{}

	void value()
	{
		//op1.value() + op.value;
	}

};

class Field
{
	typedef int it_is_a_node;

	void value()
	{
		// add non zero
	}
};

template<typename operand_type1, typename operand_type2/*, typename sfinae=typename std::enable_if<
																						std::is_same<typename operand_type1::it_is_a_node,int>::value
																						>::type*/ >
plus<operand_type1,operand_type2> operator+(const operand_type1 & op1, const operand_type2 & op2)
{
	return plus<operand_type1,operand_type2>(op1,op2);
}



#endif /* DCPSE_OP_HPP_ */
