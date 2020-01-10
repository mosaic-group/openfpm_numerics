/*
 * DCPSE_op.hpp
 *
 *  Created on: Jan 7, 2020
 *      Author: Abhinav Singh, Pietro Incardona
 */

#ifndef DCPSE_OP_HPP_
#define DCPSE_OP_HPP_
#include "DCPSE_op/DCPSE_copy/Dcpse.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

/*! \brief Subtraction operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename exp1,typename DCPSE_type>
class vector_dist_expression_op<exp1,DCPSE_type,VECT_DCPSE>
{
    //! expression 1
    const exp1 o1;

    DCPSE_type & dcp;

public:

    //! Costruct a subtraction expression out of two expressions
    inline vector_dist_expression_op(const exp1 & o1, DCPSE_type & dcp)
            :o1(o1),dcp(dcp)
    {}

    /*! \brief This function must be called before value
     *
     * it initialize the expression if needed
     *
     */
    inline void init() const
    {
        o1.init();
    }

    /*! \brief Evaluate the expression
     *
     * \param key where to evaluate the expression
     *
     * \return the result of the expression
     *
     */
    template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()))>::type > inline r_type value(const vect_dist_key_dx & key) const
    {
        return dcp.computeDifferentialOperator(key,o1);
    }
};
/*
template<typename operand_type>
class Derivative_x_node
{
	operand_type arg;

public:
	typedef int it_is_a_node;
	Derivative_x_node(operand_type &arg)
	:arg(arg)
	{}
};
*/

class Derivative_x
{

    void * dcpse;

public:

    template<typename particles_type>
    Derivative_x(particles_type & parts, unsigned int ord ,typename particles_type::stype rCut)
    {
        Point<particles_type::dims,unsigned int> p;
        p.zero();
        p.get(0) = 1;

        dcpse = new Dcpse<particles_type::dims,particles_type>(parts,p, ord, rCut,2.5);
    }

    template<typename operand_type>

    vector_dist_expression_op<operand_type,Dcpse<operand_type::vtype::dims,typename operand_type::vtype>,VECT_DCPSE> operator()(operand_type arg)
    {
        typedef Dcpse<operand_type::vtype::dims,typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type,dcpse_type,VECT_DCPSE>(arg,*(dcpse_type *)dcpse);
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

//template<typename operand_type1, typename operand_type2/*, typename sfinae=typename std::enable_if<
//																						std::is_same<typename operand_type1::it_is_a_node,int>::value
//																						>::type*/ >
//plus<operand_type1,operand_type2> operator+(const operand_type1 & op1, const operand_type2 & op2)
//{
//	return plus<operand_type1,operand_type2>(op1,op2);
//}
#endif /* DCPSE_OP_HPP_ */
