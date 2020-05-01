/*
 * FD_op.hpp
 *
 *  Created on: May 1, 2020
 *      Author: Abhinav Singh
 */

#ifndef FD_OP_HPP_
#define FD_OP_HPP_
#include "Derivative.hpp"
#include "eq.hpp"
#include "util/common.hpp"

#include "Operators/Vector/vector_dist_operators.hpp"


template <typename exp1,typename FD_type>
class vector_dist_expression_op<exp1,FD_type,VECT_FD>
{
    //! expression 1
    const exp1 o1;

    FD_type & fd;

public:

    typedef typename exp1::vtype vtype;

    //! Costruct a FD expression out of two expressions
    inline vector_dist_expression_op(const exp1 & o1, FD_type & fd)
            :o1(o1),fd(fd)
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
        return fd.computeDerivative(key,o1);
    }

    /*! \brief Return the vector on which is acting
     *
     * It return the vector used in getVExpr, to get this object
     *
     * \return the vector
     *
     */
    vtype & getVector()
    {
        return o1.getVector();
    }

    /*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    */
    const vtype & getVector() const
    {
        return o1.getVector();
    }
};

class Derivative_x
{

    void * fd;

public:

    template<typename particles_type>
    Derivative_x(particles_type & parts)
    {
        fd = new Fd<particles_type::dims,particles_type>(parts);

    }

    template<typename operand_type>

    vector_dist_expression_op<operand_type,Fd<operand_type::vtype::dims,typename operand_type::vtype>,VECT_FD> operator()(operand_type arg)
    {
        typedef Fd<operand_type::vtype::dims,typename operand_type::vtype> fd_type;

        return vector_dist_expression_op<operand_type,dcpse_type,VECT_FD>(arg,*(fd_type *)fd);
    }
};












////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





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


/*template <typename exp1,typename FD_type>
class vector_dist_expression_op<exp1,DCPSE_type,VECT_FD_V>
{
    //! expression 1
    const exp1 o1;

    FD_type (& fd)[FD_type::vtype::dims];

    static const int dims = FD_type::vtype::dims;
    typedef typename FD_type::vtype::stype stype;

public:

    typedef typename exp1::vtype vtype;

    //! Costruct a subtraction expression out of two expressions
    inline vector_dist_expression_op(const exp1 & o1, FD_type (& fd)[FD_type::vtype::dims])
            :o1(o1),dcp(dcp)
    {}

    *//*! \brief This function must be called before value
     *
     * it initialize the expression if needed
     *
     *//*
    inline void init() const
    {
        o1.init();
    }

    *//*! \brief Evaluate the expression
     *
     * \param key where to evaluate the expression
     *
     * \return the result of the expression
     *
     *//*
    template<typename r_type=VectorS<dims,stype> > inline r_type value(const vect_dist_key_dx & key) const
    {
       VectorS<dims,stype> v_grad;

       for (int i = 0 ; i < dims ; i++)
       {
           v_grad.get(i) = fd[i].computeDifferentialOperator(key,o1);
       }

        return v_grad;
    }
     *//*! \brief Return the vector on which is acting
      *
      * It return the vector used in getVExpr, to get this object
      *
      * \return the vector
      *
      *//*
    vtype & getVector()
    {
        return o1.getVector();
    }

    *//*! \brief Return the vector on which is acting
    *
    * It return the vector used in getVExpr, to get this object
    *
    * \return the vector
    *
    *//*
    const vtype & getVector() const
    {
        return o1.getVector();
    }

};*/


#endif
