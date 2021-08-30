/*
 * eq_solve_common.hpp
 *
 *  Created on: Jun 5, 2020
 *      Author: i-bird
 */

#ifndef EQ_SOLVE_COMMON_HPP_
#define EQ_SOLVE_COMMON_HPP_

template<unsigned int prp_id>
struct prop_id {};

class eq_id
{
    int id;

public:

    eq_id()
            :id(0)
    {}

    int getId()
    {
        return id;
    }

    void setId(int id)
    {
        this->id = id;
    }
};


enum options_solver
{
    STANDARD,
    LAGRANGE_MULTIPLIER
};

#endif /* EQ_SOLVE_COMMON_HPP_ */
