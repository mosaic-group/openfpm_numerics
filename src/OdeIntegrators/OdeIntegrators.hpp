//
// Created by Abhinav Singh on 01.12.20.
// in Development

#ifndef OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
#define OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
#include <boost/numeric/odeint.hpp>
#include "Operators/Vector/vector_dist_operators.hpp"

/*! \brief Resize Temporal Expression for Odeint
  *
  * \param
  *
  * \return
  *
  */
namespace boost { namespace numeric { namespace odeint {

            template<typename T>
            struct is_resizeable< vector_dist_expression<0,openfpm::vector<aggregate<T>> > >
            {
                typedef boost::true_type type;
                static const bool value = type::value;
            };

            template<typename T>
            struct vector_space_norm_inf< vector_dist_expression<0,openfpm::vector<aggregate<T>> > >
            {
                typedef T result_type;
                T operator()( const vector_dist_expression<0,openfpm::vector<aggregate<T>> > &x ) const
                {
                    return norm_inf(x).getReduction();
                }
            };

        } } }



#endif //OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
