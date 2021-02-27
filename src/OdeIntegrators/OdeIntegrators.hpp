//
// Created by Abhinav Singh on 01.12.20.
// in Development

#ifndef OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
#define OPENFPM_NUMERICS_ODEINTEGRATORS_HPP


//namespace std{
//    double abs(pair_ref<double,double> tmp);
//    double abs(const_pair_ref<double,double> tmp);
//}

#include <boost/numeric/odeint.hpp>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "OdeIntegrators/boost_vector_algebra_ofp.hpp"

namespace boost { namespace numeric { namespace odeint {

            template<typename T>
            struct is_resizeable< vector_dist_expression<0,openfpm::vector<aggregate<T>> > >
            {
                typedef boost::true_type type;
                static const bool value = type::value;
            };

        } } }



#endif //OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
