//
// Created by Abhinav Singh on 01.12.20.
// in Development

#ifndef OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
#define OPENFPM_NUMERICS_ODEINTEGRATORS_HPP


//namespace std{
//    double abs(pair_ref<double,double> tmp);
//    double abs(const_pair_ref<double,double> tmp);
//}

template<typename T, typename Sfinae = void>
struct has_state_vector: std::false_type {};
template<typename T>
struct has_state_vector<T, typename Void< typename T::is_state_vector>::type> : std::true_type
{};
namespace boost{
    template<class T,class Enabler=typename std::enable_if<has_state_vector<T>::value>::type>
    inline size_t
    size(const T& rng)
    {
        return rng.size();
    }
}

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
struct state_type_1d_ofp{
    state_type_1d_ofp(){
    }
    //Laplacian Lap;
    typedef size_t size_type;
    typedef int is_state_vector;
    aggregate<texp_v<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
    }
};


struct state_type_2d_ofp{
    state_type_2d_ofp(){
    }
    //Laplacian Lap;
    typedef size_t size_type;
    typedef int is_state_vector;
    aggregate<texp_v<double>,texp_v<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
        data.get<1>().resize(n);
    }
};


struct state_type_3d_ofp{
    state_type_3d_ofp(){
    }
    //Laplacian Lap;
    typedef size_t size_type;
    typedef int is_state_vector;
    aggregate<texp_v<double>,texp_v<double>,texp_v<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
        data.get<1>().resize(n);
        data.get<2>().resize(n);
    }
};

namespace boost {
    namespace numeric {
        namespace odeint {
            template<>
            struct is_resizeable<state_type_1d_ofp> {
            typedef boost::true_type type;
            static const bool value = type::value;
            };

            template<>
            struct is_resizeable<state_type_2d_ofp> {
                typedef boost::true_type type;
                static const bool value = type::value;
            };

            template<>
            struct is_resizeable<state_type_3d_ofp> {
                typedef boost::true_type type;
                static const bool value = type::value;
            };



            template<>
            struct vector_space_norm_inf<state_type_1d_ofp>
            {
                typedef double result_type;
            };

            template<>
            struct vector_space_norm_inf<state_type_2d_ofp>
            {
                typedef double result_type;
            };

            template<>
            struct vector_space_norm_inf<state_type_3d_ofp>
            {
                typedef double result_type;
            };

        }
    }
}



#endif //OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
