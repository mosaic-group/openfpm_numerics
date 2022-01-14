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

/*! \brief A 1d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_1d_ofp{
    state_type_1d_ofp(){
    }
    typedef size_t size_type;
    typedef size_t index_type;
    typedef int is_state_vector;
    aggregate<texp_v<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
    }
};

/*! \brief A 2d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_2d_ofp{
    state_type_2d_ofp(){
    }
    typedef size_t size_type;
    typedef size_t index_type;
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

/*! \brief A 3d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_3d_ofp{
    state_type_3d_ofp(){
    }
    typedef size_t size_type;
    typedef size_t index_type;
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

/*! \brief A 4d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_4d_ofp{
    state_type_4d_ofp(){
    }
    typedef size_t size_type;
    typedef size_t index_type;
    typedef int is_state_vector;
    aggregate<texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
        data.get<1>().resize(n);
        data.get<2>().resize(n);
        data.get<3>().resize(n);
    }
};

/*! \brief A 5d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_5d_ofp{
    state_type_5d_ofp(){
    }
    typedef size_t size_type;
    typedef size_t index_type;
    typedef int is_state_vector;
    aggregate<texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
        data.get<1>().resize(n);
        data.get<2>().resize(n);
        data.get<3>().resize(n);
        data.get<4>().resize(n);

    }
};

template<int counter, typename state_type, typename ... list>
struct state_type_ofpm_add_elements
{
//    typedef aggregate<list ..., texp_v<double>> one_more;
    typedef typename state_type_ofpm_add_elements<counter-1,state_type, state_type,list ...>::type type;
};

template<typename state_type, typename ... list>
struct state_type_ofpm_add_elements<0,state_type,list ...>
{
   typedef aggregate<list ...> type; 
};

template<int n_state, typename state_type>
struct state_type_ofpm_impl
{
    typedef size_t size_type;
    typedef typename state_type::index_type index_type;
    typedef int is_state_vector;

    typedef typename state_type_ofpm_add_elements<n_state-1,state_type, state_type>::type type_data;

    type_data data;
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
            struct is_resizeable<state_type_4d_ofp> {
                typedef boost::true_type type;
                static const bool value = type::value;
            };
            template<>
            struct is_resizeable<state_type_5d_ofp> {
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

            template<>
            struct vector_space_norm_inf<state_type_4d_ofp>
            {
                typedef double result_type;
            };

            template<>
            struct vector_space_norm_inf<state_type_5d_ofp>
            {
                typedef double result_type;
            };

        }
    }
}



#endif //OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
