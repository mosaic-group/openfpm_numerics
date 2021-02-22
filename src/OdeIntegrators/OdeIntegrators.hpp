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


template<typename tvar1, typename tvar2>
struct pair_ref;

template<typename tvar1, typename tvar2>
struct pair_val
{
    tvar1 var1;
    tvar2 var2;
    pair_val()
    {}
    pair_val(const tvar1 & k, const tvar2 & v)
            :var1(k),var2(v)
    {}
    pair_val(const pair_ref<tvar1,tvar2> & tmp)
    {
        this->operator=(tmp);
    }
    bool operator<(const pair_val & tmp) const
    {
        return var1 < tmp.var1;
    }
    bool operator>(const pair_val & tmp) const
    {
        return var1 > tmp.var1;
    }
    pair_val & operator=(const pair_ref<tvar1,tvar2> & tmp)
    {
        var1 = tmp.var1;
        var2 = tmp.var2;
        return *this;
    }
    pair_val operator+(const pair_val<tvar1,tvar2> & tmp)
    {   pair_val result;
        result.var1 = var1+tmp.var1;
        result.var2 = var2+tmp.var2;
        return result;
    }
};

template<typename tvar1, typename tvar2>
struct pair_ref
{
    tvar1 & var1;
    tvar2 & var2;
    pair_ref(tvar1 & k, tvar2 & v)
            :var1(k),var2(v)
    {}
    pair_ref(pair_ref<tvar1,tvar2> && tmp)
            :var1(tmp.var1),var2(tmp.var2)
    {}

    pair_ref(const pair_ref<tvar1,tvar2> & tmp)
            :var1(tmp.var1),var2(tmp.var2)
    {}

    pair_ref & operator=(const pair_val<tvar1,tvar2> & tmp)
    {
        var1 = tmp.var1;
        var2 = tmp.var2;
        return *this;
    }
    pair_ref & operator=(const pair_ref<tvar1,tvar2> & tmp)
    {
        var1 = tmp.var1;
        var2 = tmp.var2;
        return *this;
    }
};

template<typename tvar1, typename tvar2>
struct const_pair_ref
{
    const tvar1 & var1;
    const tvar2 & var2;
    const_pair_ref(const tvar1 & k,const tvar2 & v)
            :var1(k),var2(v)
    {}
    const_pair_ref(const_pair_ref<tvar1,tvar2> && tmp)
            :var1(tmp.var1),var2(tmp.var2)
    {}
};

template<typename tvar1, typename tvar2>
struct pair_it
{

    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t difference_type;
    typedef pair_val<tvar1,tvar2> value_type;
    typedef pair_it<tvar1,tvar2> pointer;
    typedef pair_ref<tvar1,tvar2> reference;

    tvar1 * var1;
    tvar2 * var2;

    tvar1 v1;
    tvar2 v2;
    pair_ref<tvar1,tvar2> tmp_ref;

    bool operator==(const pair_it & tmp)
    {
        return (var1 == tmp.var1 && var2 == tmp.var2);
    }
    pair_ref<tvar1,tvar2>& operator*()
    {
        tmp_ref=pair_ref<tvar1,tvar2>(*var1,*var2);
        return tmp_ref;
    }

    pair_ref<tvar1,tvar2> operator[](int i)
    {
        return pair_ref<tvar1,tvar2>(*var1,*var2);
    }
    pair_it operator+(size_t count) const
    {
        pair_it tmp(var1+count,var2+count);
        return tmp;
    }
    size_t operator-(pair_it & tmp) const
    {
        return var1 - tmp.var1;
    }
    pair_it operator-(size_t count) const
    {
        pair_it tmp(var1-count,var2-count);
        return tmp;
    }
    pair_it & operator++()
    {
        ++var1;
        ++var2;
        return *this;
    }
    pair_it operator++(int)
    {
        pair_it tmp(*this);
        ++var1;
        ++var2;
        return tmp;
    }
    pair_it & operator--()
    {
        --var1;
        --var2;
        return *this;
    }
    pair_it operator--(int)
    {
        pair_it tmp(*this);
        --var1;
        --var2;
        return tmp;
    }
    bool operator!=(const pair_it & tmp) const
    {
        return var1 != tmp.var1 && var2 != tmp.var2;
    }
    bool operator<(const pair_it & tmp) const
    {
        return var1 < tmp.var1;
    }
    pair_it<tvar1,tvar2> & operator=(pair_it<tvar1,tvar2> & tmp)
    {
        var1 = tmp.var1;
        var2 = tmp.var2;
        return *this;
    }
    pair_it()
            :tmp_ref(v1,v2)
    {}
    pair_it(const pair_it<tvar1,tvar2> & tmp)
            :var1(tmp.var1),var2(tmp.var2),tmp_ref(*tmp.var1,*tmp.var2)
    {}
    pair_it(tvar1 * var1, tvar2 * var2)
            :var1(var1),var2(var2),tmp_ref(*var1,*var2)
    {}
};

template<typename tvar1, typename tvar2>
struct const_pair_it
{
    typedef std::random_access_iterator_tag iterator_category;
    typedef size_t difference_type;
    typedef const pair_val<tvar1,tvar2> value_type;
    typedef const_pair_it<tvar1,tvar2> pointer;
    typedef const_pair_ref<tvar1,tvar2> reference;

    const tvar1 * var1;
    const tvar2 * var2;
    const tvar1 v1=0;
    const tvar2 v2=0;

    const_pair_ref<tvar1,tvar2> tmp_ref;

    bool operator==(const const_pair_it & tmp)
    {
        return (var1 == tmp.var1 && var2 == tmp.var2);
    }
    const_pair_ref<tvar1,tvar2> operator*()
    {
        //tmp_ref=const_pair_ref<tvar1,tvar2>(*var1,*var2);
        return const_pair_ref<tvar1,tvar2>(*var1,*var2);
    }

    const_pair_ref<tvar1,tvar2> operator[](int i)
    {
        return pair_ref<tvar1,tvar2>(*var1,*var2);
    }
    const_pair_it operator+(size_t count) const
    {
        const_pair_it tmp(var1+count,var2+count);
        return tmp;
    }
    size_t operator-(const_pair_it & tmp) const
    {
        return var1 - tmp.var1;
    }
    const_pair_it operator-(size_t count) const
    {
        const_pair_it tmp(var1-count,var2-count);
        return tmp;
    }
    const_pair_it & operator++()
    {
        ++var1;
        ++var2;
        return *this;
    }
    const_pair_it operator++(int)
    {
        const_pair_it tmp(*this);
        ++var1;
        ++var2;
        return tmp;
    }
    const_pair_it & operator--()
    {
        --var1;
        --var2;
        return *this;
    }
    const_pair_it operator--(int)
    {
        const_pair_it tmp(*this);
        --var1;
        --var2;
        return tmp;
    }
    bool operator!=(const const_pair_it & tmp) const
    {
        return var1 != tmp.var1 && var2 != tmp.var2;
    }
    bool operator<(const const_pair_it & tmp) const
    {
        return var1 < tmp.var1;
    }
    const_pair_it<tvar1,tvar2> & operator=(const_pair_it<tvar1,tvar2> & tmp)
    {
        var1 = tmp.var1;
        var2 = tmp.var2;
        return *this;
    }
    const_pair_it():tmp_ref(v1,v2)
    {}
    const_pair_it(const const_pair_it<tvar1,tvar2> & tmp)
            :var1(tmp.var1),var2(tmp.var2),tmp_ref(*tmp.var1,*tmp.var2)
    {}
    const_pair_it(const tvar1 * var1,const tvar2 * var2)
            :var1(var1),var2(var2),tmp_ref(*var1,*var2)
    {}
};

template<typename tvar1,typename tvar2>
pair_val<tvar1,tvar2> operator*(double factor,const pair_ref<tvar1,tvar2>& pair)
{
    return pair_val<tvar1,tvar2>(factor*pair.var1,factor*pair.var2);
}

template<typename tvar1,typename tvar2>
pair_val<tvar1,tvar2> operator*(double factor,const const_pair_ref<tvar1,tvar2>& pair)
{
    return pair_val<tvar1,tvar2>(factor*pair.var1,factor*pair.var2);
}


namespace boost { namespace numeric { namespace odeint {

            template<typename T>
            struct is_resizeable< vector_dist_expression<0,openfpm::vector<aggregate<T>> > >
            {
                typedef boost::true_type type;
                static const bool value = type::value;
            };

            /*template<typename T>
            struct vector_space_norm_inf< vector_dist_expression<0,openfpm::vector<aggregate<T>> > >
            {
                typedef T result_type;
                T operator()( const vector_dist_expression<0,openfpm::vector<aggregate<T>> > &x ) const
                {
                    return norm_inf(x).getReduction();
                }
            };*/

        } } }



#endif //OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
