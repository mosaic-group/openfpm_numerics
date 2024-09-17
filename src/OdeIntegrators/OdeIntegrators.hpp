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
    inline auto
    size(const T& rng) -> decltype(rng.size())
    {
        return rng.size();
    }
}

#include <boost/numeric/odeint.hpp>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "FiniteDifference/FD_expressions.hpp"
#include "OdeIntegrators/vector_algebra_ofp.hpp"

#ifdef __NVCC__
#include "OdeIntegrators/vector_algebra_ofp_gpu.hpp"
/*! \brief A 1d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_1d_ofp_ker{
    state_type_1d_ofp_ker(){
    }
    typedef decltype(std::declval<texp_v_gpu<double>>().getVector().toKernel()) state_kernel;
    typedef size_t size_type;
    typedef int is_state_vector;
    aggregate<state_kernel> data;

    __host__ __device__ size_t size() const
    { return data.get<0>().size(); }

};
/*! \brief A 1d Odeint and Openfpm compatible structure.
 *
 *  Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_1d_ofp_gpu{
    state_type_1d_ofp_gpu(){
    }
    typedef size_t size_type;
    typedef int is_state_vector;
    aggregate<texp_v_gpu<double>> data;

    size_t size() const
    { return data.get<0>().size(); }

    void resize(size_t n)
    {
        data.get<0>().resize(n);
    }
    state_type_1d_ofp_ker toKernel() const
    {
        state_type_1d_ofp_ker s1_ker;
        s1_ker.data.get<0>()=data.get<0>().getVector().toKernel();
        return s1_ker;
    }
};
#endif

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

// On 17.09.24 - foggia
/*! \brief A 6d Odeint and Openfpm compatible structure.
 *
 * Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_6d_ofp{
  state_type_6d_ofp(){
  }
  typedef size_t size_type;
  typedef size_t index_type;
  typedef int is_state_vector;
  aggregate<texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>> data;
  
  size_t size() const
  { return data.get<0>().size(); }
  
  void resize(size_t n)
  {
    data.get<0>().resize(n);
    data.get<1>().resize(n);
    data.get<2>().resize(n);
    data.get<3>().resize(n);
    data.get<4>().resize(n);
    data.get<5>().resize(n);
  }
};

// On 17.09.24 - foggia
/*! \brief A 7d Odeint and Openfpm compatible structure.
 *
 * Use the method this.data.get<d>() to refer to property of all the particles in the dimension d.
 *
 * d starts with 0.
 *
 */
struct state_type_7d_ofp{
  state_type_7d_ofp(){
  }
  typedef size_t size_type;
  typedef size_t index_type;
  typedef int is_state_vector;
  aggregate<texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>,texp_v<double>> data;
  
  size_t size() const
  { return data.get<0>().size(); }
  
  void resize(size_t n)
  {
    data.get<0>().resize(n);
    data.get<1>().resize(n);
    data.get<2>().resize(n);
    data.get<3>().resize(n);
    data.get<4>().resize(n);
    data.get<5>().resize(n);
    data.get<6>().resize(n);
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
    typedef FD::gdb_ext_plus_g_info<state_type::dims> size_type;
    typedef typename state_type::index_type index_type;
    typedef int is_state_vector;

    typedef typename state_type_ofpm_add_elements<n_state-1,state_type, state_type>::type type_data;

    type_data data;

    FD::gdb_ext_plus_g_info<state_type::dims> size() const
    {
        return data.template get<0>().size();
    }


    void resize(const FD::gdb_ext_plus_g_info<state_type::dims> & rsz_obj)
    {
        // to fill
    }
};


namespace boost {
    namespace numeric {
        namespace odeint {

            // FOR particles

            template<>
            struct is_resizeable<state_type_1d_ofp> {
            typedef boost::true_type type;
            static const bool value = type::value;
            };
#ifdef __NVCC__
            template<>
            struct is_resizeable<state_type_1d_ofp_gpu> {
                typedef boost::true_type type;
                static const bool value = type::value;
            };
#endif
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

	  // On 17.09.24 - foggia
	  template<>
	  struct is_resizeable<state_type_6d_ofp> {
	    typedef boost::true_type type;
	    static const bool value = type::value;
	  };

	  // On 17.09.24 - foggia
	  template<>
	  struct is_resizeable<state_type_7d_ofp> {
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

	  // On 17.09.24 - foggia
	  template<>
	  struct vector_space_norm_inf<state_type_6d_ofp>
	  {
	    typedef double result_type;
	  };

	  // On 17.09.24 - foggia
	  template<>
	  struct vector_space_norm_inf<state_type_7d_ofp>
	  {
	    typedef double result_type;
	  };

	    // For GRIDs

            template<typename state_type>
            struct is_resizeable<state_type_ofpm_impl<1,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
            };

            template<typename state_type>
            struct is_resizeable<state_type_ofpm_impl<2,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
            };

            template<typename state_type>
            struct is_resizeable<state_type_ofpm_impl<3,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
            };

            template<typename state_type>
            struct is_resizeable<state_type_ofpm_impl<4,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
            };

            template<typename state_type>
            struct is_resizeable<state_type_ofpm_impl<5,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
            };

	  // On 17.09.24 - foggia
	  template<typename state_type>
	  struct is_resizeable<state_type_ofpm_impl<6,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
	  };

	  // On 17.09.24 - foggia
	  template<typename state_type>
	  struct is_resizeable<state_type_ofpm_impl<7,state_type> > {
            typedef boost::true_type type;
            static const bool value = type::value;
	  };


/*            template<>
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
            };*/



/*      //      template<unsigned int nprp, typename state_type>
            struct vector_space_norm_inf<state_type_ofpm_impl<nprp,state_type>>
            {
                typedef double result_type;
            };*/

            template<typename state_type>
            struct vector_space_norm_inf<state_type_ofpm_impl<1,state_type>>
            {
                typedef double result_type;
            };

            template<typename state_type>
            struct vector_space_norm_inf<state_type_ofpm_impl<2,state_type>>
            {
                typedef double result_type;
            };

            template<typename state_type>
            struct vector_space_norm_inf<state_type_ofpm_impl<3,state_type>>
            {
                typedef double result_type;
            };

            template<typename state_type>
            struct vector_space_norm_inf<state_type_ofpm_impl<4,state_type>>
            {
                typedef double result_type;
            };

            template<typename state_type>
            struct vector_space_norm_inf<state_type_ofpm_impl<5,state_type>>
            {
                typedef double result_type;
            };

	  // On 17.09.24 - foggia
	  template<typename state_type>
	  struct vector_space_norm_inf<state_type_ofpm_impl<6,state_type>>
	  {
	    typedef double result_type;
	  };

	  // On 17.09.24 - foggia
	  template<typename state_type>
	  struct vector_space_norm_inf<state_type_ofpm_impl<7,state_type>>
	  {
	    typedef double result_type;
	  };


        }
    }
}



#endif //OPENFPM_NUMERICS_ODEINTEGRATORS_HPP
