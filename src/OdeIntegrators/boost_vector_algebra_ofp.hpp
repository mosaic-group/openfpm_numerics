//
// Created by Abhinav Singh on 18.02.21.
//

#ifndef OPENFPM_PDATA_BOOST_VECTOR_ALGEBRA_OFP_HPP
#define OPENFPM_PDATA_BOOST_VECTOR_ALGEBRA_OFP_HPP

namespace boost {
    namespace numeric {
        namespace odeint {

/*
 * This class template has to be overload in order to call vector_space_algebra::norm_inf
 */
 //           template< class State, class Enabler = void > struct vector_space_norm_inf;

/*
 * Example: instantiation for sole doubles and complex
 */
/*            template<>
            struct vector_space_norm_inf< double >
            {
                typedef double result_type;
                double operator()( double x ) const
                {
                    using std::abs;
                    return abs(x);
                }
            };

            template<>
            struct vector_space_norm_inf< float >
            {
                typedef float result_type;
                result_type operator()( float x ) const
                {
                    using std::abs;
                    return abs(x);
                }
            };

            template< typename T >
            struct vector_space_norm_inf< std::complex<T> >
        {
            typedef T result_type;
            result_type operator()( std::complex<T> x ) const
            {
                using std::abs;
                return abs( x );
            }
        };*/

        /* It copy one element of the chunk for each property
         *
         */
        template<typename vector_type,typename index_type,typename op_type>
        struct for_each_prop1
        {

            vector_type &v;
            index_type &p;
            op_type &op;
            /*! \brief constructor
             *
             *
             * \param src source encapsulated object
             * \param dst destination encapsulated object
             *
             */
            inline for_each_prop1(vector_type &v,index_type &p,op_type &op)
            :v(v),p(p),op(op)
            {};
            //! It call the copy function for each property
            template<typename T>
            inline void operator()(T& t) const
            {

                op(v.data.template get<T::value>().getVector().template get<0>(p));
            }
        };

        struct vector_space_algebra_ofp
        {
            template< class S1 , class Op >
            static void for_each1( S1 &s1 , Op op )
            {

                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop1<S1,size_t,Op> cp(s1,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }


            template<typename S1,typename S2,typename index_type,typename op_type>
            struct for_each_prop2
            {

                S1 &v1;
                S2 &v2;
                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop2(S1 &v1,S2 &v2,index_type &p,op_type &op)
                :v1(v1),v2(v2),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p));
                }
            };
            template< class S1 , class S2 , class Op >
            static void for_each2( S1 &s1 , S2 &s2 , Op op )
            {


                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop2<S1,S2,size_t,Op> cp(s1,s2,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }



            template<typename S1,typename S2,typename S3,typename index_type,typename op_type>
            struct for_each_prop3
            {

                S1 &v1;
                S2 &v2;
                S3 &v3;
                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop3(S1 &v1,S2 &v2,S3 &v3,index_type &p,op_type &op)
                :v1(v1),v2(v2),v3(v3),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p),v3.data.template get<T::value>().getVector().template get<0>(p));
                }
            };
            template< class S1 , class S2 , class S3 , class Op >
            static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
            {


                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop3<S1,S2,S3,size_t,Op> cp(s1,s2,s3,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }


            template<typename S1,typename S2,typename S3,typename S4,typename index_type,typename op_type>
            struct for_each_prop4
            {

                S1 &v1;
                S2 &v2;
                S3 &v3;
                S4 &v4;

                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop4(S1 &v1,S2 &v2,S3 &v3,S4 &v4,index_type &p,op_type &op)
                :v1(v1),v2(v2),v3(v3),v4(v4),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p),v3.data.template get<T::value>().getVector().template get<0>(p),v4.data.template get<T::value>().getVector().template get<0>(p));
                }
            };

            template< class S1 , class S2 , class S3 , class S4 , class Op >
            static void for_each4( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op )
            {

                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop4<S1,S2,S3,S4,size_t,Op> cp(s1,s2,s3,s4,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }


            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename index_type,typename op_type>
            struct for_each_prop5
            {

                S1 &v1;
                S2 &v2;
                S3 &v3;
                S4 &v4;
                S5 &v5;

                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop5(S1 &v1,S2 &v2,S3 &v3,S4 &v4,S5 &v5,index_type &p,op_type &op)
                        :v1(v1),v2(v2),v3(v3),v4(v4),v5(v5),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p),v3.data.template get<T::value>().getVector().template get<0>(p),v4.data.template get<T::value>().getVector().template get<0>(p),v5.data.template get<T::value>().getVector().template get<0>(p));
                }
            };

            template< class S1 , class S2 , class S3 , class S4,class S5 , class Op >
            static void for_each5( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5 , Op op )
            {

                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop5<S1,S2,S3,S4,S5,size_t,Op> cp(s1,s2,s3,s4,s5,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }


            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename index_type,typename op_type>
            struct for_each_prop6
            {

                S1 &v1;
                S2 &v2;
                S3 &v3;
                S4 &v4;
                S5 &v5;
                S6 &v6;


                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop6(S1 &v1,S2 &v2,S3 &v3,S4 &v4,S5 &v5,S6 &v6,index_type &p,op_type &op)
                        :v1(v1),v2(v2),v3(v3),v4(v4),v5(v5),v6(v6),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p),v3.data.template get<T::value>().getVector().template get<0>(p),v4.data.template get<T::value>().getVector().template get<0>(p),v5.data.template get<T::value>().getVector().template get<0>(p),v6.data.template get<T::value>().getVector().template get<0>(p));
                }
            };
            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 , class Op >
            static void for_each6( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6 , Op op )
            {
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop6<S1,S2,S3,S4,S5,S6,size_t,Op> cp(s1,s2,s3,s4,s5,s6,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }



            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename index_type,typename op_type>
            struct for_each_prop7
            {

                S1 &v1;
                S2 &v2;
                S3 &v3;
                S4 &v4;
                S5 &v5;
                S6 &v6;
                S7 &v7;


                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop7(S1 &v1,S2 &v2,S3 &v3,S4 &v4,S5 &v5,S6 &v6,S7 &v7,index_type &p,op_type &op)
                        :v1(v1),v2(v2),v3(v3),v4(v4),v5(v5),v6(v6),v7(v7),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p),v3.data.template get<T::value>().getVector().template get<0>(p),v4.data.template get<T::value>().getVector().template get<0>(p),v5.data.template get<T::value>().getVector().template get<0>(p),v6.data.template get<T::value>().getVector().template get<0>(p),v7.data.template get<T::value>().getVector().template get<0>(p));
                }
            };


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7, class Op >
            static void for_each7( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7 , Op op )
            {
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop7<S1,S2,S3,S4,S5,S6,S7,size_t,Op> cp(s1,s2,s3,s4,s5,s6,s7,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }

            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename index_type,typename op_type>
            struct for_each_prop8
            {

                S1 &v1;
                S2 &v2;
                S3 &v3;
                S4 &v4;
                S5 &v5;
                S6 &v6;
                S7 &v7;
                S8 &v8;


                index_type &p;
                op_type &op;
                /*! \brief constructor
                 *
                 *
                 * \param src source encapsulated object
                 * \param dst destination encapsulated object
                 *
                 */
                inline for_each_prop8(S1 &v1,S2 &v2,S3 &v3,S4 &v4,S5 &v5,S6 &v6,S7 &v7,S8 &v8,index_type &p,op_type &op)
                        :v1(v1),v2(v2),v3(v3),v4(v4),v5(v5),v6(v6),v7(v7),v8(v8),p(p),op(op)
                {};
                //! It call the copy function for each property
                template<typename T>
                inline void operator()(T& t) const
                {
                    op(v1.data.template get<T::value>().getVector().template get<0>(p),v2.data.template get<T::value>().getVector().template get<0>(p),v3.data.template get<T::value>().getVector().template get<0>(p),v4.data.template get<T::value>().getVector().template get<0>(p),v5.data.template get<T::value>().getVector().template get<0>(p),v6.data.template get<T::value>().getVector().template get<0>(p),v7.data.template get<T::value>().getVector().template get<0>(p),v8.data.template get<T::value>().getVector().template get<0>(p));
                }
            };


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class Op >
            static void for_each8( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8 , Op op )
            {
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_prop8<S1,S2,S3,S4,S5,S6,S7,S8,size_t,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,p,op);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);

                    ++it;
                }
            }






           template<typename vector_type,typename index_type,typename norm_result_type>
           struct for_each_norm
           {

               const vector_type &v;
               index_type &p;
               norm_result_type &n;
               /*! \brief constructor
                *
                *
                * \param src source encapsulated object
                * \param dst destination encapsulated object
                *
                */
               inline for_each_norm(const vector_type &v,index_type &p,norm_result_type &n)
               :v(v),p(p),n(n)
               {};
               //! It call the copy function for each property
               template<typename T>
               inline void operator()(T& t) const
               {
                    if(fabs(v.data.template get<T::value>().getVector().template get<0>(p)) > n)
                    {
                        n=fabs(v.data.template get<T::value>().getVector().template get<0>(p));
                    }

               }
           };

            template< class S >
            static typename boost::numeric::odeint::vector_space_norm_inf< S >::result_type norm_inf( const S &s )
            {
                typename boost::numeric::odeint::vector_space_norm_inf< S >::result_type n=0;
                auto it=s.data.template get<0>().getVector().getIterator();
                while(it.isNext()){
                    auto p=it.get();
                    //converting to boost vector ids.
                    for_each_norm<S,size_t,typename boost::numeric::odeint::vector_space_norm_inf< S >::result_type> cp(s,p,n);
                    //creating an iterator on v_ids[0] [1] [2]
                    boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s.data)::max_prop>>(cp);

                    ++it;
                }
                auto &v_cl = create_vcluster();
                v_cl.max(n);
                v_cl.execute();
                //std::max();
                //std::cout<<n<<std::endl;
                return n;
            }
        };



    } // odeint
} // numeric
} // boost








#endif //OPENFPM_PDATA_BOOST_VECTOR_ALGEBRA_OFP_HPP
