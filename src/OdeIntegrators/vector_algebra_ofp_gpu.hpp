//
// Created by Abhinav Singh on 1.03.23.
//

#ifndef OPENFPM_PDATA_VECTOR_ALGEBRA_OFP_GPU_HPP
#define OPENFPM_PDATA_VECTOR_ALGEBRA_OFP_GPU_HPP

namespace boost {
    namespace numeric {
        namespace odeint {

            template<typename S1, typename Op>
            __global__ void for_each1_ker(S1 s1, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.size())	{return;}

                for_each_prop1<S1,size_t,Op> cp(s1,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }

            template<typename S1,typename S2, typename Op>
            __global__ void for_each2_ker(S1 s1,S2 s2, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //printf("%f \n",s2.data.template get<0>().getVector().template get<0>(p));
                //converting to boost vector ids.
                for_each_prop2<S1,S2,unsigned int,Op> cp(s1,s2,p,op);
                //s1.data.template get<0>().getVector().template get<0>(p)=1.0*s2.data.template get<0>().getVector().template get<0>(p)+0.05*s3.data.template get<0>().getVector().template get<0>(p);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }

            template<typename S1,typename S2,typename S3, typename Op>
            __global__ void for_each3_ker(S1 s1,S2 s2,S3 s3, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //printf("%f \n",s2.data.template get<0>().getVector().template get<0>(p));
                //converting to boost vector ids.
                for_each_prop3<S1,S2,S3,unsigned int,Op> cp(s1,s2,s3,p,op);
                //s1.data.template get<0>().getVector().template get<0>(p)=1.0*s2.data.template get<0>().getVector().template get<0>(p)+0.05*s3.data.template get<0>().getVector().template get<0>(p);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }

            template<typename S1,typename S2,typename S3,typename S4, typename Op>
            __global__ void for_each4_ker(S1 s1,S2 s2,S3 s3,S4 s4, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //printf("%f \n",s2.data.template get<0>().getVector().template get<0>(p));
                //converting to boost vector ids.
                for_each_prop4<S1,S2,S3,S4,unsigned int,Op> cp(s1,s2,s3,s4,p,op);
                //s1.data.template get<0>().getVector().template get<0>(p)=1.0*s2.data.template get<0>().getVector().template get<0>(p)+0.05*s3.data.template get<0>().getVector().template get<0>(p);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }

            template<typename S1,typename S2,typename S3,typename S4,typename S5, typename Op>
            __global__ void for_each5_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //printf("%f \n",s2.data.template get<0>().getVector().template get<0>(p));
                //converting to boost vector ids.
                for_each_prop5<S1,S2,S3,S4,S5,unsigned int,Op> cp(s1,s2,s3,s4,s5,p,op);
                //s1.data.template get<0>().getVector().template get<0>(p)=1.0*s2.data.template get<0>().getVector().template get<0>(p)+0.05*s3.data.template get<0>().getVector().template get<0>(p);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }

            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6, typename Op>
            __global__ void for_each6_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop6<S1,S2,S3,S4,S5,S6,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7, typename Op>
            __global__ void for_each7_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop7<S1,S2,S3,S4,S5,S6,S7,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8, typename Op>
            __global__ void for_each8_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop8<S1,S2,S3,S4,S5,S6,S7,S8,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9, typename Op>
            __global__ void for_each9_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop9<S1,S2,S3,S4,S5,S6,S7,S8,S9,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9,typename S10, typename Op>
            __global__ void for_each10_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9,S10 s10, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop10<S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9,typename S10,typename S11, typename Op>
            __global__ void for_each11_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9,S10 s10,S11 s11, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop11<S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9,typename S10,typename S11,typename S12, typename Op>
            __global__ void for_each12_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9,S10 s10,S11 s11,S12 s12, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop12<S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9,typename S10,typename S11,typename S12,typename S13, typename Op>
            __global__ void for_each13_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9,S10 s10,S11 s11,S12 s12,S13 s13, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop13<S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }
            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9,typename S10,typename S11,typename S12,typename S13,typename S14, typename Op>
            __global__ void for_each14_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9,S10 s10,S11 s11,S12 s12,S13 s13,S14 s14, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop14<S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }


            template<typename S1,typename S2,typename S3,typename S4,typename S5,typename S6,typename S7,typename S8,typename S9,typename S10,typename S11,typename S12,typename S13,typename S14,typename S15, typename Op>
            __global__ void for_each15_ker(S1 s1,S2 s2,S3 s3,S4 s4,S5 s5,S6 s6,S7 s7,S8 s8,S9 s9,S10 s10,S11 s11,S12 s12,S13 s13,S14 s14,S15 s15, Op op)
            {
                unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

                if (p >= s1.data.template get<0>().size())	{return;}
                //converting to boost vector ids.
                for_each_prop15<S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,unsigned int,Op> cp(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,p,op);
                //creating an iterator on v_ids[0] [1] [2]
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(cp);
            }

        struct vector_space_algebra_ofp_gpu
        {

            template< class S1 , class Op >
            static void for_each1( S1 &s1 , Op op )
            {
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getDomainIteratorGPU();
                CUDA_LAUNCH((for_each1_ker),it,s1,op);
            }
            template< class S1 , class S2 , class Op >
            static void for_each2( S1 &s1 , S2 &s2 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                //s1.data.template get<0>().getVector().resize(s2.data.template get<0>().getVector().size());
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each2_ker),it,s1.toKernel(),s2.toKernel(),op);
            }

            template< class S1 , class S2 , class S3 , class Op >
            static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype( s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each3_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),op);

            }

            template< class S1 , class S2 , class S3 , class S4 , class Op >
            static void for_each4( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each4_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),op);
            }


            template< class S1 , class S2 , class S3 , class S4,class S5 , class Op >
            static void for_each5( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each5_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),op);
            }

            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 , class Op >
            static void for_each6( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each6_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),op);
            }


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7, class Op >
            static void for_each7( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each7_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),op);
            }

            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class Op >
            static void for_each8( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each8_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),op);
            }

            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class Op >
            static void for_each9( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9 , Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each9_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),s9.toKernel(),op);
            }


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class S10, class Op >
            static void for_each10( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9 , S10 &s10, Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each10_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),s9.toKernel(),s10.toKernel(),op);

            }


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class S10, class S11, class Op >
            static void for_each11( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9 , S10 &s10,S11 &s11, Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each11_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),op);
            }


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class S10, class S11, class S12, class Op >
            static void for_each12( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9 , S10 &s10,S11 &s11,S12 &s12, Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each12_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),s9.toKernel(),s10.toKernel(),s11.toKernel(),s12.toKernel(),op);
            }


            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class S10, class S11, class S12, class S13, class Op >
            static void for_each13( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9 , S10 &s10,S11 &s11,S12 &s12,S13 &s13, Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each13_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),s9.toKernel(),s10.toKernel(),s11.toKernel(),s12.toKernel(),s13.toKernel(),op);
            }

            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class S10, class S11, class S12, class S13, class S14, class Op >
            static void for_each14( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9, S10 &s10,S11 &s11,S12 &s12,S13 &s13,S14 &s14, Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                // ToDo : build checks, that the +-*/ operators are well defined
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each14_ker),it,it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),s9.toKernel(),s10.toKernel(),s11.toKernel(),s12.toKernel(),s13.toKernel(),s14.toKernel(),op);
            }

            template< class S1 , class S2 , class S3 , class S4,class S5,class S6 ,class S7,class S8, class S9, class S10, class S11, class S12, class S13, class S14, class S15, class Op >
            static void for_each15( S1 &s1 , S2 &s2 , S3 &s3 , S4 &s4,S5 &s5,S6 &s6,S7 &s7,S8 &s8, S9 &s9, S10 &s10,S11 &s11,S12 &s12,S13 &s13,S14 &s14,S15 &s15, Op op )
            {
                for_each_prop_resize<S1,S2> the_resize(s1,s2);
                boost::mpl::for_each_ref<boost::mpl::range_c<int,0,decltype(s1.data)::max_prop>>(the_resize);
                auto it=s1.data.template get<0>().getVector().getGPUIterator();

                CUDA_LAUNCH((for_each15_ker),it,s1.toKernel(),s2.toKernel(),s3.toKernel(),s4.toKernel(),s5.toKernel(),s6.toKernel(),s7.toKernel(),s8.toKernel(),s9.toKernel(),s10.toKernel(),s11.toKernel(),s12.toKernel(),s13.toKernel(),s14.toKernel(),s15.toKernel(),op);
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


#include <algorithm>

#include <boost/config.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint/util/unit_helper.hpp>





/*
 * Notes:
 *
 * * the results structs are needed in order to work with fusion_algebra
 */
struct ofp_operations
{

    template< class Fac1 = double >
    struct scale
    {
        const Fac1 m_alpha1;

        scale( Fac1 alpha1 ) : m_alpha1( alpha1 ) { }

        template< class T1 >
 __device__ __host__       void operator()( T1 &t1 ) const
        {
            t1 *= m_alpha1;
        }

        typedef void result_type;
    };

    template< class Fac1 = double >
    struct scale_sum1
    {
        const Fac1 m_alpha1;

        scale_sum1( Fac1 alpha1 ) : m_alpha1( alpha1 ) { }

        template< class T1 , class T2 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 ) const
        {
            t1 = m_alpha1 * t2;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 >
    struct scale_sum2
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum2( Fac1 alpha1 , Fac2 alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

        template< class T1 , class T2 , class T3 >
        BOOST_FUSION_GPU_ENABLED
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 >
    struct scale_sum3
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;

        scale_sum3( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) { }

        template< class T1 , class T2 , class T3 , class T4 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 >
    struct scale_sum4
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;

        scale_sum4( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 >
    struct scale_sum5
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;

        scale_sum5( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 , Fac5 alpha5 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 >
    struct scale_sum6
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;

        scale_sum6( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 , Fac5 alpha5 , Fac6 alpha6 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ){ }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 ,const T7 &t7) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 >
    struct scale_sum7
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;

        scale_sum7( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 >
    struct scale_sum8
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;

        scale_sum8( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9;
        }

        typedef void result_type;
    };

    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 >
    struct scale_sum9
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;

        scale_sum9( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 , class T10 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 , const T10 &t10 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9 + m_alpha9 * t10;
        }

        typedef void result_type;
    };

    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 >
    struct scale_sum10
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;

        scale_sum10( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 , Fac10 alpha10 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 , class T10 , class T11 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 , const T10 &t10 , const T11 &t11 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9 + m_alpha9 * t10 + m_alpha10 * t11;
        }

        typedef void result_type;
    };


    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 >
    struct scale_sum11
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;

        scale_sum11( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 , class T10 , class T11 , class T12 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 , const T10 &t10 , const T11 &t11 , const T12 &t12 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9 + m_alpha9 * t10 + m_alpha10 * t11 + m_alpha11 * t12;
        }

        typedef void result_type;
    };

    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 , class Fac12 = Fac11 >
    struct scale_sum12
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;
        const Fac12 m_alpha12;

        scale_sum12( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 , Fac12 alpha12 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) , m_alpha12( alpha12 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 , class T10 , class T11 , class T12 , class T13 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 , const T10 &t10 , const T11 &t11 , const T12 &t12 , const T13 &t13 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9 + m_alpha9 * t10 + m_alpha10 * t11 + m_alpha11 * t12 + m_alpha12 * t13;
        }

        typedef void result_type;
    };

    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 , class Fac12 = Fac11 , class Fac13 = Fac12 >
    struct scale_sum13
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;
        const Fac12 m_alpha12;
        const Fac13 m_alpha13;

        scale_sum13( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 , Fac12 alpha12 , Fac13 alpha13 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) , m_alpha12( alpha12 ) , m_alpha13( alpha13 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 , class T10 , class T11 , class T12 , class T13 , class T14 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 , const T10 &t10 , const T11 &t11 , const T12 &t12 , const T13 &t13 , const T14 &t14 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9 + m_alpha9 * t10 + m_alpha10 * t11 + m_alpha11 * t12 + m_alpha12 * t13 + m_alpha13 * t14;
        }

        typedef void result_type;
    };

    template< class Fac1 = double , class Fac2 = Fac1 , class Fac3 = Fac2 , class Fac4 = Fac3 , class Fac5 = Fac4 , class Fac6 = Fac5 , class Fac7 = Fac6 , class Fac8 = Fac7 , class Fac9 = Fac8 , class Fac10 = Fac9 , class Fac11 = Fac10 , class Fac12 = Fac11 , class Fac13 = Fac12 , class Fac14 = Fac13 >
    struct scale_sum14
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;
        const Fac3 m_alpha3;
        const Fac4 m_alpha4;
        const Fac5 m_alpha5;
        const Fac6 m_alpha6;
        const Fac7 m_alpha7;
        const Fac8 m_alpha8;
        const Fac9 m_alpha9;
        const Fac10 m_alpha10;
        const Fac11 m_alpha11;
        const Fac12 m_alpha12;
        const Fac13 m_alpha13;
        const Fac14 m_alpha14;

        scale_sum14( Fac1 alpha1 , Fac2 alpha2 , Fac3 alpha3 , Fac4 alpha4 ,
                Fac5 alpha5 , Fac6 alpha6 , Fac7 alpha7 , Fac8 alpha8 , Fac9 alpha9 ,
                Fac10 alpha10 , Fac11 alpha11 , Fac12 alpha12 , Fac13 alpha13 , Fac14 alpha14 )
        : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) , m_alpha3( alpha3 ) , m_alpha4( alpha4 ) , m_alpha5( alpha5 ) , m_alpha6( alpha6 ) , m_alpha7( alpha7 ) , m_alpha8( alpha8 ) , m_alpha9( alpha9 ) , m_alpha10( alpha10 ) , m_alpha11( alpha11 ) , m_alpha12( alpha12 ) , m_alpha13( alpha13 ) , m_alpha14( alpha14 ) { }

        template< class T1 , class T2 , class T3 , class T4 , class T5 , class T6 , class T7 , class T8 , class T9 , class T10 , class T11 , class T12 , class T13 , class T14 , class T15 >
 __device__ __host__       void operator()( T1 &t1 , const T2 &t2 , const T3 &t3 , const T4 &t4 , const T5 &t5 , const T6 &t6 , const T7 &t7 , const T8 &t8 , const T9 &t9 , const T10 &t10 , const T11 &t11 , const T12 &t12 , const T13 &t13 , const T14 &t14 , const T15 &t15 ) const
        {
            t1 = m_alpha1 * t2 + m_alpha2 * t3 + m_alpha3 * t4 + m_alpha4 * t5 + m_alpha5 * t6 + m_alpha6 * t7 + m_alpha7 * t8 + m_alpha8 * t9 + m_alpha9 * t10 + m_alpha10 * t11 + m_alpha11 * t12 + m_alpha12 * t13 + m_alpha13 * t14 + m_alpha14 * t15;
        }

        typedef void result_type;
    };

    template< class Fac1 = double , class Fac2 = Fac1 >
    struct scale_sum_swap2
    {
        const Fac1 m_alpha1;
        const Fac2 m_alpha2;

        scale_sum_swap2( Fac1 alpha1 , Fac2 alpha2 ) : m_alpha1( alpha1 ) , m_alpha2( alpha2 ) { }

        template< class T1 , class T2 , class T3 >
 __device__ __host__       void operator()( T1 &t1 , T2 &t2 , const T3 &t3) const
        {
            const T1 tmp( t1 );
            t1 = m_alpha1 * t2 + m_alpha2 * t3;
            t2 = tmp;
        }

        typedef void result_type;
    };

    /*
     * for usage in for_each2
     *
     * Works with boost::units by eliminating the unit
     */
    template< class Fac1 = double >
    struct rel_error
    {
        const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

        rel_error( Fac1 eps_abs , Fac1 eps_rel , Fac1 a_x , Fac1 a_dxdt )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt ) { }


        template< class T1 , class T2 , class T3 >
 __device__ __host__       void operator()( T3 &t3 , const T1 &t1 , const T2 &t2 ) const
        {
            using std::abs;
            set_unit_value( t3 , abs( get_unit_value( t3 ) ) / ( m_eps_abs + m_eps_rel * ( m_a_x * abs( get_unit_value( t1 ) ) + m_a_dxdt * abs( get_unit_value( t2 ) ) ) ) );
        }

        typedef void result_type;
    };


    /*
     * for usage in for_each3
     *
     * used in the controller for the rosenbrock4 method
     *
     * Works with boost::units by eliminating the unit
     */
    template< class Fac1 = double >
    struct default_rel_error
    {
        const Fac1 m_eps_abs , m_eps_rel ;

        default_rel_error( Fac1 eps_abs , Fac1 eps_rel )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) { }


        /*
         * xerr = xerr / ( eps_abs + eps_rel * max( x , x_old ) )
         */
        template< class T1 , class T2 , class T3 >
 __device__ __host__       void operator()( T3 &t3 , const T1 &t1 , const T2 &t2 ) const
        {
            BOOST_USING_STD_MAX();
            using std::abs;
            Fac1 x1 = abs( get_unit_value( t1 ) ) , x2 = abs( get_unit_value( t2 ) );
            set_unit_value( t3 , abs( get_unit_value( t3 ) ) / ( m_eps_abs + m_eps_rel * max BOOST_PREVENT_MACRO_SUBSTITUTION ( x1 , x2 ) ) );
        }

        typedef void result_type;
    };



    /*
     * for usage in reduce
     */

    template< class Value >
    struct maximum
    {
        template< class Fac1 , class Fac2 >
  __device__ __host__      Value operator()( Fac1 t1 , const Fac2 t2 ) const
        {
            using std::abs;
            Value a1 = abs( get_unit_value( t1 ) ) , a2 = abs( get_unit_value( t2 ) );
            return ( a1 < a2 ) ? a2 : a1 ;
        }

        typedef Value result_type;
    };


    template< class Fac1 = double >
    struct rel_error_max
    {
        const Fac1 m_eps_abs , m_eps_rel;

        rel_error_max( Fac1 eps_abs , Fac1 eps_rel )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel )
        { }

        template< class Res , class T1 , class T2 , class T3 >
__device__ __host__        Res operator()( Res r , const T1 &x_old , const T2 &x , const T3 &x_err )
        {
            BOOST_USING_STD_MAX();
            using std::abs;
            Res tmp = abs( get_unit_value( x_err ) ) / ( m_eps_abs + m_eps_rel * max BOOST_PREVENT_MACRO_SUBSTITUTION ( abs( x_old ) , abs( x ) ) );
            return max BOOST_PREVENT_MACRO_SUBSTITUTION ( r , tmp );
        }
    };


    template< class Fac1 = double >
    struct rel_error_max2
    {
        const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

        rel_error_max2( Fac1 eps_abs , Fac1 eps_rel , Fac1 a_x , Fac1 a_dxdt )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
        { }

        template< class Res , class T1 , class T2 , class T3 , class T4 >
__device__ __host__        Res operator()( Res r , const T1 &x_old , const T2 &/*x*/ , const T3 &dxdt_old , const T4 &x_err )
        {
            BOOST_USING_STD_MAX();
            using std::abs;
            Res tmp = abs( get_unit_value( x_err ) ) /
                    ( m_eps_abs + m_eps_rel * ( m_a_x * abs( get_unit_value( x_old ) ) + m_a_dxdt * abs( get_unit_value( dxdt_old ) ) ) );
            return max BOOST_PREVENT_MACRO_SUBSTITUTION ( r , tmp );
        }
    };




    template< class Fac1 = double >
    struct rel_error_l2
    {
        const Fac1 m_eps_abs , m_eps_rel;

        rel_error_l2( Fac1 eps_abs , Fac1 eps_rel )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel )
        { }

        template< class Res , class T1 , class T2 , class T3 >
__device__ __host__        Res operator()( Res r , const T1 &x_old , const T2 &x , const T3 &x_err )
        {
            BOOST_USING_STD_MAX();
            using std::abs;
            Res tmp = abs( get_unit_value( x_err ) ) / ( m_eps_abs + m_eps_rel * max BOOST_PREVENT_MACRO_SUBSTITUTION ( abs( x_old ) , abs( x ) ) );
            return r + tmp * tmp;
        }
    };




    template< class Fac1 = double >
    struct rel_error_l2_2
    {
        const Fac1 m_eps_abs , m_eps_rel , m_a_x , m_a_dxdt;

        rel_error_l2_2( Fac1 eps_abs , Fac1 eps_rel , Fac1 a_x , Fac1 a_dxdt )
        : m_eps_abs( eps_abs ) , m_eps_rel( eps_rel ) , m_a_x( a_x ) , m_a_dxdt( a_dxdt )
        { }

        template< class Res , class T1 , class T2 , class T3 , class T4 >
__device__ __host__        Res operator()( Res r , const T1 &x_old , const T2 &/*x*/ , const T3 &dxdt_old , const T4 &x_err )
        {
            using std::abs;
            Res tmp = abs( get_unit_value( x_err ) ) /
                    ( m_eps_abs + m_eps_rel * ( m_a_x * abs( get_unit_value( x_old ) ) + m_a_dxdt * abs( get_unit_value( dxdt_old ) ) ) );
            return r + tmp * tmp;
        }
    };


};


    } // odeint
} // numeric
} // boost


#endif //OPENFPM_PDATA_VECTOR_ALGEBRA_OFP_HPP
