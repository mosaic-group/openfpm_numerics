//
// Created by Abhinav Singh on 03.11.20.
//


//#include "util/util_debug.hpp"
#include <boost/math/special_functions/spherical_harmonic.hpp>

//type used for dictionary arguments of a spherical harmonic coordinates
typedef std::tuple <int, int> lm;

/*! \brief Structure required for the Sph Harmonic amplitude dictionary arguments
 *
 */
struct key_hash : public std::unary_function<lm, std::size_t>
{
    std::size_t operator()(const lm& k) const
    {
        return std::get<0>(k) ^ std::get<1>(k);
    }
};

/*! \brief Structure required for the Sph Harmonic amplitude dictionary arguments
 *
 */
struct key_equal : public std::binary_function<lm, lm, bool>
{
    bool operator()(const lm& v0, const lm& v1) const
    {
        return (
                std::get<0>(v0) == std::get<0>(v1) &&
                std::get<1>(v0) == std::get<1>(v1)
        );
    }
};




namespace openfpm {
    namespace math {

        namespace detail {
            template<class T, class Policy>
            inline T spherical_harmonic_prefix_raw(unsigned n, unsigned m, T theta, const Policy &pol) {
                BOOST_MATH_STD_USING

                if (m > n)
                    return 0;
                T prefix;
                if (m==0){
                    prefix=1.0;
                }
                else{
                    prefix = boost::math::tgamma_delta_ratio(static_cast<T>(n - m + 1), static_cast<T>(2 * m), pol);
                }
                prefix *= (2 * n + 1) / (4 * boost::math::constants::pi<T>());
                prefix = sqrt(prefix);
                return prefix;
            }

            template<class T, class Policy>
            T Y(unsigned n, int m, T theta, T phi, const Policy &pol) {
                BOOST_MATH_STD_USING  // ADL of std functions

                bool sign = false;
                if (m < 0) {
                    // Reflect and adjust sign if m < 0:
                    sign = m & 1;
                    m = abs(m);
                }
               /* if (m & 1) {
                    // Check phase if theta is outside [0, PI]:
                    T mod = boost::math::tools::fmod_workaround(theta, T(2 * constants::pi<T>()));
                    if (mod < 0)
                        mod += 2 * constants::pi<T>();
                    if (mod > constants::pi<T>())
                        sign = !sign;
                }*/
                // Get the value and adjust sign as required:
                T prefix = spherical_harmonic_prefix_raw(n, m, theta, pol);
                //T sin_theta = sin(theta);
                T x = cos(theta);
                T leg = boost::math::legendre_p(n, m, x, pol);
                if (m != 0)
                    prefix *= sqrt(2);
                prefix *= leg;
                return sign ? prefix * sin(m * phi) : prefix * cos(m * phi);
            }

            template<class T, class Policy>
            T DYdTheta(unsigned n, int m, T theta, T phi, const Policy &pol) {
                BOOST_MATH_STD_USING  // ADL of std functions

                bool sign = false;
                if (m < 0) {
                    // Reflect and adjust sign if m < 0:
                    sign = m & 1;
                    m = abs(m);
                }
/*                if (m & 1) {
                    // Check phase if theta is outside [0, PI]:
                    T mod = boost::math::tools::fmod_workaround(theta, T(2 * constants::pi<T>()));
                    if (mod < 0)
                        mod += 2 * constants::pi<T>();
                    if (mod > constants::pi<T>())
                        sign = !sign;
                }*/
                // Get the value and adjust sign as required:
                T prefix = spherical_harmonic_prefix_raw(n, m, theta, pol);
                //T sin_theta = sin(theta);
                T x = cos(theta);
                T leg2,leg1 = boost::math::legendre_p(n, m + 1, x, pol);
                if(n+m==0||n-m==-1){
                    leg2=0.0;
                }
                else{
                leg2 = boost::math::legendre_p(n, m - 1, x, pol);
                }
                if (m != 0)
                    prefix *= sqrt(2);
                prefix *= (0.5 * leg1 - 0.5 * (n + m) * (n - m + 1) * leg2);
                return sign ? prefix * sin(m * phi) : prefix * cos(m * phi);
            }


            template<class T, class Policy>
            T DYdPhi(unsigned n, int m, T theta, T phi, const Policy &pol) {
                BOOST_MATH_STD_USING  // ADL of std functions

                bool sign = false;
                if (m == 0)
                    return 0;
                if (m < 0) {
                    // Reflect and adjust sign if m < 0:
                    sign = m & 1;
                    m = abs(m);
                }
               /* if (m & 1) {
                    // Check phase if theta is outside [0, PI]:
                    T mod = boost::math::tools::fmod_workaround(theta, T(2 * constants::pi<T>()));
                    if (mod < 0)
                        mod += 2 * constants::pi<T>();
                    if (mod > constants::pi<T>())
                        sign = !sign;
                }*/
                // Get the value and adjust sign as required:
                T prefix = spherical_harmonic_prefix_raw(n, m, theta, pol);
                T sin_theta = sin(theta);
                T x = cos(theta);
                T leg = boost::math::legendre_p(n, m, x, pol);
                prefix *= sqrt(2) * leg;
                return sign ? m * prefix * cos(m * phi) : -m * prefix * sin(m * phi);
            }

        }

    template<class T1, class T2, class Policy>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    Y(unsigned n, int m, T1 theta, T2 phi, const Policy &pol) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type result_type;
        typedef typename boost::math::policies::evaluation<result_type, Policy>::type value_type;
        return boost::math::policies::checked_narrowing_cast<result_type, Policy>(
                detail::Y(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi), pol),
                "openfpm::math::Y<%1%>(unsigned, int, %1%, %1%)");
    }

    template<class T1, class T2>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    Y(unsigned n, int m, T1 theta, T2 phi) {
        return openfpm::math::Y(n, m, theta, phi, boost::math::policies::policy<>());
    }


    template<class T1, class T2, class Policy>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    DYdTheta(unsigned n, int m, T1 theta, T2 phi, const Policy &pol) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type result_type;
        typedef typename boost::math::policies::evaluation<result_type, Policy>::type value_type;
        return boost::math::policies::checked_narrowing_cast<result_type, Policy>(
                detail::DYdTheta(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi), pol),
                "openfpm::math::DYdTheta<%1%>(unsigned, int, %1%, %1%)");
    }

    template<class T1, class T2>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    DYdTheta(unsigned n, int m, T1 theta, T2 phi) {
        return openfpm::math::DYdTheta(n, m, theta, phi, boost::math::policies::policy<>());
    }

    template<class T1, class T2, class Policy>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    DYdPhi(unsigned n, int m, T1 theta, T2 phi, const Policy &pol) {
        typedef typename boost::math::tools::promote_args<T1, T2>::type result_type;
        typedef typename boost::math::policies::evaluation<result_type, Policy>::type value_type;
        return boost::math::policies::checked_narrowing_cast<result_type, Policy>(
                detail::DYdPhi(n, m, static_cast<value_type>(theta), static_cast<value_type>(phi), pol),
                "openfpm::math::DYdPhi<%1%>(unsigned, int, %1%, %1%)");
    }

    template<class T1, class T2>
    inline typename boost::math::tools::promote_args<T1, T2>::type
    DYdPhi(unsigned n, int m, T1 theta, T2 phi) {
        return openfpm::math::DYdPhi(n, m, theta, phi, boost::math::policies::policy<>());
    }

    double sph_A1(int l,int  m,double v1, double vr) {
        return 0.5 * (1 + l) * (l * v1 - vr);
    }

    double sph_A2(int l,int  m,double v1, double vr) {
        return 0.5 * ((1 + l) * (-l) * v1 + (l + 3) * vr);
    }

    double sph_B(int l, int m,double v2) {
        return v2;
    }

    double sph_A3(int l,int m,double v1, double vr) {
        if (m == 1){
            return 0.5 *l* ((1 + l)*v1 - vr)-1.5*sph_A2(l,m,v1,vr);
        }
        else{
            return 0.5 *l* ((1 + l)*v1 - vr);
        }
    }

    double sph_A4(int l,int m,double v1, double vr) {
        if (m == 1){
            return 0.5* (-l*(1 + l)*v1 + (2-l)*vr)+0.5*sph_A2(l,m,v1,vr);
        }
        else{
            return 0.5* (-l*(1 + l)*v1 + (2-l)*vr);
        }
    }

    /*! \brief Returns Stokes Solution Amplitudes for Spherical Harmonics modes l,m
     *
     *  \param nu viscosity for nu*Lap(v)=Grad(P)
     *  \param l Spherical harmonic mode l
     *  \param m Spherical harmonic mode m
     *  \param vr Amplitude of Y^hat corresponding to the boundary
     *  \param v1 Amplitude of Psi^hat corresponding to the boundary
     *  \param v2 Amplitude of Phi^hat corresponding to the boundary
     *
     *  \return std::vector containing the spherical harmonic amplitudes (ur,u1,u2,p) for the solution at r for mode l,m.
     *
     */
    std::vector<double> sph_anasol_u(double nu,int l,int m,double vr,double v1,double v2,double r) {
         double ur,u1,u2,p;
         if(l==0)
             {
             if(r==0){
                 u2=sph_B(l,m,v2);
                 return {0,0,u2,0.0};
             }
             ur=sph_A1(l,m,v1,vr)*r+sph_A2(l,m,v1,vr)/r;
             u2=sph_B(l,m,v2);
             return {ur,0,u2,0};
             }
         BOOST_MATH_STD_USING
         //double ur,u1,u2,p;
         ur=sph_A1(l,m,v1,vr)*pow(r,l+1)+sph_A2(l,m,v1,vr)*pow(r,l-1);//+sph_A3(l,m,v1,vr)*pow(r,-l)+sph_A4(l,m,v1,vr)*pow(r,-l-2);
         //std::cout<<sph_A1(l,m,v1,vr)<<":"<<sph_A2(l,m,v1,vr)<<":"<<sph_A3(l,m,v1,vr)<<":"<<sph_A4(l,m,v1,vr)<<std::endl;
         //std::cout<<pow(r,l+1)<<":"<<pow(r,l-1)<<":"<<pow(r,-l)<<":"<<pow(r,-l-2)<<std::endl;
         u1=1.0/double(l*(l+1))*((l+3.0)*sph_A1(l,m,v1,vr)*pow(r,l+1)+(l+1)*sph_A2(l,m,v1,vr)*pow(r,l-1));//-(l-2)*sph_A3(l,m,v1,vr)*pow(r,-l)-l*sph_A4(l,m,v1,vr)*pow(r,-l-2);
         u2=sph_B(l,m,v2)*(pow(r,l));//+pow(r,-l-1));
         p= nu*((4.0*l+6.0)/double(l)*sph_A1(l,m,v1,vr)*pow(r,l));//(4.0*l-2.0)/double(l+1)*sph_A3(l,m,v1,vr)*pow(r,-l-1));
         return {ur,u1,u2,p};
        }

        /*! \brief Conversion from Vector spherical Harmonic basis to Cartesian basis
         *
         *  \param r Polar coordinate radius
         *  \param theta Polar coordinate theta between 0 and Pi
         *  \param phi Polar coordinate phi between -Pi and Pi
         *  \param Vr Dictionary of amplitudes of Y^hat with arguments l,m
         *  \param V1 Dictionary of amplitudes of Psi^hat with arguments l,m
         *  \param V2 Dictionary of amplitudes of Phi^hat with arguments l,m
         *  \return std::vector containing the cartesian coordinates corresponding to vector basis summation.
         *
         */
        template<unsigned int k>
        std::vector<double> sumY(double r, double theta, double phi,const std::unordered_map<const lm,double,key_hash,key_equal> &Vr,const std::unordered_map<const lm,double,key_hash,key_equal> &V1,const std::unordered_map<const lm,double,key_hash,key_equal> &V2) {
        double Sum1 = 0.0;
        double Sum2 = 0.0;
        double Sum3 = 0.0;
        for (int l = 0; l <= k; l++) {
            for (int m = -l; m <= l; m++) {
                auto Er= Vr.find(std::make_tuple(l,m));
                auto E1= V1.find(std::make_tuple(l,m));
                auto E2= V2.find(std::make_tuple(l,m));
                Sum1 += Er->second * openfpm::math::Y(l, m, theta, phi);
                double DYdPhi=openfpm::math::DYdPhi(l, m, theta, phi);
                double DYdTheta=openfpm::math::DYdTheta(l, m, theta, phi);
                /*if (DYdPhi==0 ||E2->second==0){
                    Sum2 += E1->second * DYdTheta;
                }
                else{*/
                //assert(theta!=0);
                Sum2 += E1->second * DYdTheta -
                                E2->second / sin(theta) * DYdPhi;
                /*}*/
                Sum3 += E2->second * DYdTheta +
                        E1->second/sin(theta) * DYdPhi;

           /*     Sum3 += E2->second *sin(theta)* DYdTheta +
                        E1->second * DYdPhi;*/

               /*if (DYdPhi==0 ||E1->second==0){
                        Sum3 += E2->second * DYdTheta;
                    }
                else{
                    Sum3 += E2->second * DYdTheta +
                            E1->second/sin(theta) * DYdPhi;
                }*/

            }
        }
        double x=Sum2*cos(theta)*cos(phi)-Sum3*sin(phi)+Sum1*cos(phi)*sin(theta);
        double y=Sum3*cos(phi)+Sum2*cos(theta)*sin(phi)+Sum1*sin(phi)*sin(theta);
        double z=Sum1*cos(theta)-Sum2*sin(theta);

        return {x,y,z};

    }
        /*! \brief Conversion from Scalar spherical Harmonic basis to Cartesian basis
         *
         *  \param r viscosity for nu*Lap(v)=Grad(P)
         *  \param theta Spherical harmonic mode l
         *  \param phi Spherical harmonic mode m
         *  \param Vr Dictionary of amplitudes of Y^hat with arguments l,m
         *  \return double containing the cartesian coordinate corresponding to scalar basis summation.
         *
         */
        template<unsigned int k>
        double sumY_Scalar(double r, double theta, double phi,const std::unordered_map<const lm,double,key_hash,key_equal> &Vr) {
            double Sum1 = 0.0;
            for (int l = 0; l <= k; l++) {
                for (int m = -l; m <= l; m++) {
                    auto Er= Vr.find(std::make_tuple(l,m));
                    Sum1 += Er->second * openfpm::math::Y(l, m, theta, phi);
                }
            }
            return Sum1;

        }

/*        double PsiTheta(unsigned l, int m, double theta, double phi) {
            return DYdTheta(l, m, theta, phi);
        }

        double PsiPhi(unsigned l, int m, double theta, double phi) {
            return DYdPhi(l, m, theta, phi);
        }
        double PhiTheta(unsigned l, int m, double theta, double phi) {
            return -1/sin(theta)*DYdPhi(l, m, theta, phi);
        }
        double PhiPhi(unsigned l, int m, double theta, double phi) {
            return sin(theta)*DYdTheta(l, m, theta, phi);
        }*/


}
}