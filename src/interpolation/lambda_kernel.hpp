//
// Created by Abhinav Singh on 13.01.21.
//

#ifndef OPENFPM_PDATA_LAMBDAKERNEL_HPP
#define OPENFPM_PDATA_LAMBDAKERNEL_HPP

#include <iostream>

template<typename st>
double horner(const std::array<double,10> &v, st x)
{
    st s = 0;
    for(int i=9; i>=0; i--)
        s = v[i] + (s * x);
    return s;
}
//These needs to be Checked
constexpr std::array<double,10> c1={(1.0 * 12.0) , (0.0 * 12.0) , -(5.0 * 3.0) , (0.0 * 12.0) , (1.0 * 3.0) , -(100.0 * 4.0) , (455.0 * 3.0) , -(295.0 * 6.0) , (345.0 * 3.0) , -(115.0 * 2.0) };
constexpr std::array<double,10> c2={-(199.0 * 24.0) , (5485.0 * 6.0) , -(32975.0 * 3.0) , (28425.0 * 6.0) , -(61953.0 * 3.0) , (33175.0 * 4.0) , -(20685.0 * 3.0) , (3055.0 * 6.0) , -(1035.0 * 3.0) , (115.0 * 2.0) };
constexpr std::array<double,10> c3={(5913.0 * 24.0) , -(89235.0 * 6.0) , (297585.0 * 3.0) , -(143895.0 * 6.0) , (177871.0 * 3.0) , -(54641.0 * 4.0) , (19775.0 * 3.0) , -(1715.0 * 6.0) , (345.0 * 3.0) , -(23.0 * 2.0)};


template<typename st>
class lambda4_4kernel
{
public:
    static const int np = 6;
    static inline st value(st x, size_t i)
    {
        if (i == 0)
            return horner(c3, -x) / 24.0;
        else if (i == 1)
            return horner(c2, -x) / 24.0;
        else  if (i == 2)
            return horner(c1, -x) / 12.0;
        else if (i == 3)
            return horner(c1, x) / 12.0;
        else if (i == 4)
            return horner(c2, x) / 24.0;
        else if (i == 5)
            return horner(c3, x) / 24.0;
        return 0.0;
    }
};


template<typename st>
double horner22(const std::array<double,6> &v, st x)
{
    st s = 0;
    for(int i=5; i>=0; i--)
        s = v[i] + (s * x);
    return s;
}

constexpr std::array<double,6> c221={1.0,0.0,-1.0,-4.5,7.5,-3.0};
constexpr std::array<double,6> c222={-4.0,18.0,-29.0,21.5,-7.5,1.0};


template<typename st>
class lambda2_2kernel
{
public:
    static const int np = 6;
    static inline st value(st x, size_t i)
    {
        if (i == 0)
            return horner22(c221, -x);
        else if (i == 1)
            return horner22(c222, -x);
        return 0.0;
    }
};

#endif //OPENFPM_PDATA_LAMBDAKERNEL_HPP
