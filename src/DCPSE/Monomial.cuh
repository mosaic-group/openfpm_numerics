//
// Created by Serhii
//

#ifndef OPENFPM_PDATA_MONOMIALBASISELEMENT_CUH
#define OPENFPM_PDATA_MONOMIALBASISELEMENT_CUH

#include "Space/Shape/Point.hpp"


template<unsigned int dim>
class Monomial_gpu
{
private:
    unsigned int sum = 0;
    unsigned int exponents[dim];
    unsigned int scalar = 1;

public:
    __host__ __device__ Monomial_gpu();
    __host__ __device__ Monomial_gpu(const Monomial_gpu<dim> &other);
    __host__ __device__ Monomial_gpu(const Monomial<dim> &other);
    __host__ __device__ explicit Monomial_gpu(const Point<dim, unsigned int> &other, unsigned int scalar = 1);
    __host__ __device__ explicit Monomial_gpu(const Point<dim, long int> &other, unsigned int scalar = 1);
    __host__ __device__ explicit Monomial_gpu(const unsigned int other[dim]);

    __host__ __device__ Monomial_gpu<dim> &operator=(const Monomial_gpu<dim> &other);
    __host__ __device__ Monomial_gpu<dim> &operator=(const Monomial<dim> &other);
    __host__ __device__ bool operator==(const Monomial_gpu<dim> &other) const;
    __host__ __device__ void swap(const Monomial_gpu<dim> &other);

    __host__ __device__ unsigned int order() const;
    __host__ __device__ unsigned int getExponent(unsigned int i) const;
    __host__ __device__ void setExponent(unsigned int i, unsigned int value);
    __host__ __device__ Monomial_gpu<dim> getDerivative(const Point<dim, unsigned int> differentialOrder) const;
    __host__ __device__ unsigned int getScalar() const { return scalar; }

    template<typename T> __host__ __device__ T evaluate(const Point<dim, T> x) const;
    template<typename T> __host__ __device__ T evaluate(const T (&x)[dim]) const;

private:
    __host__ __device__ void updateSum();
};

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim>::Monomial_gpu()
{
    for (size_t i = 0; i < dim; ++i) exponents[i] = 0;
    sum = 0;
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim>::Monomial_gpu(const Point<dim, unsigned int> &other, unsigned int scalar) : scalar(scalar)
{
    for (size_t i = 0; i < other.nvals; ++i)
        exponents[i] = other.value(i);
    updateSum();
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim>::Monomial_gpu(const Point<dim, long int> &other, unsigned int scalar) : scalar(scalar)
{
    for (size_t i = 0; i < other.nvals; ++i)
        exponents[i] = other.value(i);
    updateSum();
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim>::Monomial_gpu(const unsigned int other[dim]) : Monomial_gpu(Point<dim, unsigned int>(other))
{
    for (size_t i = 0; i < dim; ++i)
       exponents[i] = other[i];
    updateSum();
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim>::Monomial_gpu(const Monomial_gpu<dim> &other)
        : sum(other.sum), scalar(other.scalar) 
{
    for (size_t i = 0; i < dim; ++i)
       exponents[i] = other.exponents[i];
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim>::Monomial_gpu(const Monomial<dim> &other)
        : sum(other.order()), scalar(other.getScalar())
{
    for (size_t i = 0; i < dim; ++i)
       exponents[i] = other.getExponent(i);
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim> &Monomial_gpu<dim>::operator=(const Monomial_gpu<dim> &other)
{
    for (size_t i = 0; i < dim; ++i)
       exponents[i] = other.exponents[i];

    sum = other.sum;
    scalar = other.scalar;
    return *this;
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim> &Monomial_gpu<dim>::operator=(const Monomial<dim> &other)
{
    for (size_t i = 0; i < dim; ++i)
       exponents[i] = other.getExponent(i);

    sum = other.order();
    scalar = other.getScalar();
    return *this;
}

template<unsigned int dim>
__host__ __device__ void Monomial_gpu<dim>::updateSum()
{
    sum = 0;
    for (unsigned int i = 0; i < dim; ++i)
        sum += exponents[i];
}

template<unsigned int dim>
__host__ __device__ unsigned int Monomial_gpu<dim>::order() const
{
    return sum;
}

template<unsigned int dim>
__host__ __device__ unsigned int Monomial_gpu<dim>::getExponent(unsigned int i) const
{
    return exponents[i];
}

template<unsigned int dim>
__host__ __device__ void Monomial_gpu<dim>::setExponent(unsigned int i, unsigned int value)
{
    exponents[i] = value;
    updateSum();
}

template<unsigned int dim>
__host__ __device__ bool Monomial_gpu<dim>::operator==
        (const Monomial_gpu<dim> &other) const
{
    bool EQ = true;

    for (size_t i = 0; i < dim; ++i)
        if (exponents[i] != other[i])
            EQ = false;

    return EQ && (scalar == other.scalar);
}

template<unsigned int dim>
template<typename T>
__host__ __device__ T Monomial_gpu<dim>::evaluate(const Point<dim, T> x) const
{
    T res = scalar;
    for (unsigned int i = 0; i < dim; ++i)
        res *= pow(x[i], getExponent(i));

    return res;
}

template<unsigned int dim>
template<typename T>
__host__ __device__ T Monomial_gpu<dim>::evaluate(const T (& x) [dim]) const
{
    T res = scalar;
    for (unsigned int i = 0; i < dim; ++i)
        res *= pow(x[i], getExponent(i));

    return res;
}

template<unsigned int dim>
__host__ __device__ Monomial_gpu<dim> Monomial_gpu<dim>::getDerivative(const Point<dim, unsigned int> differentialOrder) const
{
    unsigned int s = scalar;
    Point<dim, unsigned int> e(exponents);
    for (unsigned int i = 0; i < dim; ++i)
    {
        unsigned int origExp = e.value(i);
        int targetExp = static_cast<int>(origExp) - static_cast<int>(differentialOrder.value(i));
        for (int k = origExp; k > targetExp && k >= 0; --k)
        {
            s *= k;
        }
        e.get(i) = static_cast<unsigned int>((targetExp < 0) ? 0 : targetExp);
    }
    return Monomial_gpu(e, s);
}

template<unsigned int dim>
__host__ __device__ void Monomial_gpu<dim>::swap(const Monomial_gpu<dim> &other)
{
    sum = other.sum;
    scalar = other.scalar;
    for (size_t i = 0; i < dim; ++i)
       exponents[i] = other.exponents[i];
}


#endif //OPENFPM_PDATA_MONOMIALBASISELEMENT_CUH
