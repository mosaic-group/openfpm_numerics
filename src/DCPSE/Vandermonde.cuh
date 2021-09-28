//
// Created by Serhii on 15/05/21
// 


#ifndef OPENFPM_PDATA_VANDERMONDE_CUH
#define OPENFPM_PDATA_VANDERMONDE_CUH

#include "MonomialBasis.hpp"
#include "VandermondeRowBuilder.hpp"
#include "Support.hpp"

template<unsigned int dim, typename T, typename MonomialBasis_type = MonomialBasis<dim>>
class Vandermonde_gpu
{
private:
    size_t N;
    T* offsets;
    const Point<dim, T> point;
    const MonomialBasis_type& monomialBasis;
    T eps;

public:
    template<typename vector_type>
    __host__ __device__ Vandermonde_gpu( T* mem, size_t supportRefKey, size_t supportKeysSize,
        const size_t* supportKeys, const MonomialBasis_type &monomialBasis, const vector_type & particles)
            : offsets(mem), monomialBasis(monomialBasis), point(particles.getPos(supportRefKey))
    {
        initialize(supportRefKey, supportKeysSize, supportKeys, particles);
    }

    __host__ __device__ void getMatrix(T *M)
    {
        // Build the Vandermonde_gpu matrix, row-by-row
        VandermondeRowBuilder<dim, T, MonomialBasis_type> vrb(monomialBasis);
        unsigned int row = 0;

        for (size_t i = 0; i < N; ++i)
        {
            const T(*offset) [dim] = reinterpret_cast<const T(*) [dim]>(&offsets[i*dim]);
            vrb.buildRow_gpu(&M[i*monomialBasis.size()], row, *offset, eps);
            ++row;
        }
    }

    __host__ __device__ T getEps()
    {
        return eps;
    }


    __host__ __device__ ~Vandermonde_gpu()
    {
    }

private:


    __host__ __device__ void computeEps(T factor)
    {
        T avgNeighbourSpacing = 0;
        for (size_t i = 0; i < N; ++i)
        {
            const T(*offset) [dim] = reinterpret_cast<const T(*) [dim]>(&offsets[i*dim]);
            avgNeighbourSpacing += computeAbsSum(*offset);

        }
        avgNeighbourSpacing /= N;
        eps = factor * avgNeighbourSpacing;
        assert(eps != 0);
    }

    __host__ __device__ static T computeAbsSum(const Point<dim, T> &x)
    {
        T absSum = 0;
        for (unsigned int i = 0; i < dim; ++i)
        {
            absSum += fabs(x.value(i));
        }
        return absSum;
    }

    __host__ __device__ static T computeAbsSum(const T (& x) [dim])
    {
        T absSum = 0;
        for (unsigned int i = 0; i < dim; ++i)
        {
            absSum += fabs(x[i]);
        }
        return absSum;
    }


    template<typename vector_type>
    __host__ __device__ void initialize(size_t supportRefKey, size_t supportKeysSize, const size_t* supportKeys, const vector_type & particles)
    {
        N = supportKeysSize;

        for (int i = 0 ; i < N ; i++)
        {
            size_t otherKey = supportKeys[i];
            Point<dim,T> p = particles.getPos(supportRefKey);
            const auto& p2 = particles.getPos(otherKey);

            p -= p2;

            for (size_t j = 0; j < dim; j++) {
                offsets[i*dim+j] = p[j];
                auto diff = p[j];
            }
        }

        computeEps(2);
    }

};


#endif //OPENFPM_PDATA_VANDERMONDE_CUH
