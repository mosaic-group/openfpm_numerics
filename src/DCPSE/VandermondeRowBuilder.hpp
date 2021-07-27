//
// Created by tommaso on 22/03/19.
// Modified by Serhii
//

#ifndef OPENFPM_PDATA_VANDERMONDEROW_HPP
#define OPENFPM_PDATA_VANDERMONDEROW_HPP

#include "MonomialBasis.hpp"

template <unsigned int dim, typename T, typename MonomialBasis_type = MonomialBasis<dim>>
class VandermondeRowBuilder
{
private:
    const MonomialBasis_type& monomialBasis;

public:
    __host__ __device__ VandermondeRowBuilder(const MonomialBasis_type &monomialBasis) : monomialBasis(monomialBasis) {}

    template <typename MatrixType>
    void buildRow(MatrixType &M, unsigned int row, Point<dim, T> x, T eps);

    __host__ __device__ void buildRow_gpu(T *M, unsigned int row, const T (& x) [dim], T eps);
};

template<unsigned int dim, typename T,typename MonomialBasis_type>
template <typename MatrixType>
void VandermondeRowBuilder<dim, T, MonomialBasis_type>::buildRow(MatrixType &M, unsigned int row, Point<dim, T> x, T eps)
{
    auto& basisElements = monomialBasis.getElements();

    for (size_t col = 0; col < basisElements.size(); ++col)
    {
        Monomial<dim> m = basisElements.get(col);
        M(row, col) = m.evaluate(x);
        M(row, col) /= openfpm::math::intpowlog(eps, m.order());
    }
}

template<unsigned int dim, typename T, typename MonomialBasis_type>
__host__ __device__ void VandermondeRowBuilder<dim, T, MonomialBasis_type>::buildRow_gpu(T* M, unsigned int row, const T (& x) [dim], T eps)
{
    // M is a pointer to the row being filled
    const auto& basisElements = monomialBasis.getElements();

    for (size_t col = 0; col < basisElements.size(); ++col)
    {
        const Monomial_gpu<dim>& m = basisElements.get(col);
        M[col] = m.evaluate(x) / pow(eps, m.order());
    }
}


#endif //OPENFPM_PDATA_VANDERMONDEROW_HPP
