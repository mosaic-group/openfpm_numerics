//
// Created by Serhii
//

#ifndef OPENFPM_PDATA_DCPSERHS_CUH
#define OPENFPM_PDATA_DCPSERHS_CUH

#include "MonomialBasis.hpp"


template<unsigned int dim, typename MonomialBasis_type>
class DcpseRhs_gpu
{
private:
    const Point<dim, unsigned int>& differentialSignature;
    const MonomialBasis_type& monomialBasis;
    int sign;
public:
    __host__ __device__ DcpseRhs_gpu(const MonomialBasis_type &monomialBasis, const Point<dim, unsigned int> &differentialSignature);

    template<typename T>
    __host__ __device__ void getVector(T* b);
};


template<unsigned int dim, typename MonomialBasis_type>
__host__ __device__ DcpseRhs_gpu<dim, MonomialBasis_type>::DcpseRhs_gpu(const MonomialBasis_type &monomialBasis,
                        const Point<dim, unsigned int> &differentialSignature)
        : monomialBasis(monomialBasis), differentialSignature(differentialSignature) {
    unsigned int order = (Monomial_gpu<dim>(differentialSignature)).order();
    if (order % 2 == 0)
        sign = 1;
    else
        sign = -1;
}

template<unsigned int dim, typename MonomialBasis_type>
template<typename T>
__host__ __device__ void DcpseRhs_gpu<dim, MonomialBasis_type>::getVector(T* b)
{
    const auto& basisElements = monomialBasis.getElements();

    for (size_t i = 0; i < basisElements.size(); ++i) {
        const Monomial_gpu<dim>& dm = basisElements.get(i).getDerivative(differentialSignature);
        b[i] = sign * dm.evaluate(Point<dim, T>(0));
    }
}

#endif //OPENFPM_PDATA_DCPSERHS_CUH
