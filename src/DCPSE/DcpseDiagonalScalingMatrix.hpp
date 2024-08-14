//
// Created by tommaso on 29/03/19.
// Modified by Serhii
//

#ifndef OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
#define OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP

#include "MonomialBasis.hpp"


template <unsigned int dim, typename monomialBasis_type = MonomialBasis<dim>>
class DcpseDiagonalScalingMatrix
{
private:
    const monomialBasis_type& monomialBasis;

public:
    DcpseDiagonalScalingMatrix(const monomialBasis_type &monomialBasis) : monomialBasis(monomialBasis) {}

    template <typename T, typename MatrixType, typename verletIterator_type, typename vector_type, typename vector_type2>
    void buildMatrix(MatrixType &M, size_t p, verletIterator_type &it, T eps, vector_type & particlesSupport , vector_type2 & particlesDomain)
    {
        // Fill the diagonal matrix
        M.setZero(); // Make sure the rest of the matrix is zero!

        Point<dim,typename vector_type::stype> xp = particlesDomain.getPos(p);

        int i = 0;
        while (it.isNext())
        {
            Point<dim,typename vector_type::stype> _xp = xp;

            size_t q = it.get();
            _xp -= particlesSupport.getPos(q);

            M(i,i) = exp(- norm2(_xp) / (2.0 * eps * eps));

            ++it; ++i;
        }
    }

    template <typename T, typename vector_type, typename vector_type2>
    __host__ __device__ void buildMatrix(T* M, size_t supportRefKey, size_t supportKeysSize, const size_t* supportKeys, T eps, vector_type & particlesSupport, vector_type2 & particlesDomain)
    {
        // Check that all the dimension constraints are met
        assert(supportKeysSize >= monomialBasis.size());

        Point<dim,typename vector_type::stype> ref_p = particlesDomain.getPos(supportRefKey);

        for (size_t i = 0; i < supportKeysSize; ++i)
        {
            size_t pt = supportKeys[i];
            Point<dim,typename vector_type::stype> p = ref_p;
            p -= particlesSupport.getPos(pt);

            M[i] = exp(- norm2(p) / (2.0 * eps * eps));
        }
    }
};

#endif //OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
