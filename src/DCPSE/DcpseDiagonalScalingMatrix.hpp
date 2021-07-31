//
// Created by tommaso on 29/03/19.
// Modified by Serhii
//

#ifndef OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
#define OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP

#include "MonomialBasis.hpp"
#include "Support.hpp"


template <unsigned int dim, typename monomialBasis_type = MonomialBasis<dim>>
class DcpseDiagonalScalingMatrix
{
private:
    const monomialBasis_type& monomialBasis;

public:
    __host__ __device__ DcpseDiagonalScalingMatrix(const monomialBasis_type &monomialBasis) : monomialBasis(monomialBasis) {}

    template <typename T, typename MatrixType, typename vector_type>
    void buildMatrix(MatrixType &M, Support support, T eps, vector_type & particles)
    {
        // Check that all the dimension constraints are met
        assert(support.size() >= monomialBasis.size());
        assert(M.rows() == support.size());
        assert(M.cols() == support.size());

        Point<dim,typename vector_type::stype> ref_p = particles.getPosOrig(support.getReferencePointKey());

        // Fill the diagonal matrix
        M.setZero(); // Make sure the rest of the matrix is zero!
        const auto& support_keys = support.getKeys();
        size_t N = support_keys.size();
        for (size_t i = 0; i < N; ++i)
        {
            const auto& pt = support_keys.get(i);
        	Point<dim,typename vector_type::stype> p = ref_p;
        	p -= particles.getPosOrig(pt);

            M(i,i) = exp(- norm2(p) / (2.0 * eps * eps));
        }
    }

    template <typename T, typename vector_type>
    __host__ __device__ void buildMatrix(T* M, size_t supportRefKey, size_t supportKeysSize, const size_t* supportKeys, T eps, vector_type & particles)
    {
        // Check that all the dimension constraints are met
        assert(supportKeysSize >= monomialBasis.size());

        Point<dim,typename vector_type::stype> ref_p = particles.getPos(supportRefKey);

        for (size_t i = 0; i < supportKeysSize; ++i)
        {
            size_t pt = supportKeys[i];
            Point<dim,typename vector_type::stype> p = ref_p;
            p -= particles.getPos(pt);

            M[i] = exp(- norm2(p) / (2.0 * eps * eps));
        }
    }
};

#endif //OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
