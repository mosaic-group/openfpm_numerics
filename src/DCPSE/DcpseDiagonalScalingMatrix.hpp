//
// Created by tommaso on 29/03/19.
//

#ifndef OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
#define OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP

#include "MonomialBasis.hpp"
#include "Support.hpp"

template <unsigned int dim>
class DcpseDiagonalScalingMatrix
{
private:
    const MonomialBasis<dim> monomialBasis;

public:

    DcpseDiagonalScalingMatrix(const MonomialBasis<dim> &monomialBasis) : monomialBasis(monomialBasis) {}

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
        int i = 0;
        for (const auto& pt : support.getKeys())
        {
        	Point<dim,typename vector_type::stype> p = ref_p;
        	p -= particles.getPosOrig(pt);

            M(i,i) = exp(- norm2(p) / (2.0 * eps * eps));
            ++i;
        }
    }

};

#endif //OPENFPM_PDATA_DCPSEDIAGONALSCALINGMATRIX_HPP
