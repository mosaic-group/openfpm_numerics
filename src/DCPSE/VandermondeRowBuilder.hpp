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
    VandermondeRowBuilder(const MonomialBasis_type &monomialBasis) : monomialBasis(monomialBasis) {}

    template <typename MatrixType>
    void buildRow(MatrixType &V, unsigned int row, Point<dim, T> x);
};

template<unsigned int dim, typename T,typename MonomialBasis_type>
template <typename MatrixType>
void VandermondeRowBuilder<dim, T, MonomialBasis_type>::buildRow(MatrixType &V, unsigned int row, Point<dim, T> x)
{
    auto& basisElements = monomialBasis.getElements();

    for (size_t col = 0; col < basisElements.size(); ++col)
    {
        Monomial<dim> m = basisElements.get(col);
        V(row, col) = m.evaluate(x);
    }
}


#endif //OPENFPM_PDATA_VANDERMONDEROW_HPP
