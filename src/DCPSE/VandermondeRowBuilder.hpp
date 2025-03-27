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
    void buildRow(MatrixType &M, unsigned int row, Point<dim, T> x, const openfpm::vector<T>& epsPowPrecomp);
};

template<unsigned int dim, typename T,typename MonomialBasis_type>
template <typename MatrixType>
void VandermondeRowBuilder<dim, T, MonomialBasis_type>::buildRow(MatrixType &M, unsigned int row, Point<dim, T> x, const openfpm::vector<T>& epsPowPrecomp)
{
    auto& basisElements = monomialBasis.getElements();

    for (size_t col = 0; col < basisElements.size(); ++col)
    {
        Monomial<dim> m = basisElements.get(col);
        M(row, col) = m.evaluate(x);
        M(row, col) /= epsPowPrecomp.get(col);
    }
}


#endif //OPENFPM_PDATA_VANDERMONDEROW_HPP
