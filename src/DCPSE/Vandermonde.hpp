//
// Created by tommaso on 21/03/19.
// Edited by Abhinav Singh on 24/01/2022

#ifndef OPENFPM_PDATA_VANDERMONDE_HPP
#define OPENFPM_PDATA_VANDERMONDE_HPP

#include "MonomialBasis.hpp"
#include "VandermondeRowBuilder.hpp"

template<unsigned int dim, typename T, typename MatrixType>
class Vandermonde
{
private:
    const Point<dim, T> point;
    openfpm::vector_std<Point<dim, T>> x_pqDists;
    const MonomialBasis<dim> monomialBasis;
    T eps,HOverEpsilon,minSpacing;

public:
    template<typename verletIterator_type, typename vector_type, typename vector_type2>
    Vandermonde(
        size_t p,
        verletIterator_type &it,
        const MonomialBasis<dim> &monomialBasis,
        const vector_type & particlesSupport,
        const vector_type2 & particlesDomain,T HOverEpsilon=0.5
    ):
        monomialBasis(monomialBasis),
        HOverEpsilon(HOverEpsilon)
    {
        initialize(p, it, particlesSupport, particlesDomain);
    }


    MatrixType &getMatrix(MatrixType &M)
    {
        // Precompute eps to the power m as those don't change row-wise
        openfpm::vector<T> epsPowPrecomp(monomialBasis.size());
        auto& basisElements = monomialBasis.getElements();

        for (size_t i = 0; i < basisElements.size(); ++i)
        {
            Monomial<dim> m = basisElements.get(i);
            epsPowPrecomp.get(i) = openfpm::math::intpowlog(eps, m.order());
        }

        // Build the Vandermonde matrix, row-by-row
        VandermondeRowBuilder<dim, T> vrb(monomialBasis);
        unsigned int row = 0;

        size_t N = x_pqDists.size();
        for (size_t i = 0; i < N; ++i)
        {
            const auto& x_pq = x_pqDists.get(i);
            vrb.buildRow(M, row, x_pq, epsPowPrecomp);
            ++row;
        }
        return M;
    }

    T getEps()
    {
        return eps;
    }
    T getMinSpacing()
    {
        return minSpacing;
    }

private:


    void computeEps(T factor)
    {
        T avgNeighbourSpacing = 0;
        minSpacing=std::numeric_limits<T>::max();
        size_t N = x_pqDists.size();
        for (size_t i = 0; i < N; ++i)
        {
            const auto& x_pq = x_pqDists.get(i);
            T dist=norm(x_pq);
            avgNeighbourSpacing += computeAbsSum(x_pq);
            if(minSpacing>dist)
            {
                minSpacing=dist;
            }
        }
        avgNeighbourSpacing /= x_pqDists.size();
        eps = avgNeighbourSpacing/factor;
        assert(eps != 0);
    }

    static T computeAbsSum(const Point<dim, T> &x)
    {
        T absSum = 0;
        for (unsigned int i = 0; i < dim; ++i)
        {
            absSum += fabs(x.value(i));
        }
        return absSum;
    }

    template<typename verletIterator_type, typename vector_type, typename vector_type2>
    void initialize(size_t p, verletIterator_type &it, const vector_type & particlesSupport, vector_type2 &particlesDomain)
    {
        while (it.isNext())
        {
            size_t q = it.get();

            Point<dim,T> xp = particlesDomain.getPos(p);
            xp -= particlesSupport.getPos(q);
            x_pqDists.add(xp);

            ++it;
        }


        // First check that the number of points given is enough for building the Vandermonde matrix
        if (x_pqDists.size() < monomialBasis.size())
        {
            ACTION_ON_ERROR(std::length_error("Not enough neighbour points passed for Vandermonde matrix construction!"));
        }
        // Compute eps for this point
        // factor here. This is C factor.
        computeEps(HOverEpsilon);
    }

};



#endif //OPENFPM_PDATA_VANDERMONDE_HPP
