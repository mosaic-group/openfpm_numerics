//
// Created by tommaso on 21/03/19.
// Edited by Abhinav Singh on 24/01/2022

#ifndef OPENFPM_PDATA_VANDERMONDE_HPP
#define OPENFPM_PDATA_VANDERMONDE_HPP

#include "MonomialBasis.hpp"
#include "VandermondeRowBuilder.hpp"
#include "Support.hpp"

template<unsigned int dim, typename T, typename MatrixType>
class Vandermonde
{
private:
    const Point<dim, T> point;
    openfpm::vector_std<Point<dim, T>> offsets;
    const MonomialBasis<dim> monomialBasis;
    T eps,HOverEpsilon;

public:
/*    Vandermonde(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours,
                const MonomialBasis<dim> &monomialBasis);*/

    template<typename vector_type,
             typename vector_type2>
    Vandermonde(const Support &support,
                const MonomialBasis<dim> &monomialBasis,
                const vector_type & particlesFrom,
                const vector_type2 & particlesTo,T HOverEpsilon=0.5)    //0.5 for the test
    : point(particlesTo.getPosOrig(support.getReferencePointKey())),
                  monomialBasis(monomialBasis),HOverEpsilon(HOverEpsilon)
    {
        initialize(support,particlesFrom,particlesTo);
    }


    MatrixType &getMatrix(MatrixType &M)
    {
        // Build the Vandermonde matrix, row-by-row
        VandermondeRowBuilder<dim, T> vrb(monomialBasis);
        unsigned int row = 0;

        size_t N = offsets.size();
        for (size_t i = 0; i < N; ++i)
        {
            const auto& offset = offsets.get(i);
            vrb.buildRow(M, row, offset, eps);
            ++row;
        }
        return M;
    }

    T getEps()
    {
        return eps;
    }

private:


    void computeEps(T factor)
    {
        T avgNeighbourSpacing = 0;
        size_t N = offsets.size();
        for (size_t i = 0; i < N; ++i)
        {
            const auto& offset = offsets.get(i);
            avgNeighbourSpacing += computeAbsSum(offset);
        }
        avgNeighbourSpacing /= offsets.size();
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

    template<typename vector_type, typename vector_type2>
    void initialize(const Support &sup, const vector_type & particlesFrom, vector_type2 &particlesTo)
    {
    	auto & keys = sup.getKeys();

    	for (int i = 0 ; i < keys.size() ; i++)
    	{
    		Point<dim,T> p = particlesTo.getPosOrig(sup.getReferencePointKey());
            p -= particlesFrom.getPosOrig(keys.get(i));
            offsets.add(p);
    	}

        // First check that the number of points given is enough for building the Vandermonde matrix
        if (offsets.size() < monomialBasis.size())
        {
            ACTION_ON_ERROR(std::length_error("Not enough neighbour points passed for Vandermonde matrix construction!"));
        }
        // Compute eps for this point
        //factor here. This is C factor.
        computeEps(HOverEpsilon);
    }

};



#endif //OPENFPM_PDATA_VANDERMONDE_HPP
