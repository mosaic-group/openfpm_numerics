//
// Created by tommaso on 21/03/19.
//

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
    std::vector<Point<dim, T>> offsets;
    const MonomialBasis<dim> monomialBasis;
    T eps;

public:
/*    Vandermonde(const Point<dim, T> &point, const std::vector<Point<dim, T>> &neighbours,
                const MonomialBasis<dim> &monomialBasis);*/

    template<typename vector_type>
    Vandermonde(const Support &support,
                const MonomialBasis<dim> &monomialBasis,
                const vector_type & particles)
    : point(particles.getPosOrig(support.getReferencePointKey())),
                  monomialBasis(monomialBasis)
    {
        initialize(support,particles);
    }


    MatrixType &getMatrix(MatrixType &M)
    {
        // Build the Vandermonde matrix, row-by-row
        VandermondeRowBuilder<dim, T> vrb(monomialBasis);
        unsigned int row = 0;
        for (auto &offset : offsets)
        {
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
        for (auto &offset : offsets)
        {
            avgNeighbourSpacing += computeAbsSum(offset);
        }
        avgNeighbourSpacing /= offsets.size();
        eps = factor * avgNeighbourSpacing;
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

    template<typename vector_type>
    void initialize(const Support &sup, const vector_type & particles)
    {
    	auto & keys = sup.getKeys();

    	for (int i = 0 ; i < keys.size() ; i++)
    	{
    		Point<dim,T> p = particles.getPosOrig(sup.getReferencePointKey());
    		p -= particles.getPosOrig(keys[i]);
    		offsets.push_back(p);
    	}

        // First check that the number of points given is enough for building the Vandermonde matrix
        if (offsets.size() < monomialBasis.size())
        {
            ACTION_ON_ERROR(std::length_error("Not enough neighbour points passed for Vandermonde matrix construction!"));
        }
        // Compute eps for this point
        //factor here. This is C factor.
        computeEps(2);
    }

};



#endif //OPENFPM_PDATA_VANDERMONDE_HPP
