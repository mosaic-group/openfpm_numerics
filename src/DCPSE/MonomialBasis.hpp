//
// Created by tommaso on 20/03/19.
//

#ifndef OPENFPM_PDATA_MONOMIALBASIS_H
#define OPENFPM_PDATA_MONOMIALBASIS_H

#include "Vector/map_vector.hpp"
#include <Grid/grid_sm.hpp>
#include <Grid/iterators/grid_key_dx_iterator_sub_bc.hpp>
#include "Monomial.hpp"
#include "Monomial.cuh"


template<unsigned int dim, typename T = Monomial<dim>, template<typename, template<typename...> class...> class vector_type = openfpm::vector_std, template<typename...> class... Args>
class MonomialBasis
{
private:
    vector_type<T, Args...> basis;

public:
    MonomialBasis() {}

    MonomialBasis(const vector_type<unsigned int, Args...> &degrees, unsigned int convergenceOrder);

    MonomialBasis(unsigned int degrees[dim], unsigned int convergenceOrder);

//    explicit MonomialBasis(Point<dim, unsigned int> degrees, unsigned int convergenceOrder);

    __host__ __device__ explicit MonomialBasis(const vector_type<T, Args...> &basis) : basis(basis) {}

    __host__ __device__ MonomialBasis(const MonomialBasis &other);

    __host__ __device__ MonomialBasis &operator=(const MonomialBasis &other);

    __host__ __device__ unsigned int size() const;

    __host__ __device__ const T &getElement(size_t i) const;

    __host__ __device__ T &getElement(size_t i);

    __host__ __device__ const vector_type<T, Args...> &getElements() const;

    __host__ __device__ MonomialBasis<dim, T, vector_type, Args...> getDerivative(Point<dim, unsigned int> differentialOrder) const;

    __host__ __device__ bool operator==(const MonomialBasis &other) const;

    __host__ __device__ vector_type<T, Args...>& getBasis() { return basis; }

    template<typename charT, typename traits>
    friend std::basic_ostream<charT, traits> &
    operator<<(std::basic_ostream<charT, traits> &lhs, MonomialBasis<dim, T, vector_type, Args...> const &rhs)
    {
        lhs << "MonomialBasis: size=" << rhs.size() << ", elements={ ";
        for (const auto &el : rhs.getElements())
        {
            lhs << "(" << el << ") ";
        }
        lhs << "}" << std::endl;
        return lhs;
    }

private:
    void generateBasis(vector_type<unsigned int, Args...> m, unsigned int r);
};

//// Definitions below

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ MonomialBasis<dim, T, vector_type, Args...>::MonomialBasis(const vector_type<unsigned int, Args...> &degrees, unsigned int convergenceOrder)
{
    generateBasis(degrees, convergenceOrder);
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ MonomialBasis<dim, T, vector_type, Args...>::MonomialBasis(unsigned int *degrees, unsigned int convergenceOrder)
        : MonomialBasis(vector_type<unsigned int, Args...>(degrees, degrees + dim), convergenceOrder) {}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ MonomialBasis<dim, T, vector_type, Args...>::MonomialBasis(const MonomialBasis &other)
{
    basis = other.basis; // Here it works because both vector_type and Monomial perform a deep copy.
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ MonomialBasis<dim, T, vector_type, Args...> &MonomialBasis<dim, T, vector_type, Args...>::operator=(const MonomialBasis &other)
{
    basis = other.basis; // Here it works because both vector_type and Monomial perform a deep copy.
    return *this;
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ unsigned int MonomialBasis<dim, T, vector_type, Args...>::size() const
{
    return basis.size();
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ const T &MonomialBasis<dim, T, vector_type, Args...>::getElement(size_t i) const
{
    return basis.get(i);
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ T &MonomialBasis<dim, T, vector_type, Args...>::getElement(size_t i)
{
    return basis.get(i);
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
void MonomialBasis<dim, T, vector_type, Args...>::generateBasis(vector_type<unsigned int, Args...> m, unsigned int r)
{
    // Compute the vector of actual dimensions to iterate over
    // NOTE: each index can go up to sum(m)+r
    unsigned int mSum = 0U;
    for (size_t i = 0; i < m.size(); ++i) mSum += m.get(i);

    unsigned int orderLimit = mSum + r;
    size_t dimensions[dim];
    std::fill(dimensions, dimensions + dim, orderLimit);

    // Now initialize grid with appropriate size, then start-stop points and boundary conditions for the iterator
    grid_sm<dim, void> grid(dimensions);

    long int startV[dim] = {}; // 0-initialized
    grid_key_dx<dim, long int> start(startV);
    grid_key_dx<dim, long int> stop(dimensions);

    size_t bc[dim];
    std::fill(bc, bc + dim, NON_PERIODIC);

    grid_key_dx_iterator_sub_bc<dim> it(grid, start, stop, bc);

    // Finally compute alpha_min
    unsigned char alphaMin = static_cast<unsigned char>(!(mSum % 2)); // if mSum is even, alpha_min must be 1
    if(mSum==0)
    {
        alphaMin = 0;
    }
    //std::cout<<"AlphaMin: "<<alphaMin<<std::endl;
    //unsigned char alphaMin = 0; // we want to always have 1 in the basis

    while (it.isNext())
    {
        Point<dim, long int> p = it.get().get_k();
        T candidateBasisElement(p);
        // Filter out the elements which don't fullfil the theoretical condition for being in the vandermonde matrix
        if (candidateBasisElement.order() < orderLimit && candidateBasisElement.order() >= alphaMin)
        {
            basis.add(candidateBasisElement);
        }
        ++it;
    }
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ const vector_type<T, Args...> &MonomialBasis<dim, T, vector_type, Args...>::getElements() const
{
    return basis;
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ MonomialBasis<dim, T, vector_type, Args...> MonomialBasis<dim, T, vector_type, Args...>::getDerivative(const Point<dim, unsigned int> differentialOrder) const
{
    vector_type<T, Args...> derivatives;

    for (size_t i = 0; i < basis.size(); ++i)
    {
        // used insted of rhs ref as it does swap internally (not supported by Monomial)
        T d = basis.get(i).getDerivative(differentialOrder);
        derivatives.add(d);
    }

    return MonomialBasis<dim, T, vector_type, Args...>(derivatives);
}

template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
__host__ __device__ bool MonomialBasis<dim, T, vector_type, Args...>::operator==(const MonomialBasis &other) const
{
    return basis == other.basis;
}

//template<unsigned int dim, typename T, template<typename, template<typename...> class...> class vector_type, template<typename...> class... Args>
// __host__ __device__ //MonomialBasis<dim, T, vector_type, Args...>::MonomialBasis(Point<dim, unsigned int> degrees, unsigned int convergenceOrder)
//        : MonomialBasis(degrees.asArray(), convergenceOrder) {}

#endif //OPENFPM_PDATA_MONOMIALBASIS_H
