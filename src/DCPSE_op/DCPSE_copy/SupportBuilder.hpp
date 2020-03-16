//
// Created by tommaso on 25/03/19.
//

#ifndef OPENFPM_PDATA_SUPPORTBUILDER_HPP
#define OPENFPM_PDATA_SUPPORTBUILDER_HPP

// This is to automatically get the right support (set of particles) for doing DCPSE on a given particle.
// todo: This could be enhanced to support an algorithm for choosing the support in a smart way (to lower condition num)

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>
#include "Support.hpp"

template<typename vector_type>
class SupportBuilder
{
private:
    vector_type &domain;
    CellList<vector_type::dims, typename vector_type::stype, Mem_fast<HeapMemory, local_index_>> cellList;
    const Point<vector_type::dims, unsigned int> differentialSignature;

public:
    SupportBuilder(vector_type &domain, Point<vector_type::dims, unsigned int> differentialSignature, typename vector_type::stype rCut);

    SupportBuilder(vector_type &domain, unsigned int differentialSignature[vector_type::dims], typename vector_type::stype rCut);

    Support<vector_type::dims, typename vector_type::stype, typename vector_type::value_type> getSupport(vector_dist_iterator itPoint, unsigned int requiredSize);

private:
    size_t getCellLinId(const grid_key_dx<vector_type::dims> &cellKey);

    size_t getNumElementsInCell(const grid_key_dx<vector_type::dims> &cellKey);

    size_t getNumElementsInSetOfCells(const std::set<grid_key_dx<vector_type::dims>> &set);

    void enlargeSetOfCellsUntilSize(std::set<grid_key_dx<vector_type::dims>> &set, unsigned int requiredSize);

    std::vector<size_t> getPointsInSetOfCells(std::set<grid_key_dx<vector_type::dims>> set, vect_dist_key_dx & p, size_t requiredSupportSize);

    bool isCellKeyInBounds(grid_key_dx<vector_type::dims> key);
};

// Method definitions below

template<typename vector_type>
SupportBuilder<vector_type>::SupportBuilder(vector_type &domain, const Point<vector_type::dims, unsigned int> differentialSignature,
                                            typename vector_type::stype rCut) : domain(domain), differentialSignature(differentialSignature)
{
    cellList = domain.template getCellList<CellList<vector_type::dims, typename vector_type::stype, Mem_fast<HeapMemory, local_index_>>>(rCut);
}

template<typename vector_type>
Support<vector_type::dims, typename vector_type::stype, typename vector_type::value_type> SupportBuilder<vector_type>::getSupport(vector_dist_iterator itPoint, unsigned int requiredSize)
{
    // Get spatial position from point iterator
    vect_dist_key_dx p = itPoint.get();
    Point<vector_type::dims, typename vector_type::stype> pos = domain.getPos(p.getKey());

    // Get cell containing current point and add it to the set of cell keys
    grid_key_dx<vector_type::dims> curCellKey = cellList.getCellGrid(pos); // Here get the key of the cell where the current point is
    std::set<grid_key_dx<vector_type::dims>> supportCells;
    supportCells.insert(curCellKey);

    // Make sure to consider a set of cells providing enough points for the support
    enlargeSetOfCellsUntilSize(supportCells, requiredSize + 1); // NOTE: this +1 is because we then remove the point itself

    // Now return all the points from the support into a vector
    std::vector<size_t> supportKeys = getPointsInSetOfCells(supportCells,p,requiredSize);
    std::remove(supportKeys.begin(), supportKeys.end(), p.getKey());
    return Support<vector_type::dims, typename vector_type::stype, typename vector_type::value_type>(domain, p.getKey(), supportKeys);
}

template<typename vector_type>
size_t SupportBuilder<vector_type>::getNumElementsInCell(const grid_key_dx<vector_type::dims> &cellKey)
{
    const size_t curCellId = getCellLinId(cellKey);
    size_t numElements = cellList.getNelements(curCellId);
    return numElements;
}

template<typename vector_type>
size_t SupportBuilder<vector_type>::getNumElementsInSetOfCells(const std::set<grid_key_dx<vector_type::dims>> &set)
{
    size_t tot = 0;
    for (const auto cell : set)
    {
        tot += getNumElementsInCell(cell);
    }
    return tot;
}

template<typename vector_type>
void SupportBuilder<vector_type>::enlargeSetOfCellsUntilSize(std::set<grid_key_dx<vector_type::dims>> &set, unsigned int requiredSize)
{
    while (getNumElementsInSetOfCells(set) < 5*requiredSize)
    {
        auto tmpSet = set;
        for (const auto el : tmpSet)
        {
            for (unsigned int i = 0; i < vector_type::dims; ++i)
            {
                const auto pOneEl = el.move(i, +1);
                const auto mOneEl = el.move(i, -1);
                if (isCellKeyInBounds(pOneEl))
                {
                    set.insert(pOneEl);
                }
                if (isCellKeyInBounds(mOneEl))
                {
                    set.insert(mOneEl);
                }
            }
        }
    }
}

template<typename vector_type>
size_t SupportBuilder<vector_type>::getCellLinId(const grid_key_dx<vector_type::dims> &cellKey)
{
    mem_id id = cellList.getGrid().LinId(cellKey);
    return static_cast<size_t>(id);
}

template<typename vector_type>
std::vector<size_t> SupportBuilder<vector_type>::getPointsInSetOfCells(std::set<grid_key_dx<vector_type::dims>> set, vect_dist_key_dx & p, size_t requiredSupportSize)
{
    struct reord
    {
        typename vector_type::stype dist;
        size_t offset;

        bool operator<(const reord & p) const
        {return this->dist < p.dist;}
    };

    openfpm::vector<reord> rp;
    std::vector<size_t> points;
    Point<vector_type::dims,typename vector_type::stype> xp = domain.getPos(p);
    for (const auto cellKey : set)
    {
        const size_t cellLinId = getCellLinId(cellKey);
        const size_t elemsInCell = getNumElementsInCell(cellKey);
        for (size_t k = 0; k < elemsInCell; ++k)
        {
            size_t el = cellList.get(cellLinId, k);

            if (p.getKey() == el)   {continue;}

            Point<vector_type::dims,typename vector_type::stype> xq = domain.getPos(el);
            //points.push_back(el);

            reord pr;

            pr.dist = xp.distance(xq);
            pr.offset = el;

            rp.add(pr);
        }
    }

    rp.sort();

    for (int i = 0 ; i < requiredSupportSize ; i++)
    {
        points.push_back(rp.get(i).offset);
    }

    return points;
}

template<typename vector_type>
SupportBuilder<vector_type>::SupportBuilder(vector_type &domain, unsigned int *differentialSignature, typename vector_type::stype rCut)
        : SupportBuilder(domain, Point<vector_type::dims, unsigned int>(differentialSignature), rCut) {}

template<typename vector_type>
bool SupportBuilder<vector_type>::isCellKeyInBounds(grid_key_dx<vector_type::dims> key)
{
    const size_t *cellGridSize = cellList.getGrid().getSize();
    for (size_t i = 0; i < vector_type::dims; ++i)
    {
        if (key.value(i) < 0 || key.value(i) >= cellGridSize[i])
        {
            return false;
        }
    }
    return true;
}

#endif //OPENFPM_PDATA_SUPPORTBUILDER_HPP
