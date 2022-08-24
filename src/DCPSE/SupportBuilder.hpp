//
// Created by tommaso on 25/03/19.
//
// Modified by Abhinav and Pietro

#ifndef OPENFPM_PDATA_SUPPORTBUILDER_HPP
#define OPENFPM_PDATA_SUPPORTBUILDER_HPP

// This is to automatically get the right support (set of particles) for doing DCPSE on a given particle.
// todo: This could be enhanced to support an algorithm for choosing the support in a smart way (to lower condition num)

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>
#include "Support.hpp"
#include <utility>

enum support_options
{
    N_PARTICLES,
    RADIUS,
    LOAD,
    ADAPTIVE_SURFACE
};


template<typename vector_type,typename vector_type2>
class SupportBuilder
{
private:
    vector_type &domainFrom;
    vector_type2 &domainTo;
    decltype(std::declval<vector_type>().getCellList(0.0)) cellList;
    const Point<vector_type::dims, unsigned int> differentialSignature;
    typename vector_type::stype rCut,AvgSpacing;
    bool is_interpolation;

public:

    SupportBuilder(vector_type &domainFrom,vector_type2 &domainTo, const Point<vector_type::dims, unsigned int> differentialSignature,
                                                typename vector_type::stype rCut,
                                                bool is_interpolation)
    :domainFrom(domainFrom),
     domainTo(domainTo),
    differentialSignature(differentialSignature),
    rCut(rCut),is_interpolation(is_interpolation)
    {
        cellList = domainFrom.getCellList(rCut);
    }

    SupportBuilder(vector_type &domainFrom,vector_type2 &domainTo, unsigned int differentialSignature[vector_type::dims], typename vector_type::stype rCut, bool is_interpolation)
    : SupportBuilder(domainFrom,domainTo, Point<vector_type::dims, unsigned int>(differentialSignature), rCut) {}

    template<typename iterator_type>
    Support getSupport(iterator_type itPoint, unsigned int requiredSize, support_options opt)
    {
        // Get spatial position from point iterator
        vect_dist_key_dx p = itPoint.get();
        vect_dist_key_dx pOrig = itPoint.getOrig();
        Point<vector_type::dims, typename vector_type::stype> pos = domainTo.getPos(p.getKey());

        // Get cell containing current point and add it to the set of cell keys
        grid_key_dx<vector_type::dims> curCellKey = cellList.getCellGrid(pos); // Here get the key of the cell where the current point is
        std::set<grid_key_dx<vector_type::dims>> supportCells;
        supportCells.insert(curCellKey);

        // Make sure to consider a set of cells providing enough points for the support
        enlargeSetOfCellsUntilSize(supportCells, requiredSize + 1,opt); // NOTE: this +1 is because we then remove the point itself

        // Now return all the points from the support into a vector

        std::vector<size_t> supportKeys = getPointsInSetOfCells(supportCells,p,pOrig,requiredSize,opt);

        if (is_interpolation == false)
        {
            auto p_o = domainFrom.getOriginKey(p.getKey());
            std::remove(supportKeys.begin(), supportKeys.end(), p_o.getKey());
        }

        auto p_o = domainTo.getOriginKey(p.getKey());
        return Support(p_o.getKey(), openfpm::vector_std<size_t>(supportKeys.begin(), supportKeys.end()));
    }

    typename vector_type::stype getLastAvgspacing()
    {
        return this->AvgSpacing;
    }

private:

    size_t getCellLinId(const grid_key_dx<vector_type::dims> &cellKey)
    {
        mem_id id = cellList.getGrid().LinId(cellKey);
        return static_cast<size_t>(id);
    }

    size_t getNumElementsInCell(const grid_key_dx<vector_type::dims> &cellKey)
    {
        const size_t curCellId = getCellLinId(cellKey);
        size_t numElements = cellList.getNelements(curCellId);
        return numElements;
    }

    size_t getNumElementsInSetOfCells(const std::set<grid_key_dx<vector_type::dims>> &set)
    {
        size_t tot = 0;
        for (const auto cell : set)
        {
            tot += getNumElementsInCell(cell);
        }
        return tot;
    }

    void enlargeSetOfCellsUntilSize(std::set<grid_key_dx<vector_type::dims>> &set, unsigned int requiredSize,
            support_options opt)
    {
        if (opt==support_options::RADIUS){
            auto cell=*set.begin();
            grid_key_dx<vector_type::dims> middle;
            int n=std::ceil(rCut/cellList.getCellBox().getHigh(0));
            size_t sz[vector_type::dims];
            for (int i=0;i<vector_type::dims;i++)
            {
                sz[i]=2*n+1;
                middle.set_d(i,n);
            }
            grid_sm<vector_type::dims,void> g(sz);
            grid_key_dx_iterator<vector_type::dims> g_k(g);
            while(g_k.isNext())
            {
                auto key=g_k.get();
                key=cell+key-middle;
                if (isCellKeyInBounds(key))
                {
                    set.insert(key);
                }
                ++g_k;
            }
        }
        else{
            while (getNumElementsInSetOfCells(set) < 5.0*requiredSize) //Why 5*requiredSize? Becasue it can help with adaptive resolutions.
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
    }

    std::vector<size_t> getPointsInSetOfCells(std::set<grid_key_dx<vector_type::dims>> set,
                                                                            vect_dist_key_dx & p,
                                                                            vect_dist_key_dx & pOrig,
                                                                            size_t requiredSupportSize,
                                                                            support_options opt)
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
        Point<vector_type::dims,typename vector_type::stype> xp = domainTo.getPos(p);
        for (const auto cellKey : set)
        {
            const size_t cellLinId = getCellLinId(cellKey);
            const size_t elemsInCell = getNumElementsInCell(cellKey);
            for (size_t k = 0; k < elemsInCell; ++k)
            {
                size_t el = cellList.get(cellLinId, k);

                if (pOrig.getKey() == el && is_interpolation == false)   {continue;}

                Point<vector_type::dims,typename vector_type::stype> xq = domainFrom.getPosOrig(el);
                //points.push_back(el);

                reord pr;

                pr.dist = xp.distance(xq);
                pr.offset = el;

                rp.add(pr);
            }
        }

        if (opt == support_options::RADIUS)
        {
            for (int i = 0 ; i < rp.size() ; i++)
            {
                if (rp.get(i).dist < rCut)
                {
                    points.push_back(rp.get(i).offset);
                }
            }
    /*        #ifdef SE_CLASS1
            if (points.size()<requiredSupportSize)
            {
                std::cerr<<__FILE__<<":"<<__LINE__<<"Note that the DCPSE neighbourhood doesn't have asked no. particles (Increase the rCut or reduce the over_sampling factor)";
                std::cout<<"Particels asked (minimum*oversampling_factor): "<<requiredSupportSize<<". Particles Possible with given options:"<<points.size()<<"."<<std::endl;
            }
            #endif*/
        }
        else
        {   rp.sort();
            AvgSpacing=rp.get(1).dist;
            for (int i = 0 ; i < requiredSupportSize ; i++)
            {
                if(opt==support_options::ADAPTIVE_SURFACE && rp.get(i).dist > rCut)
                {}
                else{
                points.push_back(rp.get(i).offset);
                }
                //AvgSpacing+=rp.get(i).dist;
            }
            //AvgSpacing=AvgSpacing/requiredSupportSize;
        }

        return points;
    }

    bool isCellKeyInBounds(grid_key_dx<vector_type::dims> key)
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
};


#endif //OPENFPM_PDATA_SUPPORTBUILDER_HPP
