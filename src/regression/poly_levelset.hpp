/*
 * Regression module
 * Obtains polynomial models for data from vector_dist
 * author : sachin (sthekke@mpi-cbg.de)
 * date : 18.01.2023
 *
 */
#ifndef OPENFPM_NUMERICS_POLYLEVELSET_HPP
#define OPENFPM_NUMERICS_POLYLEVELSET_HPP

#include "Vector/map_vector.hpp"
#include "Space/Shape/Point.hpp"
#include "DMatrix/EMatrix.hpp"

#include "minter.h"


template<int spatial_dim, typename vector_type, typename MatType = EMatrixXd, typename VecType = EVectorXd>
class PolyLevelset
{
    minter::LevelsetPoly<spatial_dim, MatType, VecType> *model;
    
public:
    PolyLevelset(vector_type &vd, double tol)
    {
        constexpr int dim = vector_type::dims;
    
        MatType points(vd.size_local(), dim);
        
        auto it = vd.getDomainIterator();
        int i = 0;
        while(it.isNext())
        {
            auto key = it.get();
            for(int j = 0;j < dim;++j)
                points(i,j) = vd.getPos(key)[j];

            ++i;
            ++it;
        }

        // construct polynomial model
        model = new minter::LevelsetPoly<spatial_dim, MatType, VecType>(points, tol);
    }


    ~PolyLevelset()
    {
        if(model)
            delete model;
    }

    // TODO: Make the return types more generic
    double eval(Point<vector_type::dims, typename vector_type::stype> pos)
    {
        int dim = pos.dims;
        MatType point(1,dim);
        for(int j = 0;j < dim;++j)
            point(0,j) = pos.get(j);

        return model->eval(point)(0);
    }

    double deriv(Point<vector_type::dims, typename vector_type::stype> pos, \
        Point<vector_type::dims, int> deriv_order)
    {
        int dim = pos.dims;
        MatType point(1,dim);
        for(int j = 0;j < dim;++j)
            point(0,j) = pos.get(j);

        std::vector<int> order;
        for(int j = 0;j < dim;++j)
            order.push_back(deriv_order.get(j));

        return model->deriv_eval(point, order)(0);
    }

    Point<vector_type::dims, typename vector_type::stype> estimate_normals_at(Point<vector_type::dims, typename vector_type::stype> pos)
    {
        int dim = pos.dims;
        MatType point(1,dim);
        for(int j = 0;j < dim;++j)
            point(0,j) = pos.get(j);

        Point<vector_type::dims, typename vector_type::stype> normal;
        auto normal_minter = model->estimate_normals_at(point);

        for(int j = 0;j < dim;++j)
            normal.get(j) = normal_minter(0,j);

        return normal;
    }

    double estimate_mean_curvature_at(Point<vector_type::dims, typename vector_type::stype> pos)
    {
        int dim = pos.dims;
        MatType point(1,dim);
        for(int j = 0;j < dim;++j)
            point(0,j) = pos.get(j);

        auto mc = model->estimate_mean_curvature_at(point);

        return mc(0);
    }
};



#endif /* POLYLEVELSET_HPP_ */