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


template<int spatial_dim, typename MatType = EMatrixXd, typename VecType = EVectorXd>
class PolyLevelset
{
    minter::LevelsetPoly<spatial_dim, MatType, VecType> *model;
    
public:
    template<typename vector_type>
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

    // T : Point<vector_type::dims, typename vector_type::stype>
    template<typename T>
    double eval(T pos)
    {
        int dim = pos.dims;
        MatType point(1,dim);
        for(int j = 0;j < dim;++j)
            point(0,j) = pos.get(j);

        return model->eval(point)(0);
    }

    // T1 : Point<vector_type::dims, typename vector_type::stype>
    // T2 : Point<vector_type::dims, int>
    template<typename T1, typename T2>
    double deriv(T1 pos, T2 deriv_order)
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

    // T : Point<vector_type::dims, typename vector_type::stype>
    template<typename T>
    T estimate_normals_at(T pos)
    {
        int dim = pos.dims;
        MatType point(1,dim);
        for(int j = 0;j < dim;++j)
            point(0,j) = pos.get(j);

        T normal;
        auto normal_minter = model->estimate_normals_at(point);

        for(int j = 0;j < dim;++j)
            normal.get(j) = normal_minter(0,j);

        return normal;
    }

    // T : Point<vector_type::dims, typename vector_type::stype>
    template<typename T>
    double estimate_mean_curvature_at(T pos)
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