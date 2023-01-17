/*
 * Regression module
 * Obtains polynomial models for data from vector_dist
 * author : sachin (sthekke@mpi-cbg.de)
 * date : 28.04.2022
 *
 */
#ifndef OPENFPM_NUMERICS_REGRESSION_HPP
#define OPENFPM_NUMERICS_REGRESSION_HPP

#include "Vector/map_vector.hpp"
#include "Space/Shape/Point.hpp"
#include "DMatrix/EMatrix.hpp"

#include "minter.h"

template<typename vector_type, typename NN_type>
class RegressionDomain
{
	openfpm::vector_std<size_t> keys;
	
public:

	RegressionDomain(vector_type &vd, Point<vector_type::dims, typename vector_type::stype> pos, typename vector_type::stype size, NN_type &NN)
	{
		// Not efficient to compute cell list each time.
        //auto NN = vd.getCellList(size);
        auto Np = NN.template getNNIterator(NN.getCell(pos));
     	
        while(Np.isNext())
        {
        	auto key = Np.get();
        	if (pos.distance(vd.getPos(key)) < size)
            	keys.add(key);
            ++Np;
        }
	}

	auto getKeys()
	{
		return keys;
	}

	auto getNumParticles()
	{
		return keys.size();
	}
};

// template<int spatial_dim, typename MatType, typename VecType>
// class OpenFPMPolyModel
// {

// public:
// 	// minter::PolyModel<spatial_dim, typename MatType::BaseMatrix, typename VecType::BaseMatrix> *model = nullptr;
// 	minter::PolyModel<spatial_dim, MatType, VecType> *model = nullptr;

// 	// OpenFPMPolyModel(minter::PolyModel<spatial_dim, typename MatType::BaseMatrix, typename VecType::BaseMatrix> *mdl) : model(mdl) {}
// 	OpenFPMPolyModel(minter::PolyModel<spatial_dim, MatType, VecType> *mdl) : model(mdl) {}

// 	// TODO: Make the return types more generic
// 	double eval(Point<spatial_dim, double> pos)
// 	{
// 		int dim = pos.dims;
// 		MatType point(1,dim);
// 		for(int j = 0;j < dim;++j)
// 			point(0,j) = pos.get(j);

// 		return model->eval(point)(0);
// 	}

// 	double deriv(Point<spatial_dim, double> pos, \
// 		Point<spatial_dim, int> deriv_order)
// 	{
// 		int dim = pos.dims;
// 		MatType point(1,dim);
// 		for(int j = 0;j < dim;++j)
// 			point(0,j) = pos.get(j);

// 		std::vector<int> order;
// 		for(int j = 0;j < dim;++j)
// 			order.push_back(deriv_order.get(j));

// 		return model->deriv_eval(point, order)(0);
// 	}

// };

template<int spatial_dim, unsigned int prp_id, typename vector_type, typename MatType, typename VecType>
class RegressionModel
{

public:
	minter::PolyModel<spatial_dim, MatType, VecType> *model = nullptr;
	minter::PolyModel<spatial_dim, MatType, VecType> *deriv_model[spatial_dim];
	
	template<typename dom_type>
	RegressionModel(vector_type &vd, dom_type &domain, unsigned int poly_degree, float lp_degree)
	{
		int num_particles = domain->getNumParticles();
		int dim = vector_type::dims;

		MatType points(num_particles, dim);
		VecType values(num_particles);
		
		auto keys = domain->getKeys();
		for(int i = 0;i < num_particles;++i)
		{
			for(int j = 0;j < dim;++j)
				points(i,j) = vd.getPos(keys.get(i))[j];
			values(i) = vd.template getProp<prp_id>(keys.get(i));
		}

		// std::cout << boost::core::demangle(typeid(decltype(obj)).name()) << '\n';

		// construct polynomial model (degree 4)
    	// auto mdl = new minter::PolyModel<spatial_dim, typename MatType::BaseMatrix, typename VecType::BaseMatrix>(points, values, poly_degree, lp_degree);
    	model = new minter::PolyModel<spatial_dim, MatType, VecType>(points, values, poly_degree, lp_degree);

		// model = new minter::PolyModel<spatial_dim, MatType, VecType>(mdl);

		for(int i = 0;i < dim;++i)
			deriv_model[i] = nullptr;
	}

	
	// Constructor for all points in a proc (domain + ghost) and a specified poly_degree
	RegressionModel(vector_type &vd, unsigned int poly_degree)
	{
		int num_particles = vd.size_local_with_ghost();
		int dim = vector_type::dims;

		MatType points(num_particles, dim);
		VecType values(num_particles);
		
		auto it = vd.getDomainAndGhostIterator();
		int i = 0;
		while (it.isNext())
		{
			auto key = it.get();
			for(int j = 0;j < dim;++j)
				points(i,j) = vd.getPos(key)[j];

			values(i) = vd.template getProp<prp_id>(key);

			++it;
			++i;
		}

		// construct polynomial model (degree 4)
    	//auto mdl = new minter::PolyModel<spatial_dim, typename MatType::BaseMatrix, typename VecType::BaseMatrix>(points, values, poly_degree, 2.0);
		// auto mdl = new minter::PolyModel<spatial_dim, MatType, VecType>(points, values, poly_degree, 2.0);
		
		model = new minter::PolyModel<spatial_dim, MatType, VecType>(points, values, poly_degree, 2.0);

		for(i = 0;i < dim;++i)
			deriv_model[i] = nullptr;
	}

	// Constructor for all points in a proc (domain + ghost) within a tolerance
	RegressionModel(vector_type &vd, double tolerance)
	{
		int num_particles = vd.size_local_with_ghost();
		int dim = vector_type::dims;

		MatType points(num_particles, dim);
		VecType values(num_particles);
		
		auto it = vd.getDomainAndGhostIterator();
		int i = 0;
		while (it.isNext())
		{
			auto key = it.get();
			for(int j = 0;j < dim;++j)
				points(i,j) = vd.getPos(key)[j];

			values(i) = vd.template getProp<prp_id>(key);

			++it;
			++i;
		}

		int poly_degree = 1;
		double error = -1.0;
		// minter::PolyModel<spatial_dim, typename MatType::BaseMatrix, typename VecType::BaseMatrix> *mdl = nullptr;
		minter::PolyModel<spatial_dim, MatType, VecType> *mdl = nullptr;
		
		do
		{
			++poly_degree;
			if(mdl)
				delete mdl;

			// construct polynomial model
	    	// mdl = new minter::PolyModel<spatial_dim, typename MatType::BaseMatrix, typename VecType::BaseMatrix>(points, values, poly_degree, 2.0);
	    	mdl = new minter::PolyModel<spatial_dim, MatType, VecType>(points, values, poly_degree, 2.0);

	    	// check if linf_error is within the tolerance
	    	error = -1.0;
	    	for(i = 0;i < num_particles;++i)
	    	{
	    		double pred = mdl->eval(points.block(i,0,1,dim))(0); // evaluated for one point
	    		double err = std::abs(pred - values(i));
	    		if (err > error)
	    			error = err;
	    	}
	    	std::cout<<"Fit of degree "<<poly_degree<<" with error = "<<error<<std::endl;

	    }while(error > tolerance);

	    // model = new OpenFPMPolyModel<spatial_dim, MatType, VecType>(mdl);
	    model = mdl;
	    for(i = 0;i < dim;++i)
			deriv_model[i] = nullptr;
	    
	}

	~RegressionModel()
	{

		if(model)
			delete model;

		for(int i = 0;i < spatial_dim;++i)
		{
			if(deriv_model[i])
				delete deriv_model[i];
		}
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

	void compute_grad()
	{
		for(int i = 0;i < spatial_dim;++i)
		{
			std::vector<int> ord(spatial_dim, 0);
			ord[i] = 1;
			// deriv_model[i] = new OpenFPMPolyModel<spatial_dim, MatType, VecType>(model->model->derivative(ord));
			deriv_model[i] = model->derivative(ord);
		}
	}

	Point<vector_type::dims, typename vector_type::stype> eval_grad(Point<vector_type::dims, typename vector_type::stype> pos)
	{
		Point<vector_type::dims, typename vector_type::stype> res;

		// int dim = pos.dims;
		// typename MatType::BaseMatrix point(1,dim);
		// for(int j = 0;j < dim;++j)
		// 	point(0,j) = pos.get(j);

		if(!deriv_model[0])
			compute_grad();

		for(int i = 0;i < spatial_dim;++i)
			res.get(i) = deriv_model[i]->eval(pos);

		return res;
	}
};


template<int spatial_dim, typename vector_type, typename MatType, typename VecType>
class LevelsetPoly
{
	minter::LevelsetPoly<spatial_dim, MatType, VecType> *model;
	
public:
	LevelsetPoly(vector_type &vd, double tol)
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

		// construct polynomial model (degree 4)
    	model = new minter::LevelsetPoly<spatial_dim, MatType, VecType>(points, tol);
	}


	~LevelsetPoly()
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

	double estimate_mean_curvatures_at(Point<vector_type::dims, typename vector_type::stype> pos)
	{
		int dim = pos.dims;
		MatType point(1,dim);
		for(int j = 0;j < dim;++j)
			point(0,j) = pos.get(j);

		auto mc = model->estimate_mean_curvature_at(point);

		return mc(0);
	}
};

#endif /* REGRESSION_HPP_ */