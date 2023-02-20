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
#include "DCPSE/SupportBuilder.hpp"
#include "/Users/lschulze/Applications/openfpm_problems/hackathon/minter/include/minter.h"


template<typename vector_type_support, typename NN_type>
class RegressionSupport
{
	std::vector<size_t> keys;
	
public:

	template<typename iterator_type>
	RegressionSupport(vector_type_support &vd, iterator_type itPoint, unsigned int requiredSize, support_options opt, NN_type &cellList_in) : domain(vd),
	cellList(cellList_in)
	{
		rCut = cellList_in.getCellBox().getHigh(0);
		std::cout<<rCut<<std::endl;
		// Get spatial position from point iterator
		vect_dist_key_dx p = itPoint.get();
		vect_dist_key_dx pOrig = itPoint.getOrig();
		Point<vector_type_support::dims, typename vector_type_support::stype> pos = domain.getPos(p.getKey());

		// Get cell containing current point and add it to the set of cell keys
		grid_key_dx<vector_type_support::dims> curCellKey = cellList.getCellGrid(pos); // Here get the key of the cell where the current point is
		std::set<grid_key_dx<vector_type_support::dims>> supportCells;
        	supportCells.insert(curCellKey);

		// Make sure to consider a set of cells providing enough points for the support
		enlargeSetOfCellsUntilSize(supportCells, requiredSize + 1,
                                   opt); // NOTE: this +1 is because we then remove the point itself

        	// Now return all the points from the support into a vector
        	keys = getPointsInSetOfCells(supportCells, p, pOrig, requiredSize, opt);
	}
	
	auto getKeys()
	{
		return keys;
	}

	auto getNumParticles()
	{
		return keys.size();
	}

private:

	vector_type_support &domain;
	NN_type &cellList;
	typename vector_type_support::stype rCut;


    	void enlargeSetOfCellsUntilSize(std::set<grid_key_dx<vector_type_support::dims>> &set, unsigned int requiredSize,
                                    support_options opt) {
        	if (opt == support_options::RADIUS) {
            	auto cell = *set.begin();
            	grid_key_dx<vector_type_support::dims> middle;
            	int n = std::ceil(rCut / cellList.getCellBox().getHigh(0));
            	size_t sz[vector_type_support::dims];
            	for (int i = 0; i < vector_type_support::dims; i++) {
                	sz[i] = 2 * n + 1;
                	middle.set_d(i, n);
            	}
            	grid_sm<vector_type_support::dims, void> g(sz);
            	grid_key_dx_iterator<vector_type_support::dims> g_k(g);
            	while (g_k.isNext()) {
                	auto key = g_k.get();
                	key = cell + key - middle;
                	if (isCellKeyInBounds(key)) {
                    		set.insert(key);
                	}
                	++g_k;
            	}
        	}
	        else if (opt == support_options::AT_LEAST_N_PARTICLES) {

		int done = 0;
		int n = 0;
            	auto cell = *set.begin();
            	grid_key_dx<vector_type_support::dims> middle;
		size_t sz[vector_type_support::dims];
		
		while(true) // loop for the number of cells enlarged per dimension
		{
			std::set<grid_key_dx<vector_type_support::dims>> temp_set;
            		for (int i = 0; i < vector_type_support::dims; i++) {
                		sz[i] = 2 * n + 1;
                		middle.set_d(i, n);
            		}

            		grid_sm<vector_type_support::dims, void> g(sz);
            		grid_key_dx_iterator<vector_type_support::dims> g_k(g);
            		while (g_k.isNext()) {
                		auto key = g_k.get();
                		key = cell + key - middle;
                		if (isCellKeyInBounds(key)) {
                    			temp_set.insert(key);
                		}
                		++g_k;
            		}
			if (getNumElementsInSetOfCells(temp_set) < requiredSize) n++;
			else 
			{
				set = temp_set;
				std::cout<<"Enlarged "<<n<<" times"<<std::endl;
				break;
			}
		}
		}
		else {
            	while (getNumElementsInSetOfCells(set) <
                   5.0 * requiredSize) //Why 5*requiredSize? Becasue it can help with adaptive resolutions.
            	{
                	auto tmpSet = set;
                	for (const auto el: tmpSet) {
                    	for (unsigned int i = 0; i < vector_type_support::dims; ++i) {
                        	const auto pOneEl = el.move(i, +1);
                        	const auto mOneEl = el.move(i, -1);
                        	if (isCellKeyInBounds(pOneEl)) {
                            	set.insert(pOneEl);
                        	}
                        	if (isCellKeyInBounds(mOneEl)) {
                            	set.insert(mOneEl);
                        	}
                    	}
                	}

            	}	
		}   	
    }
    std::vector<size_t> getPointsInSetOfCells(std::set<grid_key_dx<vector_type_support::dims>> set,
                                              vect_dist_key_dx &p,
                                              vect_dist_key_dx &pOrig,
                                              size_t requiredSupportSize,
                                              support_options opt) {
        struct reord {
            typename vector_type_support::stype dist;
            size_t offset;

            bool operator<(const reord &p) const { return this->dist < p.dist; }
        };

        openfpm::vector<reord> rp;
        std::vector<size_t> points;
        Point<vector_type_support::dims, typename vector_type_support::stype> xp = domain.getPos(p);
        for (const auto cellKey: set) {
            const size_t cellLinId = getCellLinId(cellKey);
            const size_t elemsInCell = getNumElementsInCell(cellKey);
            for (size_t k = 0; k < elemsInCell; ++k) {
                size_t el = cellList.get(cellLinId, k);

                Point<vector_type_support::dims, typename vector_type_support::stype> xq = domain.getPosOrig(el);
                //points.push_back(el);

                reord pr;

                pr.dist = xp.distance(xq);
                pr.offset = el;
                rp.add(pr);
            }
        }

        if (opt == support_options::RADIUS) {
            for (int i = 0; i < rp.size(); i++) {
                if (rp.get(i).dist < rCut) {
                    points.push_back(rp.get(i).offset);
                }
            }
            /*      #ifdef SE_CLASS1
                    if (points.size()<requiredSupportSize)
                    {
                        std::cerr<<__FILE__<<":"<<__LINE__<<"Note that the DCPSE neighbourhood doesn't have asked no. particles (Increase the rCut or reduce the over_sampling factor)";
                        std::cout<<"Particels asked (minimum*oversampling_factor): "<<requiredSupportSize<<". Particles Possible with given options:"<<points.size()<<"."<<std::endl;
                    }
                    #endif*/
        }
	else if (opt == support_options::AT_LEAST_N_PARTICLES) {
	    for (int i = 0; i  < rp.size(); i++) points.push_back(rp.get(i).offset);
	    }
        else {
            rp.sort();
            for (int i = 0; i < requiredSupportSize; i++) {
                    points.push_back(rp.get(i).offset);
                }
            }
        //MinSpacing=MinSpacing/requiredSupportSize
        return points;
    }
    
    size_t getCellLinId(const grid_key_dx<vector_type_support::dims> &cellKey) {
        mem_id id = cellList.getGrid().LinId(cellKey);
        return static_cast<size_t>(id);
    }

    size_t getNumElementsInCell(const grid_key_dx<vector_type_support::dims> &cellKey) {
        const size_t curCellId = getCellLinId(cellKey);
        size_t numElements = cellList.getNelements(curCellId);
        return numElements;
    }
    size_t getNumElementsInSetOfCells(const std::set<grid_key_dx<vector_type_support::dims>> &set)
    {
	    size_t tot = 0;
	    for (const auto cell : set)
	    {
		    tot += getNumElementsInCell(cell);
	    }
	    return tot;
}

    bool isCellKeyInBounds(grid_key_dx<vector_type_support::dims> key)
    {
        const size_t *cellGridSize = cellList.getGrid().getSize();
        for (size_t i = 0; i < vector_type_support::dims; ++i)
        {
            if (key.value(i) < 0 || key.value(i) >= cellGridSize[i])
            {
                return false;
            }
        }
        return true;
    }
};


template<int spatial_dim, unsigned int prp_id, typename MatType = EMatrixXd, typename VecType = EVectorXd>
class RegressionModel
{

public:
	minter::PolyModel<spatial_dim, MatType, VecType> *model = nullptr;
	minter::PolyModel<spatial_dim, MatType, VecType> *deriv_model[spatial_dim];
	
	template<typename vector_type, typename dom_type>
	RegressionModel(vector_type &vd, dom_type &domain, unsigned int poly_degree, float lp_degree = 2.0)
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

		// construct polynomial model
    	model = new minter::PolyModel<spatial_dim, MatType, VecType>(points, values, poly_degree, lp_degree);

    	// placeholder for derivatives
		for(int i = 0;i < dim;++i)
			deriv_model[i] = nullptr;
	}

	
	// Constructor for all points in a proc (domain + ghost) and a specified poly_degree
	template<typename vector_type>
	RegressionModel(vector_type &vd, unsigned int poly_degree, float lp_degree = 2.0)
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

		// construct polynomial model 
		model = new minter::PolyModel<spatial_dim, MatType, VecType>(points, values, poly_degree, lp_degree);

		for(i = 0;i < dim;++i)
			deriv_model[i] = nullptr;
	}

	// Constructor for all points in a proc (domain + ghost) within a tolerance
	template<typename vector_type>
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
		minter::PolyModel<spatial_dim, MatType, VecType> *mdl = nullptr;
		
		do
		{
			++poly_degree;
			if(mdl)
				delete mdl;

			// construct polynomial model
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
	    	// std::cout<<"Fit of degree "<<poly_degree<<" with error = "<<error<<std::endl;

	    }while(error > tolerance);

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

	template<typename T> // Typical: Point<vector_type::dims, typename vector_type::stype> 
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

	void compute_grad()
	{
		for(int i = 0;i < spatial_dim;++i)
		{
			std::vector<int> ord(spatial_dim, 0);
			ord[i] = 1;
			deriv_model[i] = model->derivative(ord);
		}
	}

	// T: Point<vector_type::dims, typename vector_type::stype>
	template<typename T>
	T eval_grad(T pos)
	{
		T res;

		if(!deriv_model[0])
			compute_grad();

		for(int i = 0;i < spatial_dim;++i)
			res.get(i) = deriv_model[i]->eval(pos);

		return res;
	}

};


#endif /* REGRESSION_HPP_ */
