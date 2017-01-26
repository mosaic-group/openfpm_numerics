/*
 * DrawParticles.hpp
 *
 *  Created on: Jan 5, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_HPP_
#define OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_HPP_

#include "PointIterator.hpp"
#include "PointIteratorSkin.hpp"
#include "Vector/vector_dist.hpp"

class DrawParticles
{
public:

	template<unsigned int dim, typename T, typename aggr, typename Decomposition> static PointIteratorSkin<dim,T,Decomposition>
	DrawSkin(vector_dist<dim,T,aggr,Decomposition> & vd,
			 size_t (& sz)[dim],
			 Box<dim,T> & domain,
			 Box<dim,T> & sub_A,
			 Box<dim,T> & sub_B)
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		return PointIteratorSkin<dim,T,Decomposition>(vd.getDecomposition(),sz,domain,sub_A, sub_B, bc);
	}

	template<unsigned int dim, typename T, typename aggr, typename Decomposition> static PointIteratorSkin<dim,T,Decomposition>
	DrawSkin(vector_dist<dim,T,aggr,Decomposition> & vd,
			 size_t (& sz)[dim],
			 Box<dim,T> & domain,
			 openfpm::vector<Box<dim,T>> & sub_A,
			 Box<dim,T> & sub_B)
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		PointIteratorSkin<dim,T,Decomposition> it(vd.getDecomposition(),sz,domain,sub_A.get(0), sub_B, bc);

		for (size_t i = 1 ; i < sub_A.size() ; i++)
			it.addBoxA(Box<dim,T>(sub_A.get(i)));

		return it;
	}

	template<unsigned int dim, typename T, typename aggr, typename Decomposition> static PointIterator<dim,T,Decomposition>
	DrawBox(vector_dist<dim,T,aggr,Decomposition> & vd,
			 size_t (& sz)[dim],
			 Box<dim,T> & domain,
			 Box<dim,T> & sub)
	{
		return PointIterator<dim,T,Decomposition>(vd.getDecomposition(),sz,domain,sub);
	}

};


#endif /* OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_HPP_ */
