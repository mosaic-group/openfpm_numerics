/*
 * PointIteratorSkin.hpp
 *
 *  Created on: Jan 4, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_DRAW_POINTITERATORSKIN_HPP_
#define OPENFPM_NUMERICS_SRC_DRAW_POINTITERATORSKIN_HPP_

#include "Grid/Iterators/grid_dist_id_iterator_dec_skin.hpp"

#define RETURN_A 1
#define RETURN_B 2

/*! \brief this class draw particles on subset of grid-like position
 *
 \verbatim
             (Box B)
                |                    (0.6,0.6)
  +   +   +   + | +   +   +   +   +
    +---+-------+---+
  + | * | +   + | * | +   +   +   +
    |   |       |
  + | * | +   + | * | +   +   +   +
    |   |       |   |
  + | * | +   + | * | +   +   +   +
    |   |       |   |
  + | * | +   + | * | +   +   +   +
    |   |       |   |
  + | * | +   + | * | +   +   +   +
    |   |       |   |
  + | * | +   + | * | +   +   +   +
    |   +-------+   |
  + | *   *   *   * | +   +   +   +
    +---------------+
  +   +   +   +   + | +   +   +   +
(-1.0,-1.0)         |
                    |
                 (Box A)


 \endverbatim
 *
 *
 *  Suppose we have a grid 9x9 from (-1.0,-1.0 to 0.6,0.6)
 *
 *  Defined a Box A (-0.9,-0.9 to -0.1,0.5) and a box B (-0.7, -0.7 to -0.3,0.5)
 *
 *  This class will return the points indicated with *
 *
 */
template<unsigned int dim, typename T, typename Decomposition>
class PointIteratorSkin: protected grid_dist_id_iterator_dec_skin<Decomposition>
{
	//! Actual point
	Point<dim,T> ap;

	//! sub_domain (required to filter out points)
	openfpm::vector<Box<dim,T>> sub_domainA;

	//! domain
	Box<dim,T> domain;

	//! Spacing
	T sp[dim];

	static Box<dim,long int> getAB(size_t (& gs)[dim], const Box<dim,T> & dom, const Box<dim,T> & sub_domA , const Box<dim,T> & sub_domB, T (& sp)[dim], size_t AB)
	{
		for (size_t i = 0 ; i < dim ; i++)
			sp[i] = (dom.getHigh(i) - dom.getLow(i)) / (gs[i] -1);

		Box<dim,long int> box;

		for (size_t i = 0 ; i < dim ; i++)
		{
			size_t Ast = std::ceil( (sub_domA.getLow(i) - dom.getLow(i)) / sp[i] ) - 1;
			size_t Asp = std::floor( (sub_domA.getHigh(i) - dom.getLow(i)) / sp[i] ) + 1;

			size_t Bst = std::ceil( (sub_domB.getLow(i) - dom.getLow(i)) / sp[i] );
			size_t Bsp = std::floor( (sub_domB.getHigh(i) - dom.getLow(i)) / sp[i] );

			// grid_dist_id_iterator_dec_skin only work if A is contained into B
			Ast = (Ast < Bst)?Bst:Ast;
			Asp = (Asp > Bsp)?Bsp:Asp;

			if (AB == RETURN_A)
			{
				box.setLow(i,Ast);
				box.setHigh(i,Asp);
			}
			else
			{
				box.setLow(i,Bst);
				box.setHigh(i,Bsp);
			}
		}

		return box;
	}

	void calculateAp()
	{
		grid_key_dx<dim> key = grid_dist_id_iterator_dec_skin<Decomposition>::get();

		for (size_t i = 0 ; i < dim ; i++)
			ap.get(i) = key.get(i) * sp[i] + domain.getLow(i);
	}


	/*! it check that the actual point is not inside B
	 *
	 * \return true if the point is not inside B
	 *
	 */
	bool isValidPoint()
	{
		bool valid = true;

		for (size_t i = 0 ; i < sub_domainA.size() ; i++)
		{
			if (Box<dim,T>(sub_domainA.get(i)).isInside(ap) == true)
				valid = false;
		}

		return valid;
	}

public:

	/*! \brief Draw Particles
	 *
	 * \param sp grid spacing
	 *
	 */
	PointIteratorSkin( Decomposition & dec, size_t (& sz)[dim], Box<dim,T> & domain, const Box<dim,T> & sub_A, const Box<dim,T> & sub_B, size_t (& bc)[dim])
	:grid_dist_id_iterator_dec_skin<Decomposition>(dec, sz, getAB(sz,domain,sub_A,sub_B,sp,RETURN_A), getAB(sz,domain,sub_A,sub_B,sp,RETURN_B), bc),domain(domain)
	{
		sub_domainA.add(sub_A);
		calculateAp();
	}

	/*! \Return the actual point
	 *
	 *
	 */
	Point<dim,T> & get()
	{
		return ap;
	}

	/*! \brief Next point
	 *
	 * \return itself
	 *
	 */
	PointIteratorSkin & operator++()
	{
		grid_dist_id_iterator_dec_skin<Decomposition>::operator++();
		calculateAp();

		while (grid_dist_id_iterator_dec_skin<Decomposition>::isNext() && isValidPoint() == false)
		{
			grid_dist_id_iterator_dec_skin<Decomposition>::operator++();
			calculateAp();
		}

		return *this;
	}

	bool isNext()
	{
		return grid_dist_id_iterator_dec_skin<Decomposition>::isNext();
	}

	PointIteratorSkin<dim,T,Decomposition> & operator=(PointIteratorSkin<dim,T,Decomposition> & p)
	{
		grid_dist_id_iterator_dec_skin<Decomposition>::operator=(p);

		ap = p.ap;
		sub_domainA = p.sub_domainA;
		domain = p.domain;

		for (size_t i = 0 ; i < dim; i++)
			sp[i] = p.sp[i];

		return *this;
	}

	void addBoxA(const Box<dim,double> & BoxA)
	{
		sub_domainA.add(BoxA);
	}
};



#endif /* OPENFPM_NUMERICS_SRC_DRAW_POINTITERATORSKIN_HPP_ */
