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

/*! \brief A class to draw/create particles based on simple shaped
 *
 * ## Draw box example
 *
 * \snippet Draw/DrawParticles_unit_tests.hpp DrawBox_example
 *
 */
class DrawParticles
{
public:

	/*! \brief Draw particles in a box B excluding the area of a second box A (B - A)
	 *
	 * The function return an iterator over particles defined on a virtual grid in the simulation domain.
	 *
	 * \param vd particles where we are creating the particles
	 * \param sz indicate the grid size of the virtual grid.
	 * \param domain Domain where the virtual grid is defined (Must match the domain of vd)
	 * \param sub_B box contained in domain that define where the particle iterator must iterate,
	 *            particles are placed strictly inside this box
	 * \param sub_A box contained in the domain that define where the particle iterator should not
	 *              iterate (excluding area)
	 *
	 * \note Suppose to have a simulation domain of 1.5 on x and we use sz = 16. Consider now
	 * to have particles with spacing 0.1 on x. if we define a sub_A that on extend from 0.65
	 *  to 0.95 the first fluid particle
	 * is at 0.70 and the last is at 0.90
	 *
	 * \return an iterator to the selected particles
	 *
	 */
	template<unsigned int dim, typename T, typename aggr, typename layout, template<typename> class layout_base ,typename Decomposition> static PointIteratorSkin<dim,T,Decomposition>
	DrawSkin(vector_dist<dim,T,aggr,layout,layout_base,Decomposition> & vd,
			 size_t (& sz)[dim],
			 Box<dim,T> & domain,
			 Box<dim,T> & sub_A,
			 Box<dim,T> & sub_B)
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		return PointIteratorSkin<dim,T,Decomposition>(vd.getDecomposition(),sz,vd.getDecomposition().getDomain(),sub_A, sub_B, bc);
	}


	/*! \brief Draw particles in a box B excluding the areas of an array of boxes A_n
	 *
	 * The function return an iterator over particles defined on a virtual grid in the simulation domain.
	 *
	 * \param vd particles where we are creating the particles
	 * \param sz indicate the grid size of the virtual grid.
	 * \param domain Domain where the virtual grid is defined (Must match the domain of vd)
	 * \param sub_B box contained in domain that define where the particle iterator must iterate,
	 *            particles are placed strictly inside this box
	 * \param sub_A array of boxes contained in the domain that define where the particle iterator should not
	 *              iterate (excluding areas)
	 *
	 * \note Suppose to have a simulation domain of 1.5 on x and we use sz = 16. Consider now
	 * to have particles with spacing 0.1 on x. if we define a sub_A that on extend from 0.65
	 *  to 0.95 the first fluid particle
	 * is at 0.70 and the last is at 0.90
	 *
	 * \return an iterator to the selected particles
	 *
	 */
	template<unsigned int dim, typename T, typename aggr, typename layout, template <typename> class layout_base, typename Decomposition> static PointIteratorSkin<dim,T,Decomposition>
	DrawSkin(vector_dist<dim,T,aggr,layout,layout_base,Decomposition> & vd,
			 size_t (& sz)[dim],
			 Box<dim,T> & domain,
			 openfpm::vector<Box<dim,T>> & sub_A,
			 Box<dim,T> & sub_B)
	{
		size_t bc[dim];

		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		PointIteratorSkin<dim,T,Decomposition> it(vd.getDecomposition(),sz,vd.getDecomposition().getDomain(),sub_A.get(0), sub_B, bc);

		for (size_t i = 1 ; i < sub_A.size() ; i++)
			it.addBoxA(Box<dim,T>(sub_A.get(i)));

		return it;
	}

	/*! \brief Draw particles in a box
	 *
	 * The function return an iterator over particles defined on a virtual grid in the simulation domain.
	 *
	 * \param vd particles where we are creating the particles
	 * \param sz indicate the grid size of the virtual grid.
	 * \param domain Domain where the virtual grid is defined (Must match the domain of vd)
	 * \param sub box contained in domain that define where the particle iterator must iterate,
	 *            particles are placed strictly inside this box
	 *
	 * \note Suppose to have a simulation domain of 1.5 on x and we use sz = 16. Consider now
	 * to have particles with spacing 0.1 on x. if we define a sub box that on extend from 0.65
	 *  to 0.95 the first fluid particle
	 * is at 0.70 and the last is at 0.90
	 *
	 * \return an iterator to the selected particles
	 *
	 */
	template<unsigned int dim, typename T, typename aggr, typename layout, template <typename> class layout_base,  typename Decomposition> static PointIterator<dim,T,Decomposition>
	DrawBox(vector_dist<dim,T,aggr,layout,layout_base,Decomposition> & vd,
			 size_t (& sz)[dim],
			 Box<dim,T> & domain,
			 Box<dim,T> & sub)
	{
		return PointIterator<dim,T,Decomposition>(vd.getDecomposition(),sz,vd.getDecomposition().getDomain(),sub);
	}

};


#endif /* OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLES_HPP_ */
