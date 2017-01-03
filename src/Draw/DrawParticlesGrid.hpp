/*
 * DrawParticlesGrid.hpp
 *
 *  Created on: Jan 3, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLESGRID_HPP_
#define OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLESGRID_HPP_

/*! \brief this class draw particles on subset of grid-like position
 *
 \verbatim
                  Point B
                    |                 (0.6,0.6)
  +   +   +   +   + |  +   +   +   +
    +---------------+
  + | +   +   +   + | +   +   +   +
    |               |
  + | +   +   +   + | +   +   +   +
    |               |
  + | +   +   +   + | +   +   +   +
    |               |
  + | +   +   +   + | +   +   +   +
    |               |
  + | +   +   +   + | +   +   +   +
    |               |
  + | +   +   +   + | +   +   +   +
    |               |
  + | +   +   +   + | +   +   +   +
    +---------------+
  +   +   +   +   + ^ +   +   +   +
(-1.0,-1.0)         |
                    |
                   (Point A)


 \endverbatim
 *
 *
 *  Suppose we have a grid 9x9 from (-1.0,-1.0 to 0.6,0.6)
 *
 *  Defined a Box (-0.9,-0.9 to -0.1,0.5)
 *
 *  This class will return the points
 *
 *  (-0.8 -0.8)
 *  (-0.6 -0.8)
 *  (-0.4 -0.8)
 *  ...
 *  (-0.1,-0.8) Point A
 *  ...
 *  (-0.1,0.5) Point B
 *
 */
template<unsigned int dim, typename T>
class DrawParticlesGrid
{
	/*! \brief Draw Particles
	 *
	 * \param sp grid spacing
	 *
	 */
	DrawParticlesGrid(T (& sp)[dim], Box<dim,T> & domain)
	{
	}

	PointIterator getPointIterator()
	{

	}
};


#endif /* OPENFPM_NUMERICS_SRC_DRAW_DRAWPARTICLESGRID_HPP_ */
