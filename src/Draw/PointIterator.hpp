/*
 * PointIterator.hpp
 *
 *  Created on: Jan 3, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_DRAW_POINTITERATOR_HPP_
#define OPENFPM_NUMERICS_SRC_DRAW_POINTITERATOR_HPP_


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
template<unsigned int dim, typename T, typename Decomposition>
class PointIterator: protected grid_dist_id_iterator_dec<Decomposition>
{
	//! Actual point
	Point<dim,T> ap;

	//! sub_domain
	Box<dim,T> sub_domain;

	//! domain
	Box<dim,T> domain;

	//! Spacing
	T sp[dim];

	static grid_key_dx<dim> getStart(size_t (& gs)[dim], Box<dim,T> & dom, Box<dim,T> & sub_dom, T (& sp)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
			sp[i] = (dom.getHigh(i) - dom.getLow(i)) / (gs[i] -1);

		grid_key_dx<dim> pkey;

		for (size_t i = 0 ; i < dim ; i++)
			pkey.set_d(i,std::ceil( (sub_dom.getLow(i) - dom.getLow(i)) / sp[i]));

		return pkey;
	}

	static grid_key_dx<dim> getStop(size_t (& gs)[dim], Box<dim,T> & dom, Box<dim,T> & sub_dom, T (& sp)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
			sp[i] = (dom.getHigh(i) - dom.getLow(i)) / (gs[i] - 1);

		grid_key_dx<dim> pkey;

		for (size_t i = 0 ; i < dim ; i++)
			pkey.set_d(i,std::floor( (sub_dom.getHigh(i) - dom.getLow(i)) / sp[i]));

		return pkey;
	}

	void calculateAp()
	{
		grid_key_dx<dim> key = grid_dist_id_iterator_dec<Decomposition>::get();

		for (size_t i = 0 ; i < dim ; i++)
			ap.get(i) = key.get(i) * sp[i] + domain.getLow(i);
	}

public:

	/*! \brief Draw Particles
	 *
	 * \param sp grid spacing
	 *
	 */
	PointIterator( Decomposition & dec, size_t (& sz)[dim], Box<dim,T> & domain, Box<dim,T> & sub_domain)
	:grid_dist_id_iterator_dec<Decomposition>(dec, sz, getStart(sz,domain,sub_domain,sp), getStop(sz,domain,sub_domain,sp)),sub_domain(sub_domain),domain(domain)
	{
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
	PointIterator & operator++()
	{
		grid_dist_id_iterator_dec<Decomposition>::operator++();
		calculateAp();

		return *this;
	}

	bool isNext()
	{
		return grid_dist_id_iterator_dec<Decomposition>::isNext();
	}

	/*! \brief Return the real Margin of the box
	 *
	 * For example consider to have a domain (0.0,0.0) to (1.0,1.0)
	 * 11 points in each direction (spacing 0.1). if we have a sub-domain
	 * (0.15,0.15) to (0.55,0.55) getBoxMargins return the box (0.2,0.2) to
	 * (0.5,0.5). the smallest box that enclose the points that the point
	 * iterator is going to give
	 *
	 */
	Box<dim,T> getBoxMargins()
	{
		Box<dim,T> box;

		grid_key_dx<dim> start = grid_dist_id_iterator_dec<Decomposition>::getStart();
		grid_key_dx<dim> stop = grid_dist_id_iterator_dec<Decomposition>::getStop();

		for (size_t i = 0 ; i < dim ; i++)
		{
			box.setLow(i, start.get(i)*sp[i] + domain.getLow(i));
			box.setHigh(i, stop.get(i)*sp[i] + domain.getLow(i));
		}

		return box;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_DRAW_POINTITERATOR_HPP_ */
