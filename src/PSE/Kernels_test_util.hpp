/*
 * Kernels_unit_tests.hpp
 *
 *  Created on: Feb 16, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_PSE_KERNELS_TEST_UTIL_HPP_
#define OPENFPM_NUMERICS_SRC_PSE_KERNELS_TEST_UTIL_HPP_

#include "Kernels.hpp"
#include "Space/Ghost.hpp"
#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"

struct PSEError
{
	double l2_error;
	double linf_error;
};


template<typename T> T f_xex2(T x)
{
	return x*exp(-x*x);
}

template<typename T> T f_xex2(Point<1,T> & x)
{
	return x.get(0)*exp(-x.get(0)*x.get(0));
}

template<typename T> T Lapf_xex2(Point<1,T> & x)
{
	return 2.0*x.get(0)*(2.0*x.get(0)*x.get(0) - 3.0)*exp(-x.get(0)*x.get(0));
}

/*
 * Given the Number of particles, it calculate the error produced by the standard
 * second order PSE kernel to approximate the laplacian of x*e^{-x^2} in the point
 * 0.5 (1D)
 *
 * The spacing h is calculated as
 *
 * \$ h = 1.0 / N_{part} \$
 *
 * epsilon of the 1D kernel is calculated as
 *
 * \$ \epsilon=overlap * h \$
 *
 * this mean that
 *
 * \$ overlap = \frac{\epsilon}{h} \$
 *
 * While the particles that are considered in the neighborhood of 0.5 are
 * calculated as
 *
 * \$ \N_contrib = 20*eps/spacing \$
 *
 * \param N_part Number of particles in the domain
 * \param overlap overlap parameter
 *
 *
 */
template<typename T, typename Kernel> void PSE_test(size_t Npart, size_t overlap,struct PSEError & error)
{
	// The domain
	Box<1,T> box({0.0},{4.0});

	// Calculated spacing
	T spacing = box.getHigh(0) / Npart;

	// Epsilon of the particle kernel
	const T eps = overlap*spacing;

	// Laplacian PSE kernel 1 dimension, on double, second order
	Kernel lker(eps);

	// For convergence we need less particles
	Npart = static_cast<size_t>(20*eps/spacing);

	// Middle particle
	long int mp = Npart / 2;

    size_t bc[1]={NON_PERIODIC};
	Ghost<1,T> g(20*eps);

	vector_dist<1,T, aggregate<T>, CartDecomposition<1,T> > vd(Npart,box,bc,g);

	auto it2 = vd.getIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		// set the position of the particles
		vd.template getPos<0>(key)[0] = 0.448000 - ((long int)key.getKey() - mp) * spacing;
		//set the property of the particles
		vd.template getProp<0>(key) = f_xex2(vd.template getPos<0>(key)[0]);

		++it2;
	}

	vect_dist_key_dx key;
	key.setKey(mp);

    // get and construct the Cell list
    auto cl = vd.getCellList(0.5);

    // Maximum infinity norm
    double linf = 0.0;

    // PSE integration accumulator
    T pse = 0;

    // Get the position of the particle
    Point<1,T> p = vd.template getPos<0>(key);

    // Get f(x) at the position of the particle
    T prp_x = vd.template getProp<0>(key);

    // Get the neighborhood of the particles
    auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
    while(NN.isNext())
    {
    	auto nnp = NN.get();

    	// Calculate contribution given by the kernel value at position p,
    	// given by the Near particle, exclude itself
    	if (nnp != key.getKey())
    	{
    		// W(x-y)
    		T ker = lker.value(p,vd.template getPos<0>(nnp));

    		// f(y)
    		T prp_y = vd.template getProp<0>(nnp);

    		// 1.0/(eps)^2 [f(y)-f(x)]*W(x,y)*V_q
    		T prp = 1.0/eps/eps * spacing * (prp_y - prp_x);
    		pse += prp * ker;
    	}

    	// Next particle
    	++NN;
    }

	// Here we calculate the L_infinity norm or the maximum difference
	// of the analytic solution from the PSE calculated

	T sol = Lapf_xex2(p);

	if (fabs(pse - sol) > linf)
		linf = static_cast<double>(fabs(pse - sol));


    error.linf_error = (double)linf;
}


#endif /* OPENFPM_NUMERICS_SRC_PSE_KERNELS_TEST_UTIL_HPP_ */
