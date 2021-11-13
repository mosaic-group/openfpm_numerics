//
// Created by jstark on 07.06.21.
//

#ifndef OPENFPM_NUMERICS_HELPFUNCTIONSFORTESTINGNBCS_HPP
#define OPENFPM_NUMERICS_HELPFUNCTIONSFORTESTINGNBCS_HPP

#include <iostream>
#include "level_set/redistancing_Sussman/NarrowBand.hpp"

template <size_t Phi_SDF_grid, size_t Phi_SDF, size_t dPhi, size_t dPhi_magn, typename grid_in_type, typename
vd_in_type>
void get_diffusion_domain(grid_in_type & g_dist, vd_in_type & vd, const double b_low, const double b_up)
{
	const double EPSILON = std::numeric_limits<double>::epsilon();
	const double _b_low = b_low + EPSILON;
	const double _b_up  = b_up - EPSILON;
	NarrowBand<grid_in_type> narrowBand(g_dist, _b_low, _b_up); // Instantiation of NarrowBand class
	narrowBand.template get_narrow_band<Phi_SDF_grid, Phi_SDF, dPhi, dPhi_magn>(g_dist, vd);
}

/**@brief Fill pid_source with IDs of particles whose SDF is >= b_low and <= b_up.
 *
 * @param vd Input particle vector_dist of type vd_type.
 * @param b_low Lower bound for the SDF value a particle must have in order to be a source particle.
 * @param b_up Upper bound for the SDF value a particle must have in order to be a source particle.
 */
template <size_t Phi_SDF, typename vd_type>
void get_source_particle_ids(vd_type & vd, openfpm::vector<vect_dist_key_dx> & keys_source, const double b_low, const
double
b_up)
{
	const double EPSILON = std::numeric_limits<double>::epsilon();
	auto dom = vd.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		if (vd.template getProp<Phi_SDF>(key) >= b_low + EPSILON && vd.template getProp<Phi_SDF>(key) <= b_up + EPSILON)
		{
			keys_source.add(key);
		}
		++dom;
	}
}


template <size_t R, size_t U, typename vd_type>
void get_IC_Eq28(vd_type &vd, const double a)
{
	auto dom = vd.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		double r = vd.template getProp<R>(key);
		if (abs(r) < a) vd.template getProp<U>(key) = cos(M_PI * r / (2*a)) * cos(M_PI * r / (2*a));
		else vd.template getProp<U>(key) = 0;
		++dom;
	}
}

template <size_t R, size_t U, typename vd_type>
void get_IC_Eq28(vd_type &vd, const double a, const openfpm::vector<aggregate<int>> & pids)
{
	for (int i = 0; i < pids.size(); i++)
	{
		auto key = pids.get<0>(i);
		double r = vd.template getProp<R>(key);
		if (abs(r) < a) vd.template getProp<U>(key) = cos(M_PI * r / (2*a)) * cos(M_PI * r / (2*a));
		else vd.template getProp<U>(key) = 0;
	}
}

template <size_t U, typename vd_type>
void get_FS_Eq29(vd_type &vd, const double a, const double R_disk)
{
	const double R = R_disk;
	double A = 0;
	if (a < R) A = M_PI * a*a/2 - 2*R*R/M_PI;
	if (a >= R) A = M_PI * R*R/2 + a*a*( cos(M_PI * R / a) -1 ) / M_PI + a * R * sin(M_PI * R / a);
	
	auto u = getV<U>(vd);
	u = A / (M_PI * R*R);
	
	std::cout << "Analytical steady state concentration = " << A / (M_PI * R*R) << std::endl;
}


#endif //OPENFPM_NUMERICS_HELPFUNCTIONSFORTESTINGNBCS_HPP
