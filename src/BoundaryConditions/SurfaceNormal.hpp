//
// Created by jstark on 28.05.21.
//

#ifndef BOUNDARY_CONDITION_SURFACENORMAL_HPP
#define BOUNDARY_CONDITION_SURFACENORMAL_HPP

#include "cmath"

template <size_t Phi_SDF, size_t Phi_Gradient, size_t Normal, typename vd_type>
void get_surface_normal_sdf(vd_type & vd, bool unit_vector=false)
{
	auto dom = vd.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		
		// Get the magnitude of the gradient
		double sum = 0;
		for(size_t d = 0; d < vd_type::dims; d++)
		{
			sum += vd.template getProp<Phi_Gradient>(key)[d] * vd.template getProp<Phi_Gradient>(key)[d];
		}
		
		double gradient_magnitude = sqrt(sum);
		for(size_t d = 0; d < vd_type::dims; d++)
		{
			if(unit_vector)
			{
				vd.template getProp<Normal>(key)[d] = - vd.template getProp<Phi_Gradient>(key)[d] / gradient_magnitude;
			}
			else
			{
				vd.template getProp<Normal>(key)[d] =
						- vd.template getProp<Phi_Gradient>(key)[d] / gradient_magnitude
								* vd.template getProp<Phi_SDF>(key) * 2.0;
			}
		}
		++dom;
	}
}

template <size_t Phi_SDF, size_t Phi_Gradient, size_t Normal, typename vd_type>
void get_surface_normal_sdf_subset(vd_type & vd, const openfpm::vector<vect_dist_key_dx> & keys_subset, bool
unit_vector=false)
{
	for (int i = 0; i < keys_subset.size(); i++)
	{
		auto key = keys_subset.get(i);
		
		// Get the magnitude of the gradient
		double sum = 0;
		for(size_t d = 0; d < vd_type::dims; d++)
		{
			sum += vd.template getProp<Phi_Gradient>(key)[d] * vd.template getProp<Phi_Gradient>(key)[d];
		}
		
		double gradient_magnitude = sqrt(sum);
		for(size_t d = 0; d < vd_type::dims; d++)
		{
			vd.template getProp<Normal>(key)[d] = - vd.template getProp<Phi_Gradient>(key)[d] / gradient_magnitude;
		}
		if(!unit_vector)
		{
			vd.template getProp<Normal>(key) = vd.template getProp<Normal>(key) * 2.0 * vd.template getProp<Phi_SDF>
					(key);
		}
		
		std::cout << "n.norm = " << vd.template getProp<Normal>(key).norm() << ", 2xsdf = " <<
				vd.template getProp<Phi_SDF>(key) * 2.0 << std::endl;
	}
}



#endif //BOUNDARY_CONDITION_SURFACENORMAL_HPP
