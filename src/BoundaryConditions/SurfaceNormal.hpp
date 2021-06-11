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
	}
}

#if 0
template <typename key_g_type, int dim>
bool is_corner(key_g_type key_g)
{
	bool nodetype = false;
	switch(dim)
	{
		case 1:
			nodetype = (key_g.get(0) < 1 || key_g.get(0) > grid.size(0) - 1);
			break;
		case 2:
			nodetype = ()
			
	}
	
	return nodetype;
}

template <typename key_g_type, int dim>
bool is_edge(key_g_type key_g)
{

}

template <typename key_g_type, int dim>
bool is_bulk(key_g_type key_g, int border_width_in_nodes)
{
	for (int d=0; d<dim; d++)
	{
		if (key_g.get(d) < border_width_in_nodes || key_g.get(d) > grid.size(d) - border_width_in_nodes)
		{return false;}
	}
	return true;
}

template <typename key_g_type, int dim>
bool is_boundary(key_g_type key_g, int border_width_in_nodes)
{
	for (int d=0; d<dim; d++)
	{
		if (key_g.get(d) > border_width_in_nodes || key_g.get(d) < grid.size(d) - border_width_in_nodes)
		{return false;}
	}
	return true;
}



template <size_t Normal, size_t SDF, typename grid_type>
void get_surface_normal_box(grid_type & grid, bool unit_vector=false)
{
	auto dom = grid.getDomainIterator();
	while (dom.isNext())
	{
		auto key   = dom.get();
		auto key_g = grid.getGKey(key);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		if (key_g.get(d) > 0 && key_g.get(d) < grid.size(d) - 1) // if point lays not at the border of the grid
		{
			grid.template get<Phi_grad> (key) [d] = upwind_FD<Phi, Phi_0_sign>(grid, key, d);
		}
		
		else if (key_g.get(d) == 0) // if point lays at left border, use right sided kernel
		{
			grid.template get<Phi_grad> (key) [d] = forward_FD<Phi>(grid, key, d);
		}
		
		else if (key_g.get(d) >= grid.size(d) - 1) // if point lays at right border, use left sided kernel
		{
			grid.template get<Phi_grad> (key) [d] = backward_FD<Phi>(grid, key, d);
		}
		
		
		
		
		++dom;
	}
}


#endif

#endif //BOUNDARY_CONDITION_SURFACENORMAL_HPP
