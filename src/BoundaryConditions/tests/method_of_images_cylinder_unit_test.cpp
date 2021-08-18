//
// Created by Abhinav Singh & Justina Stark on 10.06.21.
//
#define BOOST_TEST_DYN_LINK
#include <iostream>
#include "config.h"
#include "util/common.hpp"
#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>

#include "Vector/vector_dist_subset.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

// For the Neumann BCs (Method of images)
#include "BoundaryConditions/MethodOfImages.hpp"


constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;
constexpr size_t CONCENTRATION = 0;
constexpr size_t NORMAL        = 1;
constexpr size_t IS_SOURCE     = 2;
constexpr size_t subset_id_real   = 0;
constexpr size_t subset_id_mirror = 1;



BOOST_AUTO_TEST_SUITE(MethodOfImagesTestSuite)
	BOOST_AUTO_TEST_CASE(mirror_cylinder_base_test) {
		const size_t sz[3] = {18, 18, 5};
		double boxsize = 20.0;
		double spacing_p = 10.0 / (sz[x]-1);
		double spacing_np = 10.0 / (sz[x]-1);
		Box<3, double> box({-2.0, -2.0, 0}, {boxsize, boxsize, (sz[z])*spacing_p});
		size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, PERIODIC};
//		Ghost<3, double> ghost(2.9*spacing_np+spacing_np/8.0);
		const double ghost_layer_thickness = 2.9*spacing_np+spacing_np/8.0;
//		const double ghost_layer_thickness = 5.0 - 3*spacing_np/4;
		Ghost<3, double> ghost(ghost_layer_thickness);
		
		typedef vector_dist_ws<3, double, aggregate<double, VectorS<3,double>, int>> vd_type;
		vd_type Particles(0,box,bc,ghost);
		Particles.setPropNames({"Concentration", "Normal", "is_source"});
		
		//Grid Like Particles
		auto it = Particles.getGridIterator(sz);
		while (it.isNext())
		{
			auto key = it.get();
			double x_coord = key.get(x) * spacing_np;
			double y_coord = key.get(y) * spacing_np;
			double z_coord = key.get(z) * spacing_p;
			
			if(sqrt((5.0-x_coord)*(5.0-x_coord)+(5.0-y_coord)*(5.0-y_coord))<5.0-7*spacing_np/6.0)
			{
				Particles.add();
				if (key.get(x)==sz[x]-1)
				{
					Particles.getLastPos()[x] = boxsize-1e-6;
				}
				else
				{
					Particles.getLastPos()[x] = x_coord;
				}
				
				if (key.get(y)==sz[y]-1)
				{
					Particles.getLastPos()[y] = boxsize-1e-6;
				}
				else
				{
					Particles.getLastPos()[y] = y_coord;
				}
				Particles.getLastPos()[z] = z_coord;
				Particles.getLastProp<NORMAL>()[x] = 0;
				Particles.getLastProp<NORMAL>()[y] = 0;
				Particles.getLastProp<NORMAL>()[z] = 0;
				Particles.getLastProp<CONCENTRATION>() = x_coord+y_coord+z_coord;
				
				Particles.getLastProp<IS_SOURCE>() = 0;
			}
			++it;
		}
		
		//Adding a layer to represent geometry better. (Mirroring should only have this layer as it is slightly inside radius and the gridlike particles have been set to 0 normal)
		int n_b=int(sz[0])*5;
		double radius = 5.0 - 3*spacing_np/4;
		//double Golden_angle=M_PI * (3.0 - sqrt(5.0));
		double Golden_angle=2.0*M_PI/double(n_b);
		int number_of_border_particles = 0;
		for (int j=0;j<int(sz[z]);j++)
		{
			for(int i=1;i<=n_b;i++)
			{
				double Golden_theta = Golden_angle * i;
				double x_coord = 5.0+cos(Golden_theta) * radius;
				double y_coord = 5.0+sin(Golden_theta) * radius;
				double z_coord = j*spacing_p;
				Particles.add();
				Particles.getLastPos()[x] = x_coord;
				Particles.getLastPos()[y] = y_coord;
				Particles.getLastPos()[z] = z_coord;

				Particles.getLastProp<NORMAL>()[x] = (x_coord-5.0)/sqrt((x_coord-5.0)*(x_coord-5.0)+(y_coord-5.0)*(y_coord-5.0));
				Particles.getLastProp<NORMAL>()[y] = (y_coord-5.0)/sqrt((x_coord-5.0)*(x_coord-5.0)+(y_coord-5.0)*(y_coord-5.0));
				Particles.getLastProp<NORMAL>()[z] = 0.0;
				Particles.getLastProp<CONCENTRATION>() = x_coord+y_coord+z_coord;
				
				Particles.getLastProp<IS_SOURCE>() = 1;

				number_of_border_particles++;
			}
		}
		auto &v_cl = create_vcluster();
		v_cl.sum(number_of_border_particles);
		v_cl.execute();
		
		if (v_cl.rank() == 0) std::cout << "Number of particles with surface normal = " << number_of_border_particles
		<< std::endl;
		Particles.map();
		Particles.ghost_get<NORMAL,IS_SOURCE>();
		//We write the particles to check if the initialization is correct.
		Particles.write("Init");
		
		//Here Mirror Particle and do method of Images and check if it matches  property 0 mirroring (x+y+z of the mirror).
		
		openfpm::vector<vect_dist_key_dx> keys_source;
		
		auto dom = Particles.getDomainIterator();
		while(dom.isNext())
		{
			auto key = dom.get();
			if (Particles.template getProp<IS_SOURCE>(key) == 1)
			{
				keys_source.add(key);
			}
			++dom;
		}
		size_t number_of_source_particles = keys_source.size();
		v_cl.sum(number_of_source_particles);
		v_cl.execute();
		size_t number_of_real_particle_no_ghost = Particles.size_local();
		size_t number_of_real_particle_with_ghost = Particles.size_local_with_ghost();
		
		/*
		if (v_cl.rank() == 0)
		{
			std::cout << "number_of_source_particles = " << number_of_source_particles << std::endl;
			std::cout << "number_of_real_particle_no_ghost = " << number_of_real_particle_no_ghost << std::endl;
			std::cout << "number_of_real_particle_with_ghost before mirroring = " << number_of_real_particle_with_ghost << std::endl;
		}
		*/
		
		
		// Apply Method of images to impose noflux Neumann Boundary Conditions
		MethodOfImages<NORMAL, vd_type> NBCs(Particles, keys_source, subset_id_real, subset_id_mirror);
		NBCs.get_mirror_particles(Particles);
		NBCs.apply_noflux<CONCENTRATION>(Particles);
		
		size_t number_of_mirror_particles = Particles.size_local() - number_of_real_particle_no_ghost;
		v_cl.sum(number_of_mirror_particles);
		v_cl.execute();
		
		if (v_cl.rank() == 0) std::cout << "Number of mirror particles = " << number_of_mirror_particles << std::endl;
		
		/*
		if (v_cl.rank() == 0)
		{
			std::cout << "number_of_real_particle_with_ghost + mirror particles = " << Particles.size_local_with_ghost() <<
					std::endl;
			std::cout << "Total number of particles expected after mirroring = " << number_of_real_particle_with_ghost +
					keys_source.size() << std::endl;
		}
		*/
		
		Particles.write("Cylinder_with_mirror_particles");
		
		BOOST_CHECK(number_of_source_particles == number_of_border_particles);
		BOOST_CHECK(number_of_mirror_particles == number_of_source_particles);
		
		for (int i = 0; i < keys_source.size(); ++i)
		{
			auto key_source = keys_source.get<0>(i); // Get key of one source particle
			auto key_mirror = NBCs.pid_mirror.get<0>(i); // Get key of corresponding mirror particle to that source
			// particle
			BOOST_CHECK(Particles.template getProp<CONCENTRATION>(key_mirror) == Particles.template getProp<CONCENTRATION>(key_source));
		}
		
	}
BOOST_AUTO_TEST_SUITE_END()

