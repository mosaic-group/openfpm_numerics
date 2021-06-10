//
// Created by Abhinav Singh on 10.06.21.
// Modified by Justina Stark on 10.06.21
//

#include "config.h"
#include "util/common.hpp"
#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
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


BOOST_AUTO_TEST_SUITE(odeInt_BASE_tests)
	BOOST_AUTO_TEST_CASE(mirror_base_test) {
		const size_t sz[3] = {18, 18, 5};
		double boxsize = 10.0   ;
		double spacing_p = 10.0 / (sz[0]-1);
		double spacing_np = 10.0 / (sz[0]-1);
		Box<3, double> box({0.0, 0.0, 0}, {boxsize, boxsize, (sz[2])*spacing_p});
		size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, PERIODIC};
		Ghost<3, double> ghost(2.9*spacing_np+spacing_np/8.0);
		
		typedef vector_dist_ws<3, double, aggregate<double, VectorS<3,double>, bool>> vd_type;
		vd_type Particles(0,box,bc,ghost);
		Particles.setPropNames({"Concentration","Normal"});
		
		//Grid Like Particles
		auto it = Particles.getGridIterator(sz);
		while (it.isNext()) {
			auto key = it.get();
			double x = key.get(0) * spacing_np;
			double y = key.get(1) * spacing_np;
			double z = key.get(2) * spacing_p;
			if(sqrt((5.0-x)*(5.0-x)+(5.0-y)*(5.0-y))<5.0-7*spacing_np/6.0)
			{
				Particles.add();
				if (key.get(0)==sz[0]-1)
				{
					Particles.getLastPos()[0] = boxsize-1e-6;
				}
				else
				{
					Particles.getLastPos()[0] = x;
				}
				
				if (key.get(1)==sz[1]-1)
				{
					Particles.getLastPos()[1] = boxsize-1e-6;
				}
				else
				{
					Particles.getLastPos()[1] = y;
				}
				Particles.getLastPos()[2] = z;
				Particles.getLastProp<1>()[0] = 0;
				Particles.getLastProp<1>()[1] = 0;
				Particles.getLastProp<1>()[2] = 0;
				Particles.getLastProp<0>() = x+y+z;
				Particles.getLastSubset(0);
				
				Particles.getLastProp<IS_SOURCE>() = false;
			}
			++it;
		}
		
		//Adding a layer to represent geometry better. (Mirroring should only have this layer as it is slightly inside radius and the gridlike particles have been set to 0 normal)
		int n_b=int(sz[0])*5;
		double radius = 5.0 - 3*spacing_np/4;
		//double Golden_angle=M_PI * (3.0 - sqrt(5.0));
		double Golden_angle=2.0*M_PI/double(n_b);
		int count = 0;
		for (int j=0;j<int(sz[2]);j++)
		{
			for(int i=1;i<=n_b;i++)
			{
				double Golden_theta = Golden_angle * i;
				double x = 5.0+cos(Golden_theta) * radius;
				double y = 5.0+sin(Golden_theta) * radius;
				double z = j*spacing_p;
				Particles.add();
				Particles.getLastPos()[0] = x;
				Particles.getLastPos()[1] = y;
				Particles.getLastPos()[2] = z;
				Particles.getLastSubset(0);
				Particles.getLastProp<1>()[0] = (x-5.0)/sqrt((x-5.0)*(x-5.0)+(y-5.0)*(y-5.0));
				Particles.getLastProp<1>()[1] = (y-5.0)/sqrt((x-5.0)*(x-5.0)+(y-5.0)*(y-5.0));
				Particles.getLastProp<1>()[2] = 0.0;
				Particles.getLastProp<0>() = x+y+z;
				
				Particles.getLastProp<IS_SOURCE>() = true;
				
				/*Particles.getLastProp<8>()[0] = 1.0 ;
				Particles.getLastProp<8>()[1] = std::atan2(sqrt(x*x+y*y),z);
				Particles.getLastProp<8>()[2] = std::atan2(y,x);*/
				count++;
			}
		}
		std::cout << "Number of particles with surface normal = " << count << std::endl;
		Particles.map();
		
		//We write the particles to check if the initialization is correct.
		Particles.deleteGhost();
		Particles.write("Init");
		Particles.ghost_get<0>();
		
//		//Now we construct the subsets based on the subset number.
//		vector_dist_subset<3, double, aggregate<double,VectorS<3,double>, bool>> Particles_bulk(Particles,0);
//
//
//		//We create aliases for referring to property and and positions.
//		auto Pos = getV<PROP_POS>(Particles);
//		auto Concentration = getV<0>(Particles);
//		auto Normal = getV<1>(Particles);
//
//		//We create aliases for referring to the subset properties.
//		auto Concentration_bulk = getV<0>(Particles_bulk);
//		auto Normal_bulk = getV<1>(Particles_bulk);
		
		//Here Mirror Particle and do method of Images and check if it matches  property 0 mirroring (x+y+z of the mirror).
		
		openfpm::vector<vect_dist_key_dx> keys_source;
		
		auto dom = Particles.getDomainIterator();
		while(dom.isNext())
		{
			auto key = dom.get();
			if (Particles.template getProp<IS_SOURCE>(key))
			{
				keys_source.add(key);
			}
			++dom;
		}
		std::cout << "keys_source.size() = " << keys_source.size() << std::endl;
		std::cout << "Particles.size_local() = " << Particles.size_local() << std::endl;
		int ps_before = Particles.size_local_with_ghost();
		std::cout << "Particles.size_local_with_ghost() before mirroring = " << ps_before << std::endl;
		MethodOfImages<NORMAL, vd_type> NBCs(Particles, keys_source, 0, 1);
		NBCs.get_mirror_particles(Particles);
		std::cout << "Particles.size_local_with_ghost() after mirroring = " << Particles.size_local_with_ghost() <<
				std::endl;
		
		std::cout << "Total number of particles expected after mirroring = " << ps_before + keys_source.size() <<
		std::endl;
		
		Particles.write("Cylinder_with_mirror_particles");
		
		
		
	}
BOOST_AUTO_TEST_SUITE_END()

