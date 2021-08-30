//
// Created by jstark on 01.06.21.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

// Include header files for testing
#include "BoundaryConditions/MethodOfImages.hpp"
#include "BoundaryConditions/SurfaceNormal.hpp"
#include "libForTesting/HelpFunctionsForTestingNBCs.hpp"

BOOST_AUTO_TEST_SUITE(MethodOfImagesTestSuite)
	
	BOOST_AUTO_TEST_CASE(DiffusionWithNoFluxNBCs)
	{
		const size_t dims     = 2;
		const double L_low    = 0.0;
		const double L_up     = 1.0;
		const double R_disk   = 0.25;
		const double D = 0.005; // Diffusion constant.
		const double a = R_disk;
		const double t0 = 0; // Must be 0 in order to fit the IC from Eq. 28
		const double tmax = 20;
		constexpr size_t Phi_SDF = 0, dPhi = 1, dPhi_magn = 2, SurfaceNormal = 3, Conc_n = 4, Conc_nplus1 = 5, Lap_conc = 6,
				Radius = 7, Conc_exact = 8, Error = 9;
		constexpr size_t x = 0, y = 1;
		
//		for (size_t N = 32; N <= 512; N *= 2)
		size_t N = 32;
		{
			auto &v_cl = create_vcluster();
			
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Some indices for the grid / grid properties
			const size_t Phi_0_grid             = 0;
			const size_t Phi_SDF_grid           = 1;
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Here we create a 2D grid that stores 2 properties:
			// Prop1: store the initial Phi;
			// Prop2: here the re-initialized Phi (signed distance function) will be written to in the re-distancing step
			size_t sz[dims] = {N, N}; // Grid size in terms of number of grid points per dimension
			Box<dims, double> box({L_low, L_low}, {L_up, L_up}); // 2D
			Ghost<dims, long int> ghost(0);
			typedef aggregate<double, double> props;
			typedef grid_dist_id<dims, double, props > grid_in_type;
			grid_in_type g_dist(sz, box, ghost);
			g_dist.setPropNames({"Phi_0", "Phi_SDF"});
			g_dist.load("/Users/jstark/Desktop/diffusion/simple_geometries/get_diffusion_space/disk_analytical_sdf"
			            "/output_circle_" + std::to_string(N) + "/grid_disk_analyticalSdf.bin"); // Load SDF of a disk
			// from hdf5 file.
			double p_spacing = g_dist.spacing(x);
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Place particles into diffusion space according to SDF stored in grid. Copy Phi_SDF, dPhi and dPhi_magn from the
			// grid onto the particles.
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			size_t bc[dims] = {NON_PERIODIC, NON_PERIODIC};
			Ghost<dims, double> ghost_vd(p_spacing * 4);
			std::cout << "Ghost layer thickness: " << p_spacing * 10 << std::endl;
			typedef aggregate<double, double[dims], double, Point<dims, double>, double, double, double, double, double, double>
					props_diffSpace;
			typedef vector_dist_ws<dims, double, props_diffSpace> vd_type;
			vd_type vd(0, box, bc, ghost_vd);
			vd.setPropNames({"Phi_SDF", "dPhi", "dPhi_magn", "SurfaceNormal", "Conc_n", "Conc_n+1", "Lap_conc", "Radius",
			                 "Conc_exact", "Error"});
			
			
			// Get diffusion domain
			double dlow = 0, dup = 2*L_up;
			get_diffusion_domain<Phi_SDF_grid, Phi_SDF, dPhi, dPhi_magn>(g_dist, vd, dlow, dup);
//	vd.write(path_output + "/vd_diffusionSpace_real", FORMAT_BINARY); // VTK file
//	vd.save(path_output + "/vd_diffusionSpace_real.bin"); // HDF5 file
			
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//	Initialize the mass
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Get polar radius from cartesian coordinates
			double center_disk [] = {0.5, 0.5};
			get_radius<Radius>(vd, center_disk);
			
			// Initialize with exact solution for diffusion on disk at t=0
			// get_Eq27<Radius, Conc_n>(vd, a, D, t0);
			get_IC_Eq28<Radius, Conc_n>(vd, a);
			
			vd.template ghost_get<Conc_n>(KEEP_PROPERTIES);
			
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Place mirror particles for noflux boundary conditions.
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			size_t mirror_width = 3; // Mirror width in number of particles
//		size_t mirror_width = std::round(0.1/p_spacing); // Mirror width in number of particles
			
			if (v_cl.rank() == 0) std::cout << "Width of mirror layer: " << mirror_width << " particles." << std::endl;
			double b_low = 0, b_up = mirror_width * p_spacing;
			openfpm::vector<vect_dist_key_dx> keys_source;
			get_source_particle_ids<Phi_SDF>(vd, keys_source, b_low, b_up);
			
			get_surface_normal_sdf_subset<Phi_SDF, dPhi, SurfaceNormal>(vd, keys_source);
			
			MethodOfImages<SurfaceNormal, vd_type> NBCs(vd, keys_source, 0, 1);
			NBCs.get_mirror_particles(vd);
			
		}
		
		
	}
BOOST_AUTO_TEST_SUITE_END()
