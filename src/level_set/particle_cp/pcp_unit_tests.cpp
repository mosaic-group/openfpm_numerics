/* Unit tests for the particle closest point method
 * Date : 29 June 2023
 * Author : lschulze
 */

#include<iostream>
#include <boost/test/unit_test_log.hpp>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <math.h>
//#include <sys/_types/_size_t.h>
#include "Vector/vector_dist.hpp"
#include "DCPSE/Dcpse.hpp"
#include "DCPSE/MonomialBasis.hpp"
#include "Draw/DrawParticles.hpp"
#include "particle_cp.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include <chrono>
#include "pcp_unit_test_helpfunctions.h"

typedef struct EllipseParameters{
    double origin[3];
    double radiusA;
    double radiusB;
    double radiusC;
    double eccentricity;
} EllipseParams;

// Generate an ellipsoid initial levelset signed distance function
template<typename particles_type, typename iterator_type, size_t sdf, size_t surf_flag, size_t ref_cp>
void initializeLSEllipsoid(particles_type &vd, iterator_type particle_it, const EllipseParams &params, double bandwidth, double perturb_factor, double H)
{
    while(particle_it.isNext())
    {
        double posx = particle_it.get().get(0);
        double posy = particle_it.get().get(1);
        double posz = particle_it.get().get(2);

	double ref_cpx = 0.0;
	double ref_cpy = 0.0;
	double ref_cpz = 0.0;

	double dist = DistancePointEllipsoid(params.radiusA, params.radiusB, params.radiusC, abs(posx), abs(posy), abs(posz), ref_cpx, ref_cpy, ref_cpz);
        if (abs(dist) < bandwidth/2.0)
	{
		posx = posx + perturb_factor*randMinusOneToOne()*H;
		posy = posy + perturb_factor*randMinusOneToOne()*H;
		posz = posz + perturb_factor*randMinusOneToOne()*H;

		dist = DistancePointEllipsoid(params.radiusA, params.radiusB, params.radiusC, abs(posx), abs(posy), abs(posz), ref_cpx, ref_cpy, ref_cpz);
		vd.add();
		vd.getLastPos()[0] = posx;
		vd.getLastPos()[1] = posy;
		vd.getLastPos()[2] = posz;

		double phi_val = sqrt(((posx - params.origin[0])/params.radiusA)*((posx - params.origin[0])/params.radiusA) + ((posy - params.origin[1])/params.radiusB)*((posy - params.origin[1])/params.radiusB) + ((posz - params.origin[2])/params.radiusC)*((posz - params.origin[2])/params.radiusC)) - 1.0;
        	vd.template getLastProp<sdf>() = phi_val;
		vd.template getLastProp<ref_cp>()[0] = return_sign(posx)*ref_cpx;
		vd.template getLastProp<ref_cp>()[1] = return_sign(posy)*ref_cpy;
		vd.template getLastProp<ref_cp>()[2] = return_sign(posz)*ref_cpz;
		vd.template getLastProp<surf_flag>() = 0;
	}

	++particle_it;
    }
}

BOOST_AUTO_TEST_SUITE( particle_closest_point_test )

BOOST_AUTO_TEST_CASE( ellipsoid )
{
	// simulation params
	constexpr int poly_order = 4;
	const double H = 1.0/64.0;
	const double perturb_factor = 0.3;

	// pcp params
	const double bandwidth = 12.0*H;
	const double regression_rcut_factor = 2.4;
	const double threshold = 1e-13;

	// initialize domain and particles
	const double l = 2.0;
	Box<3, double> domain({-l/2.0, -l/3.0, -l/3.0}, {l/2.0, l/3.0, l/3.0});
	size_t sz[3] = {(size_t)(l/H + 0.5), (size_t)((2.0/3.0)*l/H + 0.5), (size_t)((2.0/3.0)*l/H + 0.5)};
	size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
	Ghost<3, double> g(bandwidth);

	constexpr int sdf = 0;
	constexpr int cp = 1;
	constexpr int normal = 2;
	constexpr int curvature = 3;
	constexpr int surf_flag = 4;
	constexpr int ref_cp = 5;
	typedef vector_dist<3, double, aggregate<double, Point<3, double>, Point<3, double>, double, int, Point<3, double>>> particles;
	//					 |		|		|		|      |		|
	//					SDF	closest point		|  	   curvature surface flag 	|
	//								surface normal	   		reference closest point

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	openfpm::vector<std::string> names({"sdf","cp","normal","curvature","ref_cp"});
	vd.setPropNames(names);
	EllipseParameters params;
	for (int k = 0; k < 3; k++) params.origin[k] = 0.0;
	params.radiusA = 0.75;
	params.radiusB = 0.5;
	params.radiusC = 0.5;

	auto particle_it = DrawParticles::DrawBox(vd, sz, domain, domain);

	// initialize spurious sdf values and reference solution
	initializeLSEllipsoid<particles, decltype(particle_it), sdf, surf_flag, ref_cp>(vd, particle_it, params, bandwidth, perturb_factor, H);

	//vd.write("pcpunittest_init");
	// initialize and run pcp redistancing
	Redist_options rdistoptions;
	rdistoptions.minter_poly_degree = poly_order;
	rdistoptions.H = H;
	rdistoptions.r_cutoff_factor = regression_rcut_factor;
	rdistoptions.sampling_radius = 0.75*bandwidth;
	rdistoptions.tolerance = threshold;
	rdistoptions.write_cp = 1;
	rdistoptions.compute_normals = 1;
	rdistoptions.compute_curvatures = 1;
	rdistoptions.only_narrowband = 0;

	static constexpr unsigned int num_coeffs = minter_lp_degree_one_num_coeffs(3, poly_order);

	particle_cp_redistancing<particles, sdf, cp, normal, curvature, num_coeffs> pcprdist(vd, rdistoptions);
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	pcprdist.run_redistancing();
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference for pcp redistancing = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

	//vd.write("pcpunittest");
	// iterate through particles and compute error
	auto part = vd.getDomainIterator();
	double sdferr = 0.0;
	double sdfmaxerr = 0.0;
	double sdfmeanerr = 0.0;
	double cperr = 0.0;
	double cpmaxerr = 0.0;
	double normalerr = 0.0;
	double normalmaxerr = 0.0;
	double curvatureerr = 0.0;
	double curvaturemaxerr = 0.0;
	int particlecount = 0;

	while(part.isNext())
	{
		auto a = part.get();
		Point<3, double> pos = vd.getPos(a);
		double reference_sdf = return_sign(sqrt(((pos[0] - params.origin[0])/params.radiusA)*((pos[0] - params.origin[0])/params.radiusA) + ((pos[1] - params.origin[1])/params.radiusB)*((pos[1] - params.origin[1])/params.radiusB) + ((pos[2] - params.origin[2])/params.radiusC)*((pos[2] - params.origin[2])/params.radiusC)) - 1.0)*norm(vd.getProp<ref_cp>(a) - pos);
		Point<3, double> reference_normal;
		double reference_curvature;

		// compute reference SDF value w/ reference closest point and check computed SDF value
		sdferr = abs(vd.getProp<sdf>(a) - reference_sdf);
		if (sdferr > sdfmaxerr) sdfmaxerr = sdferr;
		sdfmeanerr += sdferr;
		++particlecount;

		// check computed closest point against reference closest point
		cperr = norm(vd.getProp<ref_cp>(a) - vd.getProp<cp>(a));
		if (cperr > cpmaxerr) cpmaxerr = cperr;

		// compute reference normal and check computed normal
		reference_normal[0] = 2*vd.getProp<ref_cp>(a)[0]/(params.radiusA*params.radiusA);
		reference_normal[1] = 2*vd.getProp<ref_cp>(a)[1]/(params.radiusB*params.radiusB);
		reference_normal[2] = 2*vd.getProp<ref_cp>(a)[2]/(params.radiusC*params.radiusC);
		reference_normal = return_sign(reference_sdf)*reference_normal/norm(reference_normal);
		normalerr = norm(reference_normal - vd.getProp<normal>(a));
		if (normalerr > normalmaxerr) normalmaxerr = normalerr;

		// compute reference curvature and check computed curvature
		reference_curvature = (std::abs(vd.getProp<ref_cp>(a)[0]*vd.getProp<ref_cp>(a)[0] + vd.getProp<ref_cp>(a)[1]*vd.getProp<ref_cp>(a)[1] + vd.getProp<ref_cp>(a)[2]*vd.getProp<ref_cp>(a)[2] - params.radiusA*params.radiusA  - params.radiusB*params.radiusB - params.radiusC*params.radiusC))/(2*params.radiusA*params.radiusA*params.radiusB*params.radiusB*params.radiusC*params.radiusC*std::pow(vd.getProp<ref_cp>(a)[0]*vd.getProp<ref_cp>(a)[0]/std::pow(params.radiusA, 4) + vd.getProp<ref_cp>(a)[1]*vd.getProp<ref_cp>(a)[1]/std::pow(params.radiusB, 4) + vd.getProp<ref_cp>(a)[2]*vd.getProp<ref_cp>(a)[2]/std::pow(params.radiusC, 4), 1.5));
		curvatureerr = abs(reference_curvature - vd.getProp<curvature>(a));
		if (curvatureerr > curvaturemaxerr) curvaturemaxerr = curvatureerr;

		++part;
	}
	sdfmeanerr = sdfmeanerr/particlecount;
	std::cout<<"Mean error for sdf is: "<<sdfmeanerr<<std::endl;
	std::cout<<"Maximum error for sdf is: "<<sdfmaxerr<<std::endl;
        std::cout<<"Maximum error for the closest point is: "<<cpmaxerr<<std::endl;
        std::cout<<"Maximum error for surface normal is: "<<normalmaxerr<<std::endl;
        std::cout<<"Maximum error for curvature is: "<<curvaturemaxerr<<std::endl;

    	double tolerance = 1e-7;
   	bool check;
    	if (std::abs(sdfmaxerr) < tolerance)
        	check = true;
    	else
        	check = false;

    	BOOST_TEST( check );

}

BOOST_AUTO_TEST_SUITE_END()
