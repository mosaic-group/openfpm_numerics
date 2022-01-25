//
// Created by lschulze on 2021-10-12
//
/**
 * @file particle_cp.hpp
 * @class particle_cp_redistancing
 *
 * @brief Class for reinitializing a level-set function into a signed distance function using the Closest-Point
 * method on particles.
 *
 * @details The redistancing scheme here is based on Saye, Robert. "High-order methods for computing distances
 * to implicitly defined surfaces." Communications in Applied Mathematics and Computational Science 9.1 (2014):
 * 107-141.
 * Within this file, it has been extended so that it also works on arbitrarily distributed particles. It mainly
 * consists of two steps: Firstly, an interpolation of the old signed-distance function values which yields
 * a polynomial whose zero level-set reconstructs the surface. Secondly, a subsequent optimisation in order to
 * determine the closest point on the surface to a given query point, allowing for an update of the signed
 * distance function
 *
 *
 * @author Lennart J. Schulze
 * @date October 2021
 */

#include <math.h>
#include <sys/_types/_size_t.h>
#include "Vector/vector_dist.hpp"
#include "DCPSE/Dcpse.hpp"
#include "DCPSE/MonomialBasis.hpp"
#include "Draw/DrawParticles.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include <chrono>

typedef vector_dist<2, double, aggregate<vect_dist_key_dx, int, int, int, double, double[2], double[16]>> particles_surface;
//											|			   |	  |		|	   |  		|			  |
//						key to particle in vd	surface flag num neibs	|	  sdf 	sample version	interpolation coefficients
//												flag for distinguishing closest particles to surface, for which the interpolation has to be done
struct Redist_options
{
	size_t max_iter = 1000;				// params for the Newton algorithm for the solution of the constrained optimization problem
	double tolerance = 1e-11;
	
	double H;							// interparticle spacing
	double r_cutoff_factor = 2.6;		// radius of neighbors to consider during interpolation, as factor of H
	double sampling_radius;
	std::string polynomial_degree = "bicubic";
	double barrier_coefficient = 0.0;	// coefficient for barrier term in optimisation
	double armijo_tau = 0.0;			// armijo line search parameters
	double armijo_c = 0.0;
	double support_prevent = 0.0;		// prevention parameter if support may be left
	int fail_projection = 0;			// prevention projection if support may be left
	int init_project = 0;
};

template <typename particles_in_type>
class particle_cp_redistancing
{
public:
	particle_cp_redistancing(particles_in_type & vd, Redist_options &redistOptions) : redistOptions(redistOptions),
															  	  	  	  	  	  	  vd_in(vd),
																					  vd_s(vd.getDecomposition(), 0),
																					  r_cutoff2(redistOptions.r_cutoff_factor*redistOptions.r_cutoff_factor*redistOptions.H*redistOptions.H)
	{
		// do constructor things
	}
	
	void run_redistancing()
	{
		std::cout<<"detecting..."<<std::endl;
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		detect_surface_particles();
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Time difference for detection= " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

		std::cout<<"interpolating..."<<std::endl;
		begin = std::chrono::steady_clock::now();
		interpolate_sdf_field();
		end = std::chrono::steady_clock::now();
		std::cout << "Time difference for interpolation= " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

		std::cout<<"optimizing..."<<std::endl;
		begin = std::chrono::steady_clock::now();
		find_closest_point();
		end = std::chrono::steady_clock::now();
		std::cout << "Time difference for optimization= " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
		std::cout<<"finishing..."<<std::endl;
	}
	
private:
	static constexpr size_t vdkey = 0;		//todo: this should be dispensable as well by adding the interesting properties to vd_s
	static constexpr size_t surf_flag = 1;	//todo: get rid of this flag as well, as no non-surface_flag particles should be in vd_s anyways
	static constexpr size_t num_neibs = 2;
	static constexpr size_t vd_s_close_part = 3;
	static constexpr size_t vd_s_sdf = 4;
	static constexpr size_t vd_s_sample = 5;
	static constexpr size_t interpol_coeff = 6;
	static constexpr size_t vd_in_sdf = 0;	//this is needed in the vd_in vector
	static constexpr size_t vd_in_close_part = 3;

	Redist_options redistOptions;
	particles_in_type &vd_in;
	particles_surface vd_s;
	double r_cutoff2;
	
	int return_sign(double phi)
	{
		if (phi > 0) return 1;
		if (phi < 0) return -1;
		return 0;
	}

	void detect_surface_particles()
	{
		auto NN = vd_in.getCellList(sqrt(r_cutoff2) + redistOptions.H);
		vd_in.updateCellList(NN);
		auto part = vd_in.getDomainIterator();

		while (part.isNext())
		{
			vect_dist_key_dx akey = part.get();
			int surfaceflag = 0;
			int sgn_a = return_sign(vd_in.template getProp<vd_in_sdf>(akey));
			Point<2,double> xa = vd_in.getPos(akey);
			int num_neibs_a = 0;
			double min_sdf = abs(vd_in.template getProp<vd_in_sdf>(akey));
			vect_dist_key_dx min_sdf_key = akey;
			vd_in.template getProp<vd_in_close_part>(akey) = 0;
			int isclose = 0;

			auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd_in.getPos(akey)));
			while (Np.isNext())
			{
				vect_dist_key_dx bkey = Np.get();
				int sgn_b = return_sign(vd_in.template getProp<vd_in_sdf>(bkey));
				Point<2, double> xb = vd_in.getPos(bkey);
	            Point<2,double> dr = xa - xb;
	            double r2 = norm2(dr);

	            //if ((sqrt(r2) < (1.3*redistOptions.H)) && !vd_in.template getProp<vd_in_close_part>(bkey) && (sgn_a != sgn_b)) isclose = 1;
	            if ((sqrt(r2) < (1.5*redistOptions.H)) && (sgn_a != sgn_b)) isclose = 1;
	            //int isclose = (abs(vd_in.template getProp<vd_in_sdf>(akey)) < (redistOptions.H/2.0 + 1e-8));

	            if (r2 < r_cutoff2)
	            {
	            	++num_neibs_a;

	            	//if (abs(abs(vd_in.template getProp<vd_in_sdf>(bkey)) - min_sdf) < 1e-8)
	            	if (((abs(vd_in.template getProp<vd_in_sdf>(bkey)) - min_sdf) < -1e-8) && !isclose)
	            	{
	            		//std::cout<<"diff in abs sdf values: "<<abs(vd_in.template getProp<vd_in_sdf>(bkey)) - min_sdf<<std::endl;
	            		min_sdf = abs(vd_in.template getProp<vd_in_sdf>(bkey));
	            		min_sdf_key = bkey;
	            	}
	            }

	            if ((r2 < (r_cutoff2 + 2*redistOptions.r_cutoff_factor*redistOptions.H + redistOptions.H*redistOptions.H)) && (sgn_a != sgn_b))
	            {
	            	surfaceflag = 1;
	            }
	            ++Np;
	            
			}

			if (surfaceflag)
			{

				vd_s.add();
				vd_s.getLastPos()[0] = vd_in.getPos(akey)[0];
				vd_s.getLastPos()[1] = vd_in.getPos(akey)[1];
				vd_s.template getLastProp<vd_s_sdf>() = vd_in.template getProp<vd_in_sdf>(akey);
				vd_s.template getLastProp<num_neibs>() = num_neibs_a;

				//vd_in.template getProp<vd_in_close_part>(min_sdf_key) = 1;

				//if (min_sdf_key.getKey() == akey.getKey())
				if (isclose)
				{
					vd_s.template getLastProp<vd_s_close_part>() = 1;
					vd_in.template getProp<vd_in_close_part>(akey) = 1;	// use this for the optimization step
					//std::cout<<"adding particles to interpolation problem..."<<std::endl;
				}
				else
				{
					vd_s.template getLastProp<vd_s_close_part>() = 0;
					vd_in.template getProp<vd_in_close_part>(akey) = 0;	// use this for the optimization step
				}
			}
			++part;
		}

	}
	
	void interpolate_sdf_field()
	{
		int verbose = 0;
		auto NN_s = vd_s.getCellList(sqrt(r_cutoff2));
		vd_s.updateCellList(NN_s);
		auto part = vd_s.getDomainIterator();

		while (part.isNext())
		{
			vect_dist_key_dx a = part.get();
			
			if (vd_s.template getProp<vd_s_close_part>(a) != 1)
			{
				//std::cout<<"vd_s contains particles that are not closest particles.."<<std::endl;
				++part;
				continue;
			}

			const int num_neibs_a = vd_s.template getProp<num_neibs>(a);
			Point<2, double> xa = vd_s.getPos(a);
			int neib = 0;

			MonomialBasis<2> m(4);

			if (redistOptions.polynomial_degree == "quadratic")
			{
				unsigned int degrees[2];
				degrees[0] = 1;
				degrees[1] = 1;
				m = MonomialBasis<2>(degrees, 1);
			}
			else if (redistOptions.polynomial_degree == "taylor4")
			{
				unsigned int degrees[2];
				degrees[0] = 1;
				degrees[1] = 1;
				m = MonomialBasis<2>(degrees, 3);
			}

			if(m.size() > num_neibs_a) std::cout<<"warning: less number of neighbours than required for interpolation"<<std::endl;

			// std::cout<<"m size"<<m.size()<<std::endl;
			VandermondeRowBuilder<2, double> vrb(m);

			EMatrix<double, Eigen::Dynamic,Eigen::Dynamic> V(num_neibs_a, m.size());
			EMatrix<double, Eigen::Dynamic, 1> phi(num_neibs_a, 1);

			auto Np = NN_s.template getNNIterator<NO_CHECK>(NN_s.getCell(vd_s.getPos(a)));
			while(Np.isNext())
			{
				vect_dist_key_dx b = Np.get();
				Point<2, double> xb = vd_s.getPos(b);
				Point<2, double> dr = xa - xb;
				double r2 = norm2(dr);

				if (r2 > r_cutoff2)
				{
					++Np;
					continue;
				}

				// debug given data
				//std::cout<<std::setprecision(15)<<xb[0]<<"\t"<<xb[1]<<"\t"<<vd.getProp<sdf>(b)<<std::endl;
				//xb[0] = 0.1;
				//xb[1] = 0.5;
				
				// Fill phi-vector from the right hand side
				phi[neib] = vd_s.template getProp<vd_s_sdf>(b);
				vrb.buildRow(V, neib, xb, 1.0);

				++neib;
				++Np;
			}

			EMatrix<double, Eigen::Dynamic, 1> c(m.size(), 1);

			c = V.completeOrthogonalDecomposition().solve(phi);

			for (int k = 0; k<m.size(); k++)
			{
				vd_s.template getProp<interpol_coeff>(a)[k] = c[k];
			}

			double dpdx;
			double dpdy;
			double grad_p_mag2;
			double incr;
			EMatrix<double, Eigen::Dynamic, 1> x(2, 1);
			x[0] = xa[0];
			x[1] = xa[1];
			double p = get_p(x, c, redistOptions.polynomial_degree);
			while (abs(p) > redistOptions.tolerance)
			//for(int k = 0; k < 3; k++)
			//while(incr > 0.01*redistOptions.H)
			{
				dpdx = get_dpdx(x, c, redistOptions.polynomial_degree);
				dpdy = get_dpdy(x, c, redistOptions.polynomial_degree);
				grad_p_mag2 = dpdx*dpdx + dpdy*dpdy;

				x[0] = x[0] - p*dpdx/grad_p_mag2;
				x[1] = x[1] - p*dpdy/grad_p_mag2;
				p = get_p(x, c, redistOptions.polynomial_degree);
				incr = sqrt(p*p*dpdx*dpdx/(grad_p_mag2*grad_p_mag2) + p*p*dpdy*dpdy/(grad_p_mag2*grad_p_mag2));
			}
			vd_s.template getProp<vd_s_sample>(a)[0] = x[0];
			vd_s.template getProp<vd_s_sample>(a)[1] = x[1];

			if (a.getKey() == 5424) verbose = 1;
			if (verbose)
			{
				EMatrix<double, Eigen::Dynamic, 1> xaa(2,1);
				xaa[0] = xa[0];
				xaa[1] = xa[1];
				//double curvature = get_dpdxdx(xaa, c) + get_dpdydy(xaa, c);
				// matlab debugging stuff:
				//std::cout<<std::setprecision(16)<<"A = ["<<V<<"];"<<std::endl;
				//std::cout<<"tempcond = cond(A);\ncondnumber = condnumber + tempcond;"<<std::endl;
				//std::cout<<"k = k + 1;\nif(tempcond>maxcond)\nmaxcond = tempcond;\nend"<<std::endl;
				//std::cout<<"if(tempcond<mincond)\nmincond = tempcond;\nend"<<std::endl;
				std::cout<<"number of coefficients is: "<<m.size()<<std::endl;
				std::cout<<"PHI: \n"<<phi<<std::endl;
				std::cout<<"num neibs: "<<num_neibs_a<<std::endl;
				std::cout<<"x: "<<xa[0]<<" "<<xa[1]<<std::endl;
				std::cout<<"COEFF: "<<c<<std::endl;
				//std::cout<<"KAPPA: "<<curvature<<std::endl;
				std::cout<<"Vmat\n"<<V<<std::endl;
				verbose = 0;
			}

			++part;
		}

	}
	
	void find_closest_point()
	{	// iterate over all particles, i.e. do closest point optimisation for all particles //
		auto NN_s = vd_s.getCellList(redistOptions.sampling_radius);
		auto part = vd_in.getDomainIterator(); //todo: include the sampling radius somewhere (bandwidth)
		vd_s.updateCellList(NN_s);
		int verbose = 0;
		int msize;
		if(redistOptions.polynomial_degree == "quadratic") msize = 6;
		else if(redistOptions.polynomial_degree == "bicubic") msize = 16;
		else if(redistOptions.polynomial_degree == "taylor4") msize = 15;
		else msize = 1; //todo: make this throw an error

		while (part.isNext())
		{
			vect_dist_key_dx a = part.get();

			// initialise all variables
			EMatrix<double, Eigen::Dynamic, 1> xa(2,1);
			xa[0] = vd_in.getPos(a)[0];
			xa[1] = vd_in.getPos(a)[1];
			double val;
			EMatrix<double, Eigen::Dynamic, 1> x(2,1);
			EMatrix<double, Eigen::Dynamic, 1> dx(3, 1);
			dx[0] = 1.0;
			dx[1] = 1.0;
			dx[2] = 1.0;
			double lambda = 0.0;
			double norm_dx = dx.norm();

			double p = 0;
			double dpdx = 0;
			double dpdy = 0;
			double dpdxdx = 0;
			double dpdydy = 0;
			double dpdxdy = 0;
			double dpdydx = 0;

			EMatrix<double, Eigen::Dynamic, 1> c(msize, 1);
			EMatrix<double, Eigen::Dynamic, 1> nabla_f(3, 1); // 3 = ndim + 1
			EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> H(3, 3);

			Point<2,double> xaa = vd_in.getPos(a);
			//auto Np = NN_s.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
			//vd_s.updateCellList(NN_s);
			auto Np = vd_s.getDomainIterator();
			double distance = 1000000.0;
			double dist_calc = 1000000000.0;
			decltype(a) b_min = a;

			while (Np.isNext())
			{
				vect_dist_key_dx b = Np.get();
				//Point<2,double> xbb = vd_s.getPos(b);
				Point<2, double> xbb = vd_s.template getProp<vd_s_sample>(b);

				if (!vd_s.template getProp<vd_s_close_part>(b))
				{
					//std::cout<<"vd s contains particles that are not closest in a surface neighborhood..."<<std::endl;
					++Np;
					continue;
				}
				dist_calc = abs(xbb.distance(xaa));
				//if(verbose)std::cout<<dist_calc<<std::endl;
				if (dist_calc < distance)
				{
					distance = dist_calc;
					b_min = b;
				}
				++Np;
			}

			// set x0 to the particle closest to the surface
			//x[0] = vd.getPos(vd_s.getProp<0>(b_min))[0];
			//x[1] = vd.getPos(vd_s.getProp<0>(b_min))[1];
			val = vd_s.template getProp<vd_s_sdf>(b_min);
			//x[0] = vd_s.getPos(b_min)[0];
			//x[1] = vd_s.getPos(b_min)[1];
			x[0] = vd_s.template getProp<vd_s_sample>(b_min)[0];
			x[1] = vd_s.template getProp<vd_s_sample>(b_min)[1];

			EMatrix<double, Eigen::Dynamic, 1> x00(2,1);
			x00[0] = x[0];
			x00[1] = x[1];
			EMatrix<double, Eigen::Dynamic, 1> x00x(2,1);
			x00x[0] = 0.0;
			x00x[1] = 0.0;

			// take the interpolation polynomial of the sample particle closest to the query particle
			for (int k = 0 ; k < msize ; k++)
			{
				//vd.getProp<interpol_coeff>(a)[k] = vd.getProp<interp_coeff>( vd_s.getProp<0>(b_min) )[k];
				c[k] = vd_s.template getProp<interpol_coeff>(b_min)[k];
			}

			val = get_p(x, c, redistOptions.polynomial_degree);
			//if (a.getKey() == 2620) verbose = 1;
			if (a.getKey() == 2425) verbose = 1;

			if(verbose)
				{
					std::cout<<std::setprecision(16)<<"VERBOSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:"<<std::endl;
					std::cout<<"x_poly: "<<vd_s.getPos(b_min)[0]<<", "<<vd_s.getPos(b_min)[1]<<"\nxa: "<<xa[0]<<", "<<xa[1]<<"\nx_0: "<<x[0]<<", "<<x[1]<<"\nc: "<<c<<std::endl;
					std::cout<<"phi(x_0): "<<val<<std::endl;
					std::cout<<"interpol_i(x_0) = "<<get_p(x, c, redistOptions.polynomial_degree)<<std::endl;
				}

			EMatrix<double, Eigen::Dynamic, 1> xax = x - xa;

			int k = 0;

			// barrier method initialisation
			double barrier = 0.0;
			double barrier2 = 0.0;
			double barrier3 = 0.0;

			// calculations needed specifically for k == 0
			p = get_p(x, c, redistOptions.polynomial_degree);
			dpdx = get_dpdx(x, c, redistOptions.polynomial_degree);
			dpdy = get_dpdy(x, c, redistOptions.polynomial_degree);
			lambda = (-xax[0]*dpdx - xax[1]*dpdy)/(dpdx*dpdx + dpdy*dpdy);

			nabla_f[0] = xax[0] + lambda*dpdx;
			nabla_f[1] = xax[1] + lambda*dpdy;
			nabla_f[2] = p;

			// possibly start off with a Newton-style projection towards the interface
			if (redistOptions.init_project)
			{
				while(abs(p) > redistOptions.tolerance)
				{
					dpdx = get_dpdx(x, c, redistOptions.polynomial_degree);
					dpdy = get_dpdy(x, c, redistOptions.polynomial_degree);
					p = get_p(x, c, redistOptions.polynomial_degree);
					double grad_p_mag2 = dpdx*dpdx + dpdy*dpdy;

					x[0] = x[0] - p*dpdx/grad_p_mag2;
					x[1] = x[1] - p*dpdy/grad_p_mag2;
				}
				dpdx = get_dpdx(x, c, redistOptions.polynomial_degree);
				dpdy = get_dpdy(x, c, redistOptions.polynomial_degree);
				xax = x - xa;
				lambda = (-xax[0]*dpdx - xax[1]*dpdy)/(dpdx*dpdx + dpdy*dpdy);
			}

			while((norm_dx > redistOptions.tolerance) && (k<redistOptions.max_iter))
			{
				// do optimisation //

				// gather required values //
				dpdxdx = get_dpdxdx(x, c, redistOptions.polynomial_degree);
				dpdydy = get_dpdydy(x, c, redistOptions.polynomial_degree);
				dpdxdy = get_dpdxdy(x, c, redistOptions.polynomial_degree);
				dpdydx = dpdxdy;

				// Assemble Hessian matrix, grad f has been computed at the end of the last iteration //
				H(0, 0) = 1 + lambda*dpdxdx;
				H(0, 1) = lambda*dpdydx;
				H(1, 0) = lambda*dpdxdy;
				H(1, 1) = 1 + lambda*dpdydy;
				H(0, 2) = dpdx;
				H(1, 2) = dpdy;
				H(2, 0) = dpdx;
				H(2, 1) = dpdy;
				H(2, 2) = 0;

				// barrier method //
				if (redistOptions.barrier_coefficient)
				{
					x00x = x - x00;
					barrier = 0.5*(x00x[0]*x00x[0] + x00x[1]*x00x[1] - r_cutoff2);
					//std::cout<<barrier<<std::endl;
					barrier2 = barrier*barrier;
					barrier3 = barrier2*barrier;

					// Add to gradient //
					nabla_f[0] -= redistOptions.barrier_coefficient*x00x[0]/barrier2;
					nabla_f[1] -= redistOptions.barrier_coefficient*x00x[1]/barrier2;
					// Add to Hessian matrix //
					H(0, 0) += redistOptions.barrier_coefficient*(-barrier + 2*x00x[0]*x00x[0])/barrier3;
					H(0, 1) += redistOptions.barrier_coefficient*2*x00x[0]*x00x[1]/barrier3;
					H(1, 0) += redistOptions.barrier_coefficient*2*x00x[0]*x00x[1]/barrier3;
					H(1, 1) += redistOptions.barrier_coefficient*(-barrier + 2*x00x[1]*x00x[1])/barrier3;
				}

				// compute and add increment // 
				dx = - H.inverse()*nabla_f;

				// Armijo rule
				if (redistOptions.armijo_tau)
				{
					double alpha_armijo = 1.0;
					double t_armijo = -redistOptions.armijo_c*(nabla_f[0]*dx[0] + nabla_f[1]*dx[1] + nabla_f[2]*dx[2]);
					// int k_armijo = 0;
					double f = 0.5*(x - xa).norm()*(x - xa).norm() + lambda*get_p(x, c, redistOptions.polynomial_degree);
					EMatrix<double, Eigen::Dynamic, 1> x_armijo(2,1);
					x_armijo[0] = x[0] + alpha_armijo*dx[0];
					x_armijo[1] = x[1] + alpha_armijo*dx[1];
					double lambda_armijo = lambda + alpha_armijo*dx[2];
					double f_armijo = 0.5*(xa - x_armijo).norm()*(xa - x_armijo).norm() + lambda_armijo*get_p(x_armijo, c, redistOptions.polynomial_degree);
					if (verbose) std::cout<<"start armijo\n"<<"x: "<<x[0]<<", "<<x[1]<<"\ndx: "<<dx[0]<<", "<<dx[1]<<"\nt: "<<t_armijo<<"\nf: "<<f<<"\nf_armijo: "<<f_armijo<<std::endl;
					while((f - f_armijo) < alpha_armijo*t_armijo)
					{
						alpha_armijo = redistOptions.armijo_tau*alpha_armijo;

						x_armijo[0] = x[0] + alpha_armijo*dx[0];
						x_armijo[1] = x[1] + alpha_armijo*dx[1];
						lambda_armijo = lambda + alpha_armijo*dx[2];
						f_armijo = 0.5*(xa - x_armijo).norm()*(xa - x_armijo).norm() + lambda_armijo*get_p(x_armijo, c, redistOptions.polynomial_degree);
						//std::cout<<alpha_armijo<<std::endl;
						if (verbose) std::cout<<"\narmijo++\n"<<"x: "<<x[0]<<", "<<x[1]<<"\ndx: "<<dx[0]<<", "<<dx[1]<<"\nt: "<<t_armijo<<"\nf: "<<f<<"\nf_armijo: "<<f_armijo<<"\ndiff f: "<<f-f_armijo<<"\n&&&&&&&&&&&&&&&&&&&&&&&&&"<<std::endl;
					}
					dx = alpha_armijo*dx;
				}

				// prevent Newton algorithm from leaving the support radius by projecting towards zero level-set
				if (redistOptions.fail_projection)
				{
					if ((x[0] + dx[0] - x00[0])*(x[0] + dx[0] - x00[0]) + (x[1] + dx[1] - x00[1])*(x[1] + dx[1] - x00[1]) > r_cutoff2)
					{
						EMatrix<double, Eigen::Dynamic, 1> x_interm(2,1);
						x_interm[0] = 1e16;
						x_interm[1] = 1e16;
						double dpdx_interm;
						double dpdy_interm;
						double grad_p_mag_interm2;
						double p_interm;
						int j = 0;

						while ((x_interm - x00).norm() > sqrt(r_cutoff2))
						{
							dx = redistOptions.support_prevent*dx;
							if(j == 0) x_interm = x;
							dpdx_interm = get_dpdx(x_interm, c, redistOptions.polynomial_degree);
							dpdy_interm = get_dpdy(x_interm, c, redistOptions.polynomial_degree);
							grad_p_mag_interm2 = dpdx_interm*dpdx_interm + dpdy_interm*dpdy_interm;
							p_interm = get_p(x_interm, c, redistOptions.polynomial_degree);
							x_interm[0] = x_interm[0] - redistOptions.support_prevent*p_interm*dpdx_interm/grad_p_mag_interm2;
							x_interm[1] = x_interm[1] - redistOptions.support_prevent*p_interm*dpdy_interm/grad_p_mag_interm2;
							++j;
							if (a.getKey() == 12)
							{
								//std::cout<<"p: "<<p_interm<<" r: "<<(x_interm - x00).norm()<<" j: "<<j<<std::endl;
							}
						}
						x[0] = 0.0;
						x[1] = 0.0;
						lambda = 0.0;
						dx[0] = x_interm[0];
						dx[1] = x_interm[1];
						xax = x_interm - xa;
						dx[2] = (-xax[0]*dpdx_interm - xax[1]*dpdy_interm)/grad_p_mag_interm2;
					}
				}

				// prevent Newton algorithm from leaving the support radius by scaling step size
				if(redistOptions.support_prevent)
				{
					while (((x[0] + dx[0] - x00[0])*(x[0] + dx[0] - x00[0]) + (x[1] + dx[1] - x00[1])*(x[1] + dx[1] - x00[1])) > r_cutoff2)
					{
						dx = redistOptions.support_prevent*dx;
						double rand_degree = 2*M_PI*rand()/(RAND_MAX + 1.);
						double dx_old = dx[0];
						double dy_old = dx[1];
						dx[0] = cos(rand_degree)*dx_old - sin(rand_degree)*dy_old;
						dx[1] = sin(rand_degree)*dx_old + cos(rand_degree)*dy_old;
					}
				}

				// prevent Newton algorithm from leaving the support radius by scaling step size
				while((dx[0]*dx[0] + dx[1]*dx[1]) > 0.25*r_cutoff2)
				{
					std::cout<<"invoked"<<std::endl;
					dx = 0.1*dx;
				}

				// apply increment and compute new iteration values
				x[0] = x[0] + dx[0];
				x[1] = x[1] + dx[1];
				lambda = lambda + dx[2];

				// prepare values for next iteration and update the exit criterion
				xax = x - xa;
				p = get_p(x, c, redistOptions.polynomial_degree);
				dpdx = get_dpdx(x, c, redistOptions.polynomial_degree);
				dpdy = get_dpdy(x, c, redistOptions.polynomial_degree);
				nabla_f[0] = xax[0] + lambda*dpdx;
				nabla_f[1] = xax[1] + lambda*dpdy;
				nabla_f[2] = p;

				//norm_dx = dx.norm(); // alternative criterion: incremental tolerance
				norm_dx = nabla_f.norm();

				++k;

				if(verbose)
				{
//					std::cout<<"dx: "<<dx[0]<<", "<<dx[1]<<std::endl;
//					std::cout<<"H:\n"<<H<<"\nH_inv:\n"<<H.inverse()<<std::endl;
//					std::cout<<"x:"<<x[0]<<", "<<x[1]<<"\nc:"<<std::endl;
//					std::cout<<c<<std::endl;
//					std::cout<<"dpdx: "<<dpdx<<std::endl;
//					std::cout<<"k = "<<k<<std::endl;
//					std::cout<<"x_k = "<<x[0]<<", "<<x[1]<<std::endl;
					std::cout<<x[0]<<", "<<x[1]<<", "<<norm_dx<<std::endl;
				}
			}
	//		debug optimisation

	//		std::cout<<"old sdf: "<<vd.getProp<sdf>(a)<<std::endl;
			vd_in.template getProp<vd_in_sdf>(a) = return_sign(vd_in.template getProp<vd_in_sdf>(a))*xax.norm();
			//std::cout<<"new sdf: "<<vd.getProp<sdf>(a)<<std::endl;
			if(verbose)
			{
				std::cout<<"x_final: "<<x[0]<<", "<<x[1]<<std::endl;
				std::cout<<"p(x_final) :"<<get_p(x, c, redistOptions.polynomial_degree)<<std::endl;
				std::cout<<"nabla p(x_final)"<<get_dpdx(x, c, redistOptions.polynomial_degree)<<", "<<get_dpdy(x, c, redistOptions.polynomial_degree)<<std::endl;
				std::cout<<"lambda: "<<lambda<<std::endl;
			}
			if (k == redistOptions.max_iter) std::cout<<"Warning: Newton algorithm has reached maximum number of iterations, does not converge for particle "<<a.getKey()<<std::endl;
			//std::cout<<"c:\n"<<vd.getProp<interpol_coeff>(a)<<"\nmin_sdf: "<<vd.getProp<min_sdf>(a)<<" , min_sdf_x: "<<vd.getProp<min_sdf_x>(a)<<std::endl;
			verbose = 0;
			if (((abs(x00[0] - x[0]))>sqrt(r_cutoff2))||((abs(x00[1] - x[1]))>sqrt(r_cutoff2)))
			{
				std::cout<<"straying out of local neighborhood.."<<std::endl;
				std::cout<<"computed sdf: "<<vd_in.template getProp<vd_in_sdf>(a)<<std::endl;
				std::cout<<"analytical value: "<<vd_in.template getProp<8>(a)<<", diff: "<<abs(vd_in.template getProp<vd_in_sdf>(a) - vd_in.template getProp<8>(a))<<std::endl;
			}

			++part;
		}
		
	}

	//todo: template these
	inline double get_p(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c, std::string polynomialdegree)
	{
		const double x = xvector[0];
		const double y = xvector[1];
		if (polynomialdegree == "bicubic")
		{
				return(c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x + c[4]*y + c[5]*x*y + c[6]*x*x*y + c[7]*x*x*x*y + c[8]*y*y + c[9]*x*y*y + c[10]*x*x*y*y + c[11]*x*x*x*y*y +c[12]*y*y*y + c[13]*x*y*y*y + c[14]*x*x*y*y*y + c[15]*x*x*x*y*y*y);
		}
		if (polynomialdegree == "quadratic")
		{
				return(c[0] + c[1]*x + c[2]*x*x + c[3]*y + c[4]*x*y + c[5]*y*y);
		}
		if (polynomialdegree == "taylor4")
		{
				return(c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x + c[4]*x*x*x*x + c[5]*y + c[6]*x*y + c[7]*x*x*y + c[8]*x*x*x*y + c[9]*y*y + c[10]*y*y*x + c[11]*x*x*y*y + c[12]*y*y*y + c[13]*y*y*y*x + c[14]*y*y*y*y);
		}
	}

	inline double get_dpdx(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c, std::string polynomialdegree)
	{	
		const double x = xvector[0];
		const double y = xvector[1];
		if (polynomialdegree == "bicubic")
		{
				return(c[1] + 2*c[2]*x + 3*c[3]*x*x + c[5]*y + 2*c[6]*x*y + 3*c[7]*x*x*y + c[9]*y*y + 2*c[10]*x*y*y + 3*c[11]*x*x*y*y + c[13]*y*y*y + 2*c[14]*x*y*y*y + 3*c[15]*x*x*y*y*y);
		}
		if (polynomialdegree == "quadratic")
		{
				return(c[1] + 2*c[2]*x + c[4]*y);
		}
		if (polynomialdegree == "taylor4")
		{
				return(c[1] + 2*c[2]*x + 3*c[3]*x*x + 4*c[4]*x*x*x + c[6]*y + 2*c[7]*x*y + 3*c[8]*x*x*y + c[10]*y*y + 2*c[11]*x*y*y + c[13]*y*y*y);
		}
	}

	inline double get_dpdy(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c, std::string polynomialdegree)
	{	
		const double x = xvector[0];
		const double y = xvector[1];
		if (polynomialdegree == "bicubic")
		{
				return(c[4] + c[5]*x + c[6]*x*x + c[7]*x*x*x + 2*c[8]*y + 2*c[9]*x*y + 2*c[10]*x*x*y + 2*c[11]*x*x*x*y + 3*c[12]*y*y + 3*c[13]*x*y*y + 3*c[14]*x*x*y*y + 3*c[15]*x*x*x*y*y);
		}
		if (polynomialdegree == "quadratic")
		{
				return(c[3] + c[4]*x + 2*c[5]*y);
		}
		if (polynomialdegree == "taylor4")
		{
				return(c[5] + c[6]*x + c[7]*x*x + c[8]*x*x*x + 2*c[9]*y + 2*c[10]*x*y + 2*c[11]*x*x*y + 3*c[12]*y*y + 3*c[13]*y*y*x + 4*c[14]*y*y*y);
		}
	}

	inline double get_dpdxdx(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c, std::string polynomialdegree)
	{
		const double x = xvector[0];
		const double y = xvector[1];
		if (polynomialdegree == "bicubic")
		{
				return(2*c[2] + 6*c[3]*x + 2*c[6]*y + 6*c[7]*x*y + 2*c[10]*y*y + 6*c[11]*y*y*x + 2*c[14]*y*y*y + 6*c[15]*y*y*y*x);
		}
		if (polynomialdegree == "quadratic")
		{
				return(2*c[2]);
		}
		if (polynomialdegree == "taylor4")
		{
				return(2*c[2] + 6*c[3]*x + 12*c[4]*x*x + 2*c[7]*y + 6*c[8]*x*y + 2*c[11]*y*y);
		}
	}

	inline double get_dpdydy(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c, std::string polynomialdegree)
	{
		const double x = xvector[0];
		const double y = xvector[1];
		if (polynomialdegree == "bicubic")
		{
				return(2*c[8] + 2*c[9]*x + 2*c[10]*x*x + 2*c[11]*x*x*x + 6*c[12]*y + 6*c[13]*x*y + 6*c[14]*x*x*y + 6*c[15]*x*x*x*y);
		}
		if (polynomialdegree == "quadratic")
		{
				return(2*c[5]);
		}
		if (polynomialdegree == "taylor4")
		{
				return(2*c[9] + 2*c[10]*x + 2*c[11]*x*x + 6*c[12]*y + 6*c[13]*x*y + 12*c[14]*y*y);
		}
	}

	inline double get_dpdxdy(EMatrix<double, Eigen::Dynamic, 1> xvector, EMatrix<double, Eigen::Dynamic, 1> c, std::string polynomialdegree)
	{
		const double x = xvector[0];
		const double y = xvector[1];
		if (polynomialdegree == "bicubic")
		{
				return(c[5] + 2*c[6]*x + 3*c[7]*x*x + 2*c[9]*y + 4*c[10]*x*y + 6*c[11]*x*x*y + 3*c[13]*y*y + 6*c[14]*x*y*y + 9*c[15]*x*x*y*y);
		}
		if (polynomialdegree == "quadratic")
		{
				return(c[4]);
		}
		if (polynomialdegree == "taylor4")
		{
				return(c[6] + 2*c[7]*x + 3*c[8]*x*x + 2*c[10]*y +4*c[11]*x*y + 3*c[13]*y*y);
		}
	}
	
};
