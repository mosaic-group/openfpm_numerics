/*
 * petsc_solver_report_unit_tests.hpp
 *
 *  Created on: Jun 23, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_REPORT_UNIT_TESTS_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_REPORT_UNIT_TESTS_HPP_


#include "Grid/grid_dist_id.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/Laplacian.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "Solvers/petsc_solver.hpp"
#include "petsc_solver_AMG_report.hpp"

struct poisson_nn_helm
{
	//! Number of dimansions of the problem
    static const unsigned int dims = 3;

    //! We have only one scalar unknown
    static const unsigned int nvar = 1;

    //! specify the boundary conditions
    static const bool boundary[];

    //! type of the space
    typedef float stype;

    //! Grid that store the solution
    typedef grid_dist_id<3,float,aggregate<float>> b_grid;

    //! type of sparse grid that store the Matrix A
    typedef SparseMatrix<double,int,PETSC_BASE> SparseMatrix_type;

    //! type of vector that store the solution
    typedef Vector<double,PETSC_BASE> Vector_type;

    //! When the grid is not staggered use NORMAL_GRID
    static const int grid_type = NORMAL_GRID;
};

const bool poisson_nn_helm::boundary[] = {PERIODIC,PERIODIC,PERIODIC};

BOOST_AUTO_TEST_SUITE( mg_solvers_report_test )


BOOST_AUTO_TEST_CASE( laplacian_3D_int_zero_mg_report )
{
	constexpr unsigned int phi = 0;
	typedef Field<phi,poisson_nn_helm> phi_f;

	Vcluster & v_cl = create_vcluster();
	if (v_cl.getProcessingUnits() != 3)
		return;

    Ghost<3,long int> g(2);
	Ghost<3,long int> stencil_max(2);
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	periodicity<3> p({PERIODIC,PERIODIC,PERIODIC});

	typedef Lap<phi_f,poisson_nn_helm,CENTRAL_SYM> poisson;
	grid_dist_id<3,float,aggregate<float>> psi({64,64,64},domain,g,p);
	grid_dist_id<3,float,aggregate<float>> psi2(psi.getDecomposition(),{64,64,64},g);

	// Fill B

	size_t center_x = psi.size(0) / 2;
	size_t center_y = psi.size(1) / 2;
	size_t center_z = psi.size(2) / 2;
	auto it = psi.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();
		auto gkey = it.getGKey(key);

		float sx = (float)(gkey.get(0))-center_x;
		float sy = (float)(gkey.get(1))-center_y;
		float sz = (float)(gkey.get(2))-center_z;

		float gs = 100.0*exp(-((sx*sx)+(sy*sy)+(sz*sz))/100.0);

		psi.get<0>(key) = sin(2*M_PI*sx/psi.size(0))*sin(2*M_PI*sy/psi.size(1))*sin(2*M_PI*sz/psi.size(2))*gs;
		psi2.get<0>(key) = sin(2*M_PI*sx/psi.size(0))*sin(2*M_PI*sy/psi.size(1))*sin(2*M_PI*sz/psi.size(2))*gs;

		++it;
	}

	FDScheme<poisson_nn_helm> fd(stencil_max, domain, psi);

	fd.template impose_dit<0>(poisson(),psi,psi.getDomainIterator());

	petsc_AMG_report rep;

	rep.setTergetAMGAccuracy(10e-4);

	rep.try_solve(fd.getA(),fd.getB());
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_REPORT_UNIT_TESTS_HPP_ */
