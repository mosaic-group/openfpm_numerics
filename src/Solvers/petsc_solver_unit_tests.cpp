/*
 * petsc_solver_unit_tests.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_UNIT_TESTS_CPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_UNIT_TESTS_CPP_

#include "config.h"
#ifdef HAVE_PETSC

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "petsc_solver_report_unit_tests.hpp"

#include "Grid/grid_dist_id.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/Laplacian.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "Solvers/petsc_solver.hpp"


BOOST_AUTO_TEST_SUITE( mg_solvers_test )


BOOST_AUTO_TEST_CASE( laplacian_3D_int_zero_mg )
{
	constexpr unsigned int phi = 0;
	typedef Field<phi,poisson_nn_helm> phi_f;

	Vcluster<> & v_cl = create_vcluster();
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

	petsc_solver<double> solver;
	solver.setSolver(KSPBCGS);
	solver.setAbsTol(0.01);
	solver.setMaxIter(500);

	////////

	solver.setPreconditioner(PCHYPRE_BOOMERAMG);
	solver.setPreconditionerAMG_nl(6);
	solver.setPreconditionerAMG_maxit(3);
	solver.setPreconditionerAMG_relax("SOR/Jacobi");
	solver.setPreconditionerAMG_cycleType("V",6,6);
	solver.setPreconditionerAMG_coarsen("HMIS");
	solver.setPreconditionerAMG_coarsenNodalType(0);

	timer tm_solve;
	tm_solve.start();
	auto x_ = solver.solve(fd.getA(),fd.getB());
	tm_solve.stop();

	fd.template copy<phi>(x_,psi);
	psi.write("AMG_psi");

	#ifdef HAVE_OSX
	bool check = compare("AMG_psi_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk","test/AMG_psi_" + std::to_string(v_cl.getProcessUnitID()) + "_test_osx.vtk");
	#elif __GNUC__ == 6
	bool check = compare("AMG_psi_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk","test/AMG_psi_" + std::to_string(v_cl.getProcessUnitID()) + "_test_GCC6.vtk");
	#elif __GNUC__ == 4
	bool check = compare("AMG_psi_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk","test/AMG_psi_" + std::to_string(v_cl.getProcessUnitID()) + "_test_GCC4.vtk");
	#else
	bool check = true;
	#endif

	BOOST_REQUIRE_EQUAL(check,true);


	// Resolve
	timer tm_solve2;
	tm_solve2.start();
	auto x2_ = solver.solve(fd.getB());
	tm_solve2.stop();

	fd.template copy<phi>(x_,psi2);
	psi2.write("AMG_psi2");

	#ifdef HAVE_OSX
	check = compare("AMG_psi2_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk","test/AMG_psi2_" + std::to_string(v_cl.getProcessUnitID()) + "_test_osx.vtk");
	#elif __GNUC__ == 6
	check = compare("AMG_psi2_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk","test/AMG_psi2_" + std::to_string(v_cl.getProcessUnitID()) + "_test_GCC6.vtk");
	#elif __GNUC__ == 4
    	check = compare("AMG_psi2_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk","test/AMG_psi2_" + std::to_string(v_cl.getProcessUnitID()) + "_test_GCC4.vtk");
	#else
	check = true;
	#endif
	BOOST_REQUIRE_EQUAL(check,true);
}

BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_UNIT_TESTS_CPP_ */
