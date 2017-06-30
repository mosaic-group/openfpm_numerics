/*
 * petsc_solver_report.hpp
 *
 *  Created on: Jun 22, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_AMG_REPORT_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_AMG_REPORT_HPP_

#include <fstream>

/*! \brief It contain information about the performance of the AMG
 *
 *
 */
struct AMG_time_err_coars
{
	//! time to setup
	double time_setup;

	//! time to solve
	double time_solve;

	//! residual norm
	double residual;

	//! assign a score to the solver
	double score;
};

/*! \brief Class to test AMG solvers
 *
 *
 */
class petsc_AMG_report
{
	//! Number of sweeps for test
	size_t num_sweep_test = 10;

	//! Target accuracy for the preconditioner
	double target_accuracy = 1;

	//! It contains the performance of several AMG methods coarsener
	openfpm::vector<AMG_time_err_coars> perf_amg_coars;

	//! It contain the performance of several AMG methods sweep configuration
	openfpm::vector<AMG_time_err_coars> perf_amg_sweep_sym;

	//! It contain the performance of several AMG methods sweep configuration
	openfpm::vector<AMG_time_err_coars> perf_amg_sweep_asym;

	//! List of methods
	openfpm::vector<std::string> method;

	//! List of tested cycles
	openfpm::vector<size_t> v_cycle_tested_sym;

	//! List of tested cycles in case of asymmetric
	openfpm::vector<std::pair<size_t,size_t>> v_cycle_tested_asym;

	//! starting point of each coarsening test for the asymmetric test
	openfpm::vector<size_t> asym_tests;

	/*! \brief benchmark the selected setting for the preconditioner
	 *
	 * \param A sparse matrix
	 * \param b right-hand-side
	 * \param solver to use for benchmarking
	 * \param perf_amg performance of the AMG
	 *
	 */
	void benchmark(SparseMatrix<double,int,PETSC_BASE> & A,
			       const Vector<double,PETSC_BASE> & b,
				   petsc_solver<double> & solver,
				   openfpm::vector<AMG_time_err_coars> & perf_amg)
	{
		timer tm_solve;
		tm_solve.start();
		auto x_ = solver.solve(A,b);
		tm_solve.stop();

		solError serr = solver.get_residual_error(A,x_,b);

		// In order to measure the time to solve we solve again the system

		timer tm_solve2;
		tm_solve2.start();
		solver.solve(x_,b);
		tm_solve2.stop();

		double time1 = tm_solve.getwct();
		double time2 = tm_solve2.getwct();

		Vcluster & v_cl = create_vcluster();
		v_cl.max(time1);
		v_cl.max(time2);

		// Save the result
		AMG_time_err_coars tmp;
		tmp.time_setup = time1 - time2;
		tmp.time_solve = time2;
		tmp.residual = serr.err_inf;

		perf_amg.add(tmp);

		if (v_cl.getProcessUnitID() == 0)
			std::cout << "Time1: " << time1 << "    Time2: " << time2 << std::endl;
	}

	/*! \brief test the corasener for this problem
	 *
	 * \param A matrix to invert
	 * \param b right-hand-side
	 *
	 */
	void test_coarsener(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & b)
	{
		Vcluster & v_cl = create_vcluster();

		petsc_solver<double> pts;

		petsc_solver<double> solver;
		solver.setSolver(KSPRICHARDSON);
		solver.setAbsTol(0.01);
		solver.setMaxIter(1);

		//////// Firts we check several coarsener ////

		solver.setPreconditioner(PCHYPRE_BOOMERAMG);
		solver.setPreconditionerAMG_nl(6);
		solver.setPreconditionerAMG_maxit(3);
		solver.setPreconditionerAMG_relax("SOR/Jacobi");
		solver.setPreconditionerAMG_cycleType("V",6,6);
		solver.setPreconditionerAMG_coarsenNodalType(0);
		solver.log_monitor();


		// Reset the Matrix A
		A.getMatrixTriplets();
		if (v_cl.getProcessUnitID() == 0)	{std::cout << "Benchmarking HMIS coarsener" << std::endl;}
		solver.setPreconditionerAMG_coarsen("HMIS");
		benchmark(A,b,solver,perf_amg_coars);
		method.add("HMIS");

		// Reset the Matrix A
		A.getMatrixTriplets();
		if (v_cl.getProcessUnitID() == 0)	{std::cout << "Benchmarking PMIS coarsener" << std::endl;}
		solver.setPreconditionerAMG_coarsen("PMIS");
		benchmark(A,b,solver,perf_amg_coars);
		method.add("PMIS");

		// Reset the Matrix A
		A.getMatrixTriplets();
		if (v_cl.getProcessUnitID() == 0)	{std::cout << "Benchmarking Falgout coarsener" << std::endl;}
		solver.setPreconditionerAMG_coarsen("Falgout");
		benchmark(A,b,solver,perf_amg_coars);
		method.add("Falgout");

		// Reset the Matrix A
		A.getMatrixTriplets();
		if (v_cl.getProcessUnitID() == 0)	{std::cout << "Benchmarking modifiedRuge-Stueben coarsener" << std::endl;}
		solver.setPreconditionerAMG_coarsen("modifiedRuge-Stueben");
		benchmark(A,b,solver,perf_amg_coars);
		method.add("modifiedRuge-Stueben");

		// Reset the Matrix A
		A.getMatrixTriplets();
		if (v_cl.getProcessUnitID() == 0)	{std::cout << "Benchmarking Ruge-Stueben coarsener" << std::endl;}
		solver.setPreconditionerAMG_coarsen("Ruge-Stueben");
		benchmark(A,b,solver,perf_amg_coars);
		method.add("Ruge-Stueben");

		// Reset the Matrix A
		A.getMatrixTriplets();
		if (v_cl.getProcessUnitID() == 0)	{std::cout << "Benchmarking CLJP coarsener" << std::endl;}
		solver.setPreconditionerAMG_coarsen("CLJP");
		benchmark(A,b,solver,perf_amg_coars);
		method.add("CLJP");

	}

	/*! \brief Score the solver
	 *
	 * \param t_solve time to solve
	 * \param t_m accuracy reache by the solver
	 *
	 *
	 */
	double score_solver(double t_solve, double t_m)
	{
		return 1.0/t_solve * std::min(1.0,target_accuracy/t_m);
	}

	/*! \brief Write the report for coarsening
	 *
	 * \param cg Google chart
	 *
	 */
	void write_report_coars(GoogleChart & cg)
	{
		std::cout << "Composing coarsening report" << std::endl;

		openfpm::vector<std::string> x;
		openfpm::vector<openfpm::vector<double>> y;
		openfpm::vector<std::string> yn;
		openfpm::vector<double> xd;

		for (size_t i = 0 ; i < perf_amg_coars.size() ; i++)
			x.add(method.get(i));

		// Each colum can have multiple data set (in this case 4 dataset)
		// Each dataset can have a name
		yn.add("Norm infinity");
		yn.add("time to setup");
		yn.add("time to solve");

		// Each colums can have multiple data-set
		for (size_t i = 0 ; i < perf_amg_coars.size() ; i++)
			y.add({perf_amg_coars.get(i).residual,perf_amg_coars.get(i).time_setup,perf_amg_coars.get(i).time_solve});

		// Google charts options
		GCoptions options;

		options.title = std::string("Overview of different coarsener on a V cycle with 6 relaxation SOR/Jacobi sweeps for each level up and down");
		options.yAxis = std::string("Error/time");
		options.xAxis = std::string("Method");
		options.stype = std::string("bars");
		options.more = std::string("vAxis: {scaleType: 'log'}");

		std::stringstream ss;

		cg.addHTML("<h2>Coarsener</h2>");
		cg.AddHistGraph(x,y,yn,options);
	}

	/*! \brief Write the report for coarsening
	 *
	 * \param cg Google chart
	 * \param coars type of coarsening
	 *
	 */
	void write_report_cycle(GoogleChart & cg, openfpm::vector<std::string> & coars)
	{
		cg.addHTML("<h2>V cycle sweep</h2>");
		for(size_t k = 0 ; k < perf_amg_sweep_sym.size() / num_sweep_test ; k++)
		{
			openfpm::vector<double> x;
			openfpm::vector<openfpm::vector<double>> y;
			openfpm::vector<std::string> yn;
			openfpm::vector<double> xd;

			x.clear();
			for (size_t i = 0 ; i < num_sweep_test ; i++)
				x.add(v_cycle_tested_sym.get(i));

			// Each colum can have multiple data set (in this case 4 dataset)
			// Each dataset can have a name
			yn.add("error");
			yn.add("time");
			yn.add("score");

			// Each colums can have multiple data-set
			for (size_t i = 0 ; i < num_sweep_test ; i++)
			{
				double score = score_solver(perf_amg_sweep_sym.get(k*num_sweep_test+i).time_solve,perf_amg_sweep_sym.get(k*num_sweep_test+i).residual);
				y.add({perf_amg_sweep_sym.get(k*num_sweep_test+i).residual,perf_amg_sweep_sym.get(k*num_sweep_test+i).time_solve,score});
			}

			// Google charts options
			GCoptions options;

			options.title = std::string("Coarsener: " + coars.get(k)) + std::string(": Overview of V cycle with n sweeps symmetrically up and down");
			options.yAxis = std::string("error/time/score");
			options.xAxis = std::string("sweeps");
			options.stype = std::string("bars");
			options.more = std::string("vAxis: {scaleType: 'log'}");

			std::stringstream ss;

			cg.AddHistGraph(x,y,yn,options);
		}
	}

	/*! \brief Write the report for coarsening
	 *
	 * \param cg Google chart
	 * \param coars type of coarsening
	 *
	 */
	void write_report_cycle_asym(GoogleChart & cg, openfpm::vector<std::string> & coars)
	{
		cg.addHTML("<h2>V cycle sweep asymmetric test</h2>");
		for(size_t k = 0 ; k < coars.size() ; k++)
		{
			openfpm::vector<std::string> x;
			openfpm::vector<openfpm::vector<double>> y;
			openfpm::vector<std::string> yn;
			openfpm::vector<double> xd;

			size_t num_sweep_test = asym_tests.get(k+1) - asym_tests.get(k);
			size_t sweep_start = asym_tests.get(k);

			x.clear();
			for (size_t i = 0 ; i < num_sweep_test ; i++)
			{x.add(std::string("(" + std::to_string(v_cycle_tested_asym.get(i+sweep_start).first) + "," + std::to_string(v_cycle_tested_asym.get(i+sweep_start).second) + ")"));}

			// Each colum can have multiple data set (in this case 4 dataset)
			// Each dataset can have a name
			yn.add("error");
			yn.add("time");
			yn.add("score");

			// Each colums can have multiple data-set
			for (size_t i = 0 ; i < num_sweep_test ; i++)
			{
				double score = score_solver(perf_amg_sweep_asym.get(sweep_start+i).time_solve,perf_amg_sweep_asym.get(sweep_start+i).residual);
				y.add({perf_amg_sweep_asym.get(sweep_start+i).residual,perf_amg_sweep_asym.get(sweep_start+i).time_solve,score});
			}

			// Google charts options
			GCoptions options;

			options.title = std::string("Coarsener: " + coars.get(k)) + std::string(": Overview of V cycle with n sweeps asymmetrically distributed");
			options.yAxis = std::string("error/time/score");
			options.xAxis = std::string("sweeps");
			options.stype = std::string("bars");
			options.more = std::string("vAxis: {scaleType: 'log'}");

			std::stringstream ss;

			cg.AddHistGraph(x,y,yn,options);
		}
	}

	/*! \brief test the corasener for this problem
	 *
	 * \param A matrix to invert
	 * \param b right-hand-side
	 * \param coarsener_to_test coarsener to test
	 *
	 */
	void test_cycle_type(SparseMatrix<double,int,PETSC_BASE> & A,
			             const Vector<double,PETSC_BASE> & b,
						 openfpm::vector<std::string> & coarsener_to_test)
	{
		Vcluster & v_cl = create_vcluster();

		petsc_solver<double> pts;

		petsc_solver<double> solver;
		solver.setSolver(KSPRICHARDSON);
		solver.setAbsTol(0.01);
		solver.setMaxIter(1);

		//////// we check different sweep configuration for the
		//////// fastest and most accurate coarsener

		for (size_t k = 0 ; k < coarsener_to_test.size(); k++)
		{
			for (size_t i = 1 ; i <= num_sweep_test ; i++)
			{

				solver.setPreconditioner(PCHYPRE_BOOMERAMG);
				solver.setPreconditionerAMG_nl(6);
				solver.setPreconditionerAMG_maxit(3);
				solver.setPreconditionerAMG_relax("SOR/Jacobi");
				solver.setPreconditionerAMG_cycleType("V",i,i);
				solver.setPreconditionerAMG_coarsenNodalType(0);
				solver.setPreconditionerAMG_coarsen(coarsener_to_test.get(k));
				solver.log_monitor();

				A.getMatrixTriplets();
				solver.setPreconditionerAMG_coarsen(coarsener_to_test.get(k).c_str());
				benchmark(A,b,solver,perf_amg_sweep_sym);

				double score = score_solver(perf_amg_sweep_sym.get(k*num_sweep_test+(i-1)).time_solve,perf_amg_sweep_sym.get(k*num_sweep_test+(i-1)).residual);
				perf_amg_sweep_sym.last().score = score;

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Coarsener: " << coarsener_to_test.get(k) << "   Sweep: " << i << "  Score: " << score << std::endl;

				v_cycle_tested_sym.add(i);
			}
		}

	}

	/*! \brief test best asymmetric combination of sweeps
	 *
	 * \param A matrix to invert
	 * \param b right-hand-side
	 * \param ids_ts set of sweep configuration to test
	 * \param coarsener_to_test coarsener to test
	 *
	 */
	void test_cycle_asym(SparseMatrix<double,int,PETSC_BASE> & A,
			             const Vector<double,PETSC_BASE> & b,
						 openfpm::vector<size_t> & ids_ts,
						 openfpm::vector<std::string> & coarsener_to_test)
	{
		Vcluster & v_cl = create_vcluster();

		petsc_solver<double> pts;

		petsc_solver<double> solver;
		solver.setSolver(KSPRICHARDSON);
		solver.setAbsTol(0.01);
		solver.setMaxIter(1);

		//////// we check different sweep configuration for the
		//////// fastest and most accurate coarsener

		for (size_t k = 0 ; k < coarsener_to_test.size(); k++)
		{
			asym_tests.add(perf_amg_sweep_asym.size());
			size_t num_tot_sweep = 2*ids_ts.get(k);
			for (size_t i = 0 ; i <= num_tot_sweep ; i++)
			{

				solver.setPreconditioner(PCHYPRE_BOOMERAMG);
				solver.setPreconditionerAMG_nl(6);
				solver.setPreconditionerAMG_maxit(3);
				solver.setPreconditionerAMG_relax("SOR/Jacobi");
				solver.setPreconditionerAMG_cycleType("V",i,num_tot_sweep-i);
				solver.setPreconditionerAMG_coarsenNodalType(0);
				solver.setPreconditionerAMG_coarsen(coarsener_to_test.get(k));
				solver.log_monitor();

				A.getMatrixTriplets();
				solver.setPreconditionerAMG_coarsen(coarsener_to_test.get(k).c_str());

				benchmark(A,b,solver,perf_amg_sweep_asym);

				double score = score_solver(perf_amg_sweep_asym.last().time_solve,perf_amg_sweep_asym.last().residual);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Coarsener: " << coarsener_to_test.get(k) << "   Sweep: (" << i << "," << num_tot_sweep-i << ")  Score: " << score << std::endl;

				v_cycle_tested_asym.add(std::pair<size_t,size_t>(i,num_tot_sweep-i));
			}

		}

		asym_tests.add(perf_amg_sweep_asym.size());
	}

	/*! Write the report on file
	 *
	 * \param coars set of coarsener tested
	 *
	 */
	void write_report(openfpm::vector<std::string> & coars)
	{
		if (create_vcluster().getProcessUnitID() != 0)
			return;

		GoogleChart cg;
		write_report_coars(cg);

		// Write V cycles performances
		write_report_cycle(cg,coars);

		write_report_cycle_asym(cg,coars);

		cg.write("gc_AMG_preconditioners.html");
	}

	/*! \brief take the most accurate and the fastest AMG
	 *
	 * \param coars vector with the best coarsener
	 *
	 */
	void best_coarsener(openfpm::vector<std::string> & coars)
	{
		int fastest = 0;
		double best_time = perf_amg_coars.get(0).time_solve;
		int accurate = 0;
		double best_accuracy = perf_amg_coars.get(0).residual;

		for (size_t i = 1 ; i < perf_amg_coars.size() ; i++)
		{
			if (perf_amg_coars.get(i).time_solve < best_time)
			{
				fastest = i;
				best_time = perf_amg_coars.get(i).time_solve;
			}

			if (perf_amg_coars.get(i).residual < best_accuracy)
			{
				accurate = i;
				best_accuracy = perf_amg_coars.get(i).residual;
			}
		}

		coars.add(method.get(fastest));
		coars.add(method.get(accurate));
	}

	/*! \brief Return the best scoring solver
	 *
	 * \param perf where to search for the best score
	 * \param mn number of method tested
	 * \param sw_accu best sweep for the best in accuracy
	 * \param sw_fast best sweep for the best in performance
	 *
	 */
	void best_score(openfpm::vector<AMG_time_err_coars> & perf,
					openfpm::vector<size_t> & sw_optimal,
					size_t nm)
	{
		size_t page = perf.size() / nm;

		for (size_t k = 0 ; k < nm ; k++)
		{
			int score_id = page*k;
			double best_score = perf.get(score_id).score;

			for (size_t i = page*k ; i < page*k+page ; i++)
			{
				if (perf.get(i).score > best_score)
				{
					score_id = i;
					best_score = perf.get(i).score;
				}
			}

			sw_optimal.add(v_cycle_tested_sym.get(score_id));
		}
	}

public:

	/*! \brief Set the target accuracy to score the AMG solver
	 *
	 * The score is calculated as \f$  \frac{1}{t_s} max(1,t_a/t_m) \f$
	 *
	 * where \f$ t_s \f$ is the time to solve the system \f$ t_a \f$ is the
	 * target accuracy \f$ t_m \f$ is the time of the method
	 *
	 * \param t_a target accuracy
	 *
	 */
	void setTergetAMGAccuracy(double t_a)
	{
		target_accuracy = t_a;
	}

	/*! \brief Set the target accuracy to score the AMG solver
	 *
	 * Set the maximum number of sweep for testing
	 *
	 * \param max_sw max number of sweep
	 *
	 */
	void setMaxSweep(int max_sw)
	{
		num_sweep_test = max_sw;
	}

	/*! \brief Try to use AMG pre-conditioner and check how they
	 *         they perform
	 *
	 * \param A Sparse matrix
	 * \param b right-hand side
	 *
	 *
	 */
	void try_solve(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & b)
	{
		test_coarsener(A,b);

		openfpm::vector<std::string> coa_methods;
		best_coarsener(coa_methods);

		test_cycle_type(A,b,coa_methods);
		openfpm::vector<size_t> sweep_optimal;
		best_score(perf_amg_sweep_sym,sweep_optimal,coa_methods.size());

		test_cycle_asym(A,b,sweep_optimal,coa_methods);

		write_report(coa_methods);
	}
};


#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_AMG_REPORT_HPP_ */
