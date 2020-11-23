/*
 * petsc_solver.hpp
 *
 *  Created on: Apr 26, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_

#include "config.h"

#ifdef HAVE_PETSC

#include "Vector/Vector.hpp"
#include <petscksp.h>
#include <petsctime.h>
#include "Plot/GoogleChart.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include <sstream>
#include <iomanip>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

#define UMFPACK_NONE 0
#define REL_DOWN 1
#define REL_UP 2
#define REL_ALL 3

#define PCHYPRE_BOOMERAMG "petsc_boomeramg"

enum AMG_type
{
	NONE_AMG,
	HYPRE_AMG,
	PETSC_AMG,
	TRILINOS_ML
};


/*! \brief In case T does not match the PETSC precision compilation create a
 *         stub structure
 *
 * \param T precision
 *
 */
template<typename T>
class petsc_solver
{
public:

	/*! \brief Solve the linear system.
	 *
	 * In this case just return an error
	 *
	 * \param A sparse matrix
	 * \param b right-hand-side
	 *
	 * \return the solution
	 *
	 */
	template<typename impl> static Vector<T> solve(const SparseMatrix<T,impl> & A, const Vector<T> & b)
	{
		std::cerr << "Error Petsc only suppor double precision" << "/n";
	}
};

#define SOLVER_NOOPTION 0
#define SOLVER_PRINT_RESIDUAL_NORM_INFINITY 1
#define SOLVER_PRINT_DETERMINANT 2

//! It contain statistic of the error of the calculated solution
struct solError
{
	//! infinity norm of the error
	PetscReal err_inf;

	//! L1 norm of the error
	PetscReal err_norm;

	//! Number of iterations
	PetscInt its;
};


/*! \brief This class is able to do Matrix inversion in parallel with PETSC solvers
 *
 * # Example of use of this class #
 *
 * \snippet eq_unit_test.hpp lid-driven cavity 2D
 *
 */
template<>
class petsc_solver<double>
{
	//! contain the infinity norm of the residual at each iteration
	struct itError
	{
		//! Iteration
		PetscInt it;

		//! error norm at iteration it
		PetscReal err_norm;
	};

	/*! \brief It contain the benchmark information for each solver
	 *
	 *
	 */
	struct solv_bench_info
	{
		//! Solver Krylov name
		std::string method;

		//! Method name (short form)
		std::string smethod;

		//! time to converge in milliseconds
		double time;

		//! Solution error
		solError err;

		//! Convergence per iteration
		openfpm::vector<itError> res;
	};

	//! indicate if the preconditioner is set
	bool is_preconditioner_set = false;

	//! KSP Maximum number of iterations
	PetscInt maxits;

	//! Main parallel solver
	KSP ksp;

	//! Temporal variable used for calculation of static members
	size_t tmp;

	//! The full set of solvers
	openfpm::vector<std::string> solvs;


	//! It contain the solver benchmark results
	openfpm::vector<solv_bench_info> bench;

	//! Type of the algebraic multi-grid preconditioner
	AMG_type atype = NONE_AMG;

	//! Block size
	int block_sz = 0;

	/*! \brief Calculate the residual error at time t for one method
	 *
	 * \param t time
	 * \param slv solver
	 *
	 * \return the residual
	 *
	 */
	static double calculate_it(double t, solv_bench_info & slv)
	{
		double s_int = slv.time / slv.res.size();

		// Calculate the discrete point in time
		size_t pt = std::floor(t / s_int);

		if (pt < slv.res.size())
			return slv.res.get(pt).err_norm;

		return slv.res.last().err_norm;
	}

	/*! \brief Here we write the benchmark report
	 *
	 *
	 */
	void write_bench_report()
	{
		if (create_vcluster().getProcessUnitID() != 0)
			return;

		openfpm::vector<std::string> x;
		openfpm::vector<openfpm::vector<double>> y;
		openfpm::vector<std::string> yn;
		openfpm::vector<double> xd;

		for (size_t i = 0 ; i < bench.size() ; i++)
			x.add(bench.get(i).smethod);

		// Each colum can have multiple data set (in this case 4 dataset)
		// Each dataset can have a name
		yn.add("Norm infinity");
		yn.add("Norm average");
		yn.add("Number of iterations");

		// Each colums can have multiple data-set

		for (size_t i = 0 ; i < bench.size() ; i++)
			y.add({bench.get(i).err.err_inf,bench.get(i).err.err_norm,(double)bench.get(i).err.its});

		// Google charts options
		GCoptions options;

		options.title = std::string("Residual error");
		options.yAxis = std::string("Error");
		options.xAxis = std::string("Method");
		options.stype = std::string("bars");
		options.more = std::string("vAxis: {scaleType: 'log'}");

		GoogleChart cg;

		std::stringstream ss;

		cg.addHTML("<h2>Robustness</h2>");
		cg.AddHistGraph(x,y,yn,options);
		cg.addHTML("<h2>Speed</h2>");

		y.clear();
		yn.clear();

		yn.add("Time in millseconds");
		options.title = std::string("Time");

		for (size_t i = 0 ; i < bench.size() ; i++)
			y.add({bench.get(i).time});

		cg.AddHistGraph(x,y,yn,options);

		// Convergence in time

		x.clear();
		y.clear();
		yn.clear();
		xd.clear();

		// Get the maximum in time across all the solvers

		double max_time = 0.0;

		for (size_t i = 0 ; i < bench.size() ; i++)
		{
			if (bench.get(i).time > max_time)	{max_time = bench.get(i).time;}
			yn.add(bench.get(i).smethod);
		}

		size_t n_int = maxits;

		// calculate dt

		double dt = max_time / n_int;

		// Add

		// For each solver we have a convergence plot

		for (double t = dt ; t <= max_time + 0.05 * max_time ; t += dt)
		{
			y.add();
			xd.add(t);
			for (size_t j = 0 ; j < bench.size() ; j++)
				y.last().add(calculate_it(t,bench.get(j)));
		}

		std::stringstream ss2;
		ss2 << "<h2>Convergence with time</h2><br>" << std::endl;

		options.title = std::string("Residual norm infinity");
		options.yAxis = std::string("residual");
		options.xAxis = std::string("Time in milliseconds");
		options.lineWidth = 1.0;
		options.more = std::string("vAxis: {scaleType: 'log'},hAxis: {scaleType: 'log'}");

		cg.addHTML(ss2.str());
		cg.AddLinesGraph(xd,y,yn,options);

		ss << "<h2>Solvers Tested</h2><br><br>" << std::endl;
		for (size_t i = 0 ; i < bench.size() ; i++)
			ss << bench.get(i).smethod << " = " << bench.get(i).method << "<br>" << std::endl;

		cg.addHTML(ss.str());

		cg.write("gc_solver.html");
	}

	/*! \brief It set up the solver based on the provided options
	 *
	 * \param A_ Matrix
	 * \param x_ solution
	 * \param b_ right-hand-side
	 *
	 */
	void pre_solve_impl(const Mat & A_, const Vec & b_, Vec & x_)
	{
		PETSC_SAFE_CALL(KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED));

		if (atype == HYPRE_AMG)
		{
			PC pc;

			// We set the pre-conditioner
			PETSC_SAFE_CALL(KSPGetPC(ksp,&pc));

			PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
			PCFactorSetShiftAmount(pc, PETSC_DECIDE);
		}
	}

	/*! \brief Print a progress bar on standard out
	 *
	 *
	 */
	void print_progress_bar()
	{
		auto & v_cl = create_vcluster();

		if (v_cl.getProcessUnitID() == 0)
			std::cout << "-----------------------25-----------------------50-----------------------75----------------------100" << std::endl;
	}


	/*! \brief procedure print the progress of the solver in benchmark mode
	 *
	 * \param ksp Solver
	 * \param it Iteration number
	 * \param res resudual
	 * \param data custom pointer to data
	 *
	 * \return always zero
	 *
	 */
	static PetscErrorCode monitor_progress_residual(KSP ksp,PetscInt it,PetscReal res,void* data)
	{
		petsc_solver<double> * pts = (petsc_solver *)data;

		pts->progress(it);

		itError erri;
		erri.it = it;

		Mat A;
		Vec Br,v,w;

		PETSC_SAFE_CALL(KSPGetOperators(ksp,&A,NULL));
		PETSC_SAFE_CALL(MatCreateVecs(A,&w,&v));
		PETSC_SAFE_CALL(KSPBuildResidual(ksp,v,w,&Br));
		PETSC_SAFE_CALL(KSPBuildResidual(ksp,v,w,&Br));
		PETSC_SAFE_CALL(VecNorm(Br,NORM_INFINITY,&erri.err_norm));
		itError err_fill;

		size_t old_size = pts->bench.last().res.size();
		pts->bench.last().res.resize(it+1);

		if (old_size > 0)
			err_fill = pts->bench.last().res.get(old_size-1);
		else
			err_fill = erri;

		for (long int i = old_size ; i < (long int)it ; i++)
			pts->bench.last().res.get(i) = err_fill;

	    // Add the error per iteration
		pts->bench.last().res.get(it) = erri;

		return 0;
	}

	/*! \brief This function print an "*" showing the progress of the solvers
	 *
	 * \param it iteration number
	 *
	 */
	void progress(PetscInt it)
	{
		PetscInt n_max_it;

		PETSC_SAFE_CALL(KSPGetTolerances(ksp,NULL,NULL,NULL,&n_max_it));

		auto & v_cl = create_vcluster();

		if (std::floor(it * 100.0 / n_max_it) > tmp)
		{
			size_t diff = it * 100.0 / n_max_it - tmp;
			tmp = it * 100.0 / n_max_it;
			if (v_cl.getProcessUnitID() == 0)
			{
				for (size_t k = 0 ; k <  diff ; k++)	{std::cout << "*";}
				std::cout << std::flush;
			}
		}
	}

	/*! \brief Print statistic about the solution error and method used
	 *
	 * \param err structure that contain the solution errors
	 *
	 */
	void print_stat(solError & err)
	{
		auto & v_cl = create_vcluster();

		if (v_cl.getProcessUnitID() == 0)
		{
			KSPType typ;
			KSPGetType(ksp,&typ);

			std::cout << "Method: " << typ << " " << " pre-conditoner: " << PCJACOBI << "  iterations: " << err.its << std::endl;
			std::cout << "Norm of error: " << err.err_norm << "   Norm infinity: " << err.err_inf << std::endl;
		}
	}


	/*! \brief It convert the KSP type into a human read-able string
	 *
	 * \param solv solver (short form)
	 *
	 * \return the name of the solver in long form
	 *
	 */
	std::string to_string_method(const std::string & solv)
	{
		if (solv == std::string(KSPBCGS))
		{
			return std::string("BiCGStab (Stabilized version of BiConjugate Gradient Squared)");
		}
		else if (solv == std::string(KSPIBCGS))
		{
			return std::string("IBiCGStab (Improved Stabilized version of BiConjugate Gradient Squared)");
		}
		else if (solv == std::string(KSPFBCGS))
		{
			return std::string("FBiCGStab (Flexible Stabilized version of BiConjugate Gradient Squared)");
		}
		else if (solv == std::string(KSPFBCGSR))
		{
			return std::string("FBiCGStab-R (Modified Flexible Stabilized version of BiConjugate Gradient Squared)");
		}
		else if (solv == std::string(KSPBCGSL))
		{
			return std::string("BiCGStab(L) (Stabilized version of BiConjugate Gradient Squared with L search directions)");
		}
		else if (solv == std::string(KSPGMRES))
		{
			return std::string("GMRES(M) (Generalized Minimal Residual method with restart)");
		}
		else if (solv == std::string(KSPFGMRES))
		{
			return std::string("FGMRES(M) (Flexible Generalized Minimal Residual with restart)");
		}
		else if (solv == std::string(KSPLGMRES))
		{
			return std::string("LGMRES(M) (Augment Generalized Minimal Residual method with restart)");
		}
		else if (solv == std::string(KSPPGMRES))
		{
			return std::string("PGMRES(M) (Pipelined Generalized Minimal Residual method)");
		}
		else if (solv == std::string(KSPGCR))
		{
			return std::string("GCR (Generalized Conjugate Residual)");
		}

		return std::string("UNKNOWN");
	}

	/*! \brief Allocate a new benchmark slot for a method
	 *
	 * \param str Method name
	 *
	 */
	void new_bench(const std::string & str)
	{
		// Add a new benchmark result
		bench.add();
		bench.last().method = to_string_method(str);
		bench.last().smethod = std::string(str);
	}

	/*! \brief Copy the solution if better
	 *
	 * \param res Residual of the solution
	 * \param sol solution
	 * \param best_res residual of the best solution
	 * \param best_sol best solution
	 *
	 */
	void copy_if_better(double res, Vec & sol, double & best_res, Vec & best_sol)
	{
		if (res < best_res)
		{
			VecCopy(sol,best_sol);
			best_res = res;
		}
	}

	/*! \brief Try to solve the system x=inv(A)*b using all the Krylov solvers with simple Jacobi pre-conditioner
	 *
	 *   It try to solve the system using JACOBI pre-conditioner and all the Krylov solvers available at the end it write
	 *   a report
	 *
	 * \param A_ Matrix
	 * \param b_ vector of coefficents
	 * \param x_ solution
	 *
	 */
	void try_solve_simple(const Mat & A_, const Vec & b_, Vec & x_)
	{
		Vec best_sol;
		PETSC_SAFE_CALL(VecDuplicate(x_,&best_sol));

		// Best residual
		double best_res = std::numeric_limits<double>::max();

		// Create a new VCluster
		auto & v_cl = create_vcluster();

		destroyKSP();

		// for each solver
		for (size_t i = 0 ; i < solvs.size() ; i++)
		{
			initKSPForTest();

			// Here we solve without preconditioner
			PC pc;

			// We set the pre-conditioner to none
			PETSC_SAFE_CALL(KSPGetPC(ksp,&pc));
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_type",PCNONE));

			setSolver(solvs.get(i).c_str());

			// Setup for BCGSL, GMRES
			if (solvs.get(i) == std::string(KSPBCGSL))
			{
				// we try from 2 to 6 as search direction
				for (size_t j = 2 ; j < 6 ; j++)
				{
					new_bench(solvs.get(i));

					if (v_cl.getProcessUnitID() == 0)
						std::cout << "L = " << j << std::endl;
					bench.last().smethod += std::string("(") + std::to_string(j) + std::string(")");
					searchDirections(j);
					bench_solve_simple(A_,b_,x_,bench.last());

					copy_if_better(bench.last().err.err_inf,x_,best_res,best_sol);
				}
			}
			else if (solvs.get(i) == std::string(KSPGMRES) ||
					 solvs.get(i) == std::string(std::string(KSPFGMRES)) ||
					 solvs.get(i) == std::string(std::string(KSPLGMRES)) )
			{
				// we try from 2 to 6 as search direction
				for (size_t j = 50 ; j < 300 ; j += 50)
				{
					// Add a new benchmark result
					new_bench(solvs.get(i));

					if (v_cl.getProcessUnitID() == 0)
						std::cout << "Restart = " << j << std::endl;
					bench.last().smethod += std::string("(") + std::to_string(j) + std::string(")");
					setRestart(j);
					bench_solve_simple(A_,b_,x_,bench.last());

					copy_if_better(bench.last().err.err_inf,x_,best_res,best_sol);
				}
			}
			else
			{
				// Add a new benchmark result
				new_bench(solvs.get(i));

				bench_solve_simple(A_,b_,x_,bench.last());

				copy_if_better(bench.last().err.err_inf,x_,best_res,best_sol);
			}

			destroyKSP();
		}

		write_bench_report();

		// Copy the best solution to x
		PETSC_SAFE_CALL(VecCopy(best_sol,x_));
		PETSC_SAFE_CALL(VecDestroy(&best_sol));
	}


	/*! \brief Benchmark solve simple solving x=inv(A)*b
	 *
	 * \param A_ Matrix A
	 * \param b_ vector b
	 * \param x_ solution x
	 * \param bench structure that store the benchmark information
	 *
	 */
	void bench_solve_simple(const Mat & A_, const Vec & b_, Vec & x_, solv_bench_info & bench)
	{
		// timer for benchmark
		timer t;
		t.start();
		// Enable progress monitor
		tmp = 0;
		PETSC_SAFE_CALL(KSPMonitorCancel(ksp));
		PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor_progress_residual,this,NULL));
		print_progress_bar();
		solve_simple(A_,b_,x_);

		// New line
		if (create_vcluster().getProcessUnitID() == 0)
			std::cout << std::endl;

		t.stop();
		bench.time = t.getwct() * 1000.0;

		// calculate error statistic about the solution and print
		solError err = statSolutionError(A_,b_,x_,ksp);
		print_stat(err);

		bench.err = err;

		// New line
		if (create_vcluster().getProcessUnitID() == 0)
			std::cout << std::endl;
	}

	/*! \brief solve simple use a Krylov solver + Simple selected Parallel Pre-conditioner
	 *
	 * \param A_ SparseMatrix
	 * \param b_ right-hand-side
	 * \param x_ solution
	 *
	 */
	void solve_simple(const Mat & A_, const Vec & b_, Vec & x_)
	{
//		PETSC_SAFE_CALL(KSPSetType(ksp,s_type));

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;

		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));

		// We set the Matrix operators
		PETSC_SAFE_CALL(KSPSetOperators(ksp,A_,A_));

		// if we are on on best solve set-up a monitor function

		PETSC_SAFE_CALL(KSPSetFromOptions(ksp));
		PETSC_SAFE_CALL(KSPSetUp(ksp));

		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));
	}

	/*! \brief solve simple use a Krylov solver + Simple selected Parallel Pre-conditioner
	 *
	 * \param b_ right hand side
	 * \param x_ solution
	 *
	 */
	void solve_simple(const Vec & b_, Vec & x_)
	{
		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));
	}

	/*! \brief Calculate statistic on the error solution
	 *
	 * \param A_ Matrix of the system
	 * \param b_ Right hand side of the matrix
	 * \param x_ Solution
	 * \param ksp Krylov solver
	 *
	 * \return the solution error
	 *
	 */
	static solError statSolutionError(const Mat & A_, const Vec & b_, Vec & x_, KSP ksp)
	{
		solError err;

		err = getSolNormError(A_,b_,x_);

		PetscInt its;
		PETSC_SAFE_CALL(KSPGetIterationNumber(ksp,&its));
		err.its = its;

		return err;
	}


	/*! \brief initialize the KSP object
	 *
	 *
	 */
	void initKSP()
	{
		PETSC_SAFE_CALL(KSPCreate(PETSC_COMM_WORLD,&ksp));
	}

	/*! \brief initialize the KSP object for solver testing
	 *
	 *
	 */
	void initKSPForTest()
	{
		PETSC_SAFE_CALL(KSPCreate(PETSC_COMM_WORLD,&ksp));

		setMaxIter(maxits);
	}

	/*! \brief Destroy the KSP object
	 *
	 *
	 */
	void destroyKSP()
	{
		KSPDestroy(&ksp);
	}

	/*! \brief Return the norm error of the solution
	 *
	 * \param x_ the solution
	 * \param b_ the right-hand-side
	 * \param ksp Krylov solver
	 *
	 * \return the solution error
	 *
	 */
	static solError getSolNormError(const Vec & b_, const Vec & x_,KSP ksp)
	{
		Mat A_;
		Mat P_;
		KSPGetOperators(ksp,&A_,&P_);

		return getSolNormError(A_,b_,x_);
	}

	/*! \brief Return the norm error of the solution
	 *
	 * \param A_ the matrix that identity the linear system
	 * \param x_ the solution
	 * \param b_ the right-hand-side
	 *
	 * \return the solution error
	 *
	 */
	static solError getSolNormError(const Mat & A_, const Vec & b_, const Vec & x_)
	{
		PetscScalar neg_one = -1.0;

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;

		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));

		/*
		      Here we calculate the residual error
	    */
		PetscReal norm;
		PetscReal norm_inf;

		// Get a vector r for the residual
		Vector<double,PETSC_BASE> r(row,row_loc);
		Vec & r_ = r.getVec();

		PETSC_SAFE_CALL(MatMult(A_,x_,r_));
		PETSC_SAFE_CALL(VecAXPY(r_,neg_one,b_));

		PETSC_SAFE_CALL(VecNorm(r_,NORM_1,&norm));
		PETSC_SAFE_CALL(VecNorm(r_,NORM_INFINITY,&norm_inf));

		solError err;
		err.err_norm = norm / col;
		err.err_inf = norm_inf;
		err.its = 0;

		return err;
	}

public:

	//! Type of the solution object
	typedef Vector<double,PETSC_BASE> return_type;

	~petsc_solver()
	{
		PETSC_SAFE_CALL(KSPDestroy(&ksp));
	}

	petsc_solver()
	:maxits(300),tmp(0)
	{
		initKSP();

		// Add the solvers

		solvs.add(std::string(KSPBCGS));
//		solvs.add(std::string(KSPIBCGS)); <--- Produce invalid argument
//		solvs.add(std::string(KSPFBCGS));
//		solvs.add(std::string(KSPFBCGSR)); <--- Nan production problems
		solvs.add(std::string(KSPBCGSL));
		solvs.add(std::string(KSPGMRES));
		solvs.add(std::string(KSPFGMRES));
		solvs.add(std::string(KSPLGMRES));
//		solvs.add(std::string(KSPPGMRES)); <--- Take forever
//		solvs.add(std::string(KSPGCR));

		setSolver(KSPGMRES);
	}

	/*! \brief Print the preconditioner used by the solver
	 *
	 *
	 */
	void print_preconditioner()
	{
		auto & v_cl = create_vcluster();

		if (v_cl.getProcessUnitID() == 0)
		{
			PC pc;
			PCType type_pc;

			// We set the pre-conditioner
			PETSC_SAFE_CALL(KSPGetPC(ksp,&pc));
			PETSC_SAFE_CALL(PCGetType(pc,&type_pc));

			std::cout << "The precoditioner is: " << type_pc << std::endl;
		}
	}

	/*! \brief Add a test solver
	 *
	 * The try solve function use the most robust solvers in PETSC, if you want to add
	 * additionally solver like KSPIBCGS,KSPFBCGSR,KSPPGMRES, use addTestSolver(std::string(KSPIBCGS))
	 *
	 * \param solver additional solver solver to test
	 *
	 */
	void addTestSolver(std::string & solver)
	{
		solvs.add(solver);
	}

	/*! \brief Remove a test solver
	 *
	 * The try solve function use the most robust solvers in PETSC, if you want to remove
	 * a solver like use removeTestSolver(std::string(KSPIBCGS))
	 *
	 * \param solver remove solver to test
	 *
	 */
	void removeTestSolver(const std::string & solver)
	{
		// search for the test solver and remove it

		for (size_t i = 0 ; i < solvs.size() ; i++)
		{
			if (solvs.get(i) == solver)
				return;
		}
	}


	/*! \brief Set the Petsc solver
	 *
	 * \see KSPType in PETSC manual for a list of all PETSC solvers
	 *
	 */
	void log_monitor()
	{
		PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-ksp_monitor",0));
		PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_print_statistics","2"));
	}

	/*! \brief Set the Petsc solver
	 *
	 * \see KSPType in PETSC manual for a list of all PETSC solvers
	 *
	 * \param type petsc solver type
	 *
	 */
	void setSolver(KSPType type)
	{
		PetscOptionsSetValue(NULL,"-ksp_type",type);
	}

	/*! \brief Set the relative tolerance as stop criteria
	 *
	 * \see PETSC manual KSPSetTolerances for an explanation
	 *
	 * \param rtol_ Relative tolerance
	 *
	 */
	void setRelTol(PetscReal rtol_)
	{
		PetscReal rtol;
		PetscReal abstol;
		PetscReal dtol;
		PetscInt maxits;

		PETSC_SAFE_CALL(KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits));
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol_,abstol,dtol,maxits));
	}

	/*! \brief Set the absolute tolerance as stop criteria
	 *
	 * \see PETSC manual KSPSetTolerances for an explanation
	 *
	 * \param abstol_ Absolute tolerance
	 *
	 */
	void setAbsTol(PetscReal abstol_)
	{
		PetscReal rtol;
		PetscReal abstol;
		PetscReal dtol;
		PetscInt maxits;

		PETSC_SAFE_CALL(KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits));
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol,abstol_,dtol,maxits));
	}

	/*! \brief Set the divergence tolerance
	 *
	 * \see PETSC manual KSPSetTolerances for an explanation
	 *
	 * \param dtol_ tolerance
	 *
	 */
	void setDivTol(PetscReal dtol_)
	{
		PetscReal rtol;
		PetscReal abstol;
		PetscReal dtol;
		PetscInt maxits;

		PETSC_SAFE_CALL(KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits));
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol,abstol,dtol_,maxits));
	}

	/*! \brief Set the maximum number of iteration for Krylov solvers
	 *
	 * \param n maximum number of iterations
	 *
	 */
	void setMaxIter(PetscInt n)
	{
		PetscReal rtol;
		PetscReal abstol;
		PetscReal dtol;
		PetscInt maxits;

		PETSC_SAFE_CALL(KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits));
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol,abstol,dtol,n));

		this->maxits = n;
	}

	/*! For the BiCGStab(L) it define the number of search directions
	 *
	 * BiCG Methods can fail for base break-down (or near break-down). BiCGStab(L) try
	 * to avoid this problem. Such method has a parameter L (2 by default) Bigger is L
	 * more the system will try to avoid the breakdown
	 *
	 * \param l Increasing L should reduce the probability of failure of the solver because of break-down of the base
	 *
	 */
	void searchDirections(PetscInt l)
	{
		PetscOptionsSetValue(NULL,"-ksp_bcgsl_ell",std::to_string(l).c_str());
	}

	/*! \brief For GMRES based method,  the number of Krylov directions to orthogonalize against
	 *
	 * \param n number of directions
	 *
	 */
	void setRestart(PetscInt n)
	{
		PetscOptionsSetValue(NULL,"-ksp_gmres_restart",std::to_string(n).c_str());
	}

	/*! \brief Set the preconditioner of the linear solver
	 *
	 * The preconditoner that can be set are the same as PETSC
	 *
	 * For a full list please visit
	 * Please visit: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType
	 *
	 * An exception is PCHYPRE with BOOMERAMG in this case use
	 * PCHYPRE_BOOMERAMG. Many preconditioners has default values, but many
	 * times the default values are not good. Here we list some interesting case
	 *
	 * ## Algebraic-multi-grid based preconditioners##
	 *
	 * Parameters for this type of preconditioner can be set using
	 * setPreconditionerAMG_* functions.
	 *
	 * * Number of levels set by setPreconditionerAMG_nl
	 * * Maximum number of cycles setPreconditionerAMG_maxit
	 * * Smooth or relax method setPreconditionerAMG_smooth,setPreconditionerAMG_relax
	 * * Cycle type and number of sweep (relaxation steps) when going up or down setPreconditionerAMG_cycleType
	 * * Coarsening method setPreconditionerAMG_coarsen
	 * * interpolation schemes setPreconditionerAMG_interp
	 *
	 *
	 * \param type of the preconditioner
	 *
	 */
	void setPreconditioner(PCType type)
	{
		is_preconditioner_set = true;

		if (std::string(type) == PCHYPRE_BOOMERAMG)
		{
			PC pc;

			// We set the pre-conditioner
			PETSC_SAFE_CALL(KSPGetPC(ksp,&pc));

			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_type",PCHYPRE));

		    PETSC_SAFE_CALL(PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO));
		    PETSC_SAFE_CALL(PCFactorSetShiftAmount(pc, PETSC_DECIDE));
#ifdef PETSC_HAVE_HYPRE
		    PETSC_SAFE_CALL(PCHYPRESetType(pc, "boomeramg"));
#endif
		    atype = HYPRE_AMG;
		}
		else
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_type",type));
		}
	}

	/*! \brief Set the number of levels for the algebraic-multigrid preconditioner
	 *
	 * In case you select an algebraic preconditioner like PCHYPRE or PCGAMG you can
	 * set the number of levels using this function
	 *
	 * \param nl number of levels
	 *
	 */
	void setPreconditionerAMG_nl(int nl)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_max_levels",std::to_string(nl).c_str()));
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}

	/*! \brief Set the maximum number of V or W cycle for algebraic-multi-grid
	 *
	 * \param nit number of levels
	 *
	 */
	void setPreconditionerAMG_maxit(int nit)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_max_iter",std::to_string(nit).c_str()));
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}


	/*! \brief Set the relaxation method for the algebraic-multi-grid preconditioner
	 *
	 * Possible values for relazation can be
	 * "Jacobi","sequential-Gauss-Seidel","seqboundary-Gauss-Seidel",
	 * "SOR/Jacobi","backward-SOR/Jacobi",hybrid chaotic Gauss-Seidel (works only with OpenMP),
	 * "symmetric-SOR/Jacobi","l1scaled-SOR/Jacobi","Gaussian-elimination","CG","Chebyshev",
	 * "FCF-Jacobi","l1scaled-Jacobi"
	 *
	 *
	 *
	 *
	 * Every smooth operator can have additional parameters to be set in order to
	 * correctly work.
	 *
	 *
	 * \param type of relax method
	 * \param k where is applied REL_UP,REL_DOWN,REL_ALL (default is all)
	 *
	 */
	void setPreconditionerAMG_relax(const std::string & type, int k = REL_ALL)
	{
		if (atype == HYPRE_AMG)
		{
			if (k == REL_ALL)
			{PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_all",type.c_str()));}
			else if (k == REL_UP)
			{PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_up",type.c_str()));}
			else if (k == REL_DOWN)
			{PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_relax_type_down",type.c_str()));}
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}

	/*! \brief It set the type of cycle and optionally the number of sweep
	 *
	 * This function set the cycle type for the multigrid methods.
	 * Possible values are:
	 * * V cycle
	 * * W cycle
	 *
	 * Optionally you can set the number of sweep or relaxation steps
	 * on each grid when going up and when going down
	 *
	 * \param cycle_type cycle type
	 * \param sweep_up
	 * \param sweep_dw
	 * \param sweep_crs speep at the coarse level
	 *
	 */
	void setPreconditionerAMG_cycleType(const std::string & cycle_type, int sweep_up = -1, int sweep_dw = -1, int sweep_crs = -1)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_cycle_type",cycle_type.c_str()));

			if (sweep_dw != -1)
			{PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_down",std::to_string(sweep_up).c_str()));}

			if (sweep_up != -1)
			{PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_up",std::to_string(sweep_dw).c_str()));}

			if (sweep_crs != -1)
			{PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_coarse",std::to_string(sweep_crs).c_str()));}
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}

	/*! \brief Set the coarsening method for the algebraic-multi-grid preconditioner
	 *
	 * Possible values can be
	 * "CLJP","Ruge-Stueben","modifiedRuge-Stueben","Falgout", "PMIS", "HMIS"
	 *
	 * \warning in case of big problem use PMIS or HMIS otherwise the AMG can hang (take a lot of time) in setup
	 *
	 * \param type of the preconditioner smoothing operator
	 *
	 */
	void setPreconditionerAMG_coarsen(const std::string & type)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_coarsen_type",type.c_str()));
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}

	/*! \brief Set the interpolation method for the algebraic-multi-grid preconditioner
	 *
	 * Possible values can be
	 * "classical", "direct", "multipass", "multipass-wts", "ext+i",
     * "ext+i-cc", "standard", "standard-wts", "FF", "FF1"
	 *
	 *
	 * \param type of the interpolation scheme
	 *
	 */
	void setPreconditionerAMG_interp(const std::string & type)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_interp_type",type.c_str()));
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}

	/*! \brief Set the block coarsening norm type.
	 *
	 * The use of this function make sanse if you specify
	 * the degree of freedom for each  node using setBlockSize
	 *
	 * * 0 (default each variable treat independently)
	 * * 1 Frobenius norm
     * * 2 sum of absolute values of elements in each block
     * * 3 largest element in each block (not absolute value)
     * * 4 row-sum norm
	 * * 6 sum of all values in each block
	 *
	 * In case the matrix represent a system of equations in general
	 *
	 * \param norm type
	 *
	 */
	void setPreconditionerAMG_coarsenNodalType(int norm)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_nodal_coarsen",std::to_string(norm).c_str()));
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}

	/*! \brief Indicate the number of levels in the ILU(k) for the Euclid smoother
	 *
	 * \param k number of levels for the Euclid smoother
	 *
	 */
	void setPreconditionerAMG_interp_eu_level(int k)
	{
		if (atype == HYPRE_AMG)
		{
			PETSC_SAFE_CALL(PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_eu_level",std::to_string(k).c_str()));
		}
		else
		{
			std::cout << __FILE__ << ":" << __LINE__ << "Warning: the selected preconditioner does not support this option" << std::endl;
		}
	}



	/*! \brief Set how many degree of freedom each node has
	 *
	 * In case you are solving a system of equations this function,
	 * help in setting the degree of freedom each grid point has.
	 * Setting this parameter to the number of variables for each
	 * grid point it should improve the convergenve of the solvers
	 * in particular using algebraic-multi-grid
	 *
	 * \param block_sz number of degree of freedom
	 *
	 */
	void setBlockSize(int block_sz)
	{
		this->block_sz = block_sz;
	}

	/*! \brief Here we invert the matrix and solve the system
	 *
	 *  \warning umfpack is not a parallel solver, this function work only with one processor
	 *
	 *  \note if you want to use umfpack in a NON parallel, but on a distributed data, use solve with triplet
	 *
	 *	\tparam impl Implementation of the SparseMatrix
	 *
	 * \param A sparse matrix
	 * \param b vector
	 * \param initial_guess true if x has the initial guess
	 *
	 * \return the solution
	 *
	 */
	Vector<double,PETSC_BASE> solve(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & b, bool initial_guess = false)
	{
		Mat & A_ = A.getMat();
		const Vec & b_ = b.getVec();

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;

		PETSC_SAFE_CALL(KSPSetInitialGuessNonzero(ksp,PETSC_FALSE));
		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));

		Vector<double,PETSC_BASE> x(row,row_loc);
		Vec & x_ = x.getVec();

		pre_solve_impl(A_,b_,x_);
		solve_simple(A_,b_,x_);

		x.update();

		return x;
	}

	/*! \brief Return the KSP solver
	 *
	 * In case you want to do fine tuning of the KSP solver before
	 * solve your system whith this function you can retrieve the
	 * KSP object
	 *
	 * \return the Krylov solver
	 *
	 */
	KSP getKSP()
	{
		return ksp;
	}

	/*! \brief It return the resiual error
	 *
	 * \param A Sparse matrix
	 * \param x solution
	 * \param b right-hand-side
	 *
	 * \return the solution error norms
	 *
	 */
	solError get_residual_error(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & x, const Vector<double,PETSC_BASE> & b)
	{
		return getSolNormError(A.getMat(),b.getVec(),x.getVec());
	}

	/*! \brief It return the resiual error
	 *
	 * \param x solution
	 * \param b right-hand-side
	 *
	 * \return the solution error norms
	 *
	 */
	solError get_residual_error(const Vector<double,PETSC_BASE> & x, const Vector<double,PETSC_BASE> & b)
	{
		return getSolNormError(b.getVec(),x.getVec(),ksp);
	}

	/*! \brief Here we invert the matrix and solve the system
	 *
	 * \param A sparse matrix
	 * \param b vector
	 * \param x solution and initial guess
	 *
	 * \return true if succeed
	 *
	 */
	bool solve(SparseMatrix<double,int,PETSC_BASE> & A, Vector<double,PETSC_BASE> & x, const Vector<double,PETSC_BASE> & b)
	{
		Mat & A_ = A.getMat();
		const Vec & b_ = b.getVec();
		Vec & x_ = x.getVec();

		PETSC_SAFE_CALL(KSPSetInitialGuessNonzero(ksp,PETSC_TRUE));

		pre_solve_impl(A_,b_,x_);
		solve_simple(A_,b_,x_);
		x.update();

		return true;
	}

	/*! \brief Here we invert the matrix and solve the system
	 *
	 * \param b vector
	 *
	 * \return true if succeed
	 *
	 */
	Vector<double,PETSC_BASE> solve(const Vector<double,PETSC_BASE> & b)
	{
		const Vec & b_ = b.getVec();

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt row_loc;

		PETSC_SAFE_CALL(KSPSetInitialGuessNonzero(ksp,PETSC_FALSE));
		PETSC_SAFE_CALL(VecGetSize(b_,&row));
		PETSC_SAFE_CALL(VecGetLocalSize(b_,&row_loc));

		Vector<double,PETSC_BASE> x(row,row_loc);
		Vec & x_ = x.getVec();

		solve_simple(b_,x_);
		x.update();

		return x;
	}

	/*! \brief Here we invert the matrix and solve the system
	 *
	 * In this call we are interested in solving the system
	 * with multiple right-hand-side and the same Matrix.
	 * We do not set the Matrix again and this give us the
	 * possibility to-skip the preconditioning setting that
	 * in some case like Algebraic-multi-grid can be expensive
	 *
	 * \param b vector
	 * \param x solution and initial guess
	 *
	 * \return true if succeed
	 *
	 */
	bool solve(Vector<double,PETSC_BASE> & x, const Vector<double,PETSC_BASE> & b)
	{
		const Vec & b_ = b.getVec();
		Vec & x_ = x.getVec();

		solve_simple(b_,x_);
		x.update();

		return true;
	}

	/*! \brief this function give you the possibility to set PETSC options
	 *
	 * this function call PetscOptionsSetValue
	 *
	 * \param name the name of the option
	 * \param value the value of the option
	 *
	 */
	void setPetscOption(const char * name, const char * value)
	{
		PetscOptionsSetValue(NULL,name,value);
	}

	/*! \brief Try to solve the system using all the solvers and generate a report
	 *
	 * In this mode the system will try different Solvers, Preconditioner and
	 * combination of solvers in order to find the best solver in speed, and
	 * precision. As output it will produce a performance report
	 *
	 * \param A Matrix to invert
	 * \param b right hand side
	 * \return the solution
	 *
	 */
	Vector<double,PETSC_BASE> try_solve(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & b)
	{
		Mat & A_ = A.getMat();
		const Vec & b_ = b.getVec();

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;

		PETSC_SAFE_CALL(KSPSetInitialGuessNonzero(ksp,PETSC_FALSE));
		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));

		Vector<double,PETSC_BASE> x(row,row_loc);
		Vec & x_ = x.getVec();

		pre_solve_impl(A_,b_,x_);
		try_solve_simple(A_,b_,x_);

		x.update();

		return x;
	}
};

#endif

#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_ */
