/*
 * petsc_solver.hpp
 *
 *  Created on: Apr 26, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_

#ifdef HAVE_PETSC

#include "Vector/Vector.hpp"
#include "Eigen/UmfPackSupport"
#include <Eigen/SparseLU>
#include <petscksp.h>
#include <petsctime.h>
#include "Plot/GoogleChart.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"

#define UMFPACK_NONE 0


template<typename T>
class petsc_solver
{
public:

	template<typename impl> static Vector<T> solve(const SparseMatrix<T,impl> & A, const Vector<T> & b)
	{
		std::cerr << "Error Petsc only suppor double precision" << "/n";
	}
};

#define SOLVER_NOOPTION 0
#define SOLVER_PRINT_RESIDUAL_NORM_INFINITY 1
#define SOLVER_PRINT_DETERMINANT 2

// It contain statistic of the error of the calculated solution
struct solError
{
	// infinity norm of the error
	PetscReal err_inf;

	// L1 norm of the error
	PetscReal err_norm;

	// Number of iterations
	PetscInt its;
};

size_t cpu_rank;

/*! \brief This class is able to do Matrix inversion in parallel with PETSC solvers
 *
 * #Example of how to solve a lid driven cavity 3D with a PETSC solver in paralle
 * \snippet
 *
 */
template<>
class petsc_solver<double>
{
	//It contain statistic of the error at each iteration
	struct itError
	{
		PetscInt it;
		PetscReal err_inf;
		PetscReal err_norm;

		itError(const solError & e)
		{
			it = e.its;
			err_inf = e.err_inf;
			err_norm = e.err_norm;
		}

		// Default constructor
		itError()	{}
	};

	// It contain the benchmark
	struct solv_bench_info
	{
		// Method name
		std::string method;

		//Method name (short form)
		std::string smethod;

		// time to converge in milliseconds
		double time;

		// Solution error
		solError err;

		// Convergence per iteration
		openfpm::vector<itError> res;
	};

	// KSP Maximum number of iterations
	PetscInt maxits;

	// Main parallel solver
	KSP ksp;

	// Temporal variable used for calculation of static members
	size_t tmp;

	// The full set of solvers
	openfpm::vector<std::string> solvs;

	bool try_solve = false;

	// It contain the solver benchmark results
	openfpm::vector<solv_bench_info> bench;

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
			return slv.res.get(pt).err_inf;

		return slv.res.last().err_inf;
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

	/*! \brief procedure to collect the absolute residue at each iteration
	 *
	 * \param ksp Solver
	 * \param it Iteration number
	 * \param res resudual
	 * \param custom pointer to data
	 *
	 */
	static PetscErrorCode monitor(KSP ksp,PetscInt it,PetscReal res,void* data)
	{
		petsc_solver<double> * pts = static_cast<petsc_solver<double> *>(data);

		Mat A;
		Mat P;
		Vec x;
		Vec b;

		PETSC_SAFE_CALL(KSPGetOperators(ksp,&A,&P));
		PETSC_SAFE_CALL(KSPGetSolution(ksp,&x));
		PETSC_SAFE_CALL(KSPGetRhs(ksp,&b));

		solError err = statSolutionError(A,b,x,ksp);
		itError erri(err);
		erri.it = it;

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


		pts->progress(it);

		return 0;
	}

	/*! \brief procedure print the progress of the solver in benchmark mode
	 *
	 * \param ksp Solver
	 * \param it Iteration number
	 * \param res resudual
	 * \param custom pointer to data
	 *
	 */
	static PetscErrorCode monitor_progress(KSP ksp,PetscInt it,PetscReal res,void* data)
	{
		petsc_solver<double> * pts = (petsc_solver *)data;

		pts->progress(it);

		return 0;
	}

	/*! \brief This function print an "*" showing the progress of the solvers
	 *
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
	 * \param x solution
	 *
	 */
	void try_solve_simple(Mat & A_, const Vec & b_, Vec & x_)
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
	void bench_solve_simple(Mat & A_, const Vec & b_, Vec & x_, solv_bench_info & bench)
	{
		// timer for benchmark
		timer t;
		t.start();
		// Enable progress monitor
		tmp = 0;
		PETSC_SAFE_CALL(KSPMonitorCancel(ksp));
		PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor_progress,this,NULL));
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

		// Solve getting the residual per iteration
		tmp = 0;
		PETSC_SAFE_CALL(KSPMonitorCancel(ksp));
		PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor,this,NULL));
		print_progress_bar();
		solve_simple(A_,b_,x_);

		// New line
		if (create_vcluster().getProcessUnitID() == 0)
			std::cout << std::endl;
	}

	/*! \brief solve simple use a Krylov solver + Simple selected Parallel Pre-conditioner
	 *
	 */
	void solve_simple(Mat & A_, const Vec & b_, Vec & x_)
	{
//		PETSC_SAFE_CALL(KSPSetType(ksp,s_type));

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;

		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));

		PC pc;

		// We set the Matrix operators
		PETSC_SAFE_CALL(KSPSetOperators(ksp,A_,A_));

		// We set the pre-conditioner
		PETSC_SAFE_CALL(KSPGetPC(ksp,&pc));
		PETSC_SAFE_CALL(PCSetType(pc,PCJACOBI));

		// if we are on on best solve set-up a monitor function

		if (try_solve == true)
		{
			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
		}

		PETSC_SAFE_CALL(KSPSetFromOptions(ksp));

		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));
	}

	/*! \brief Calculate statistic on the error solution
	 *
	 * \brief A Matrix of the system
	 * \brief b Right hand side of the matrix
	 * \brief x Solution
	 *
	 */
	static solError statSolutionError(Mat & A_, const Vec & b_, Vec & x_, KSP ksp)
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

		// Disable convergence check
		PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
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
	 * \param A_ the matrix that identity the linear system
	 * \param x_ the solution
	 * \param b_ the right-hand-side
	 *
	 */
	static solError getSolNormError(const Mat & A_, const Vec & x_, const Vec & b_)
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

		return err;
	}

public:

	typedef Vector<double,PETSC_BASE> return_type;

	~petsc_solver()
	{
		PETSC_SAFE_CALL(KSPDestroy(&ksp));
	}

	petsc_solver()
	:maxits(300)
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

	/*! \brief Add a test solver
	 *
	 * The try solve function use the most robust solvers in PETSC, if you want to add
	 * additionally solver like KSPIBCGS,KSPFBCGSR,KSPPGMRES, use addTestSolver(std::string(KSPIBCGS))
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

	/*! \brief Set the solver to test all the solvers and generate a report
	 *
	 * In this mode the system will try different Solvers, Preconditioner and
	 * combination of solvers in order to find the best solver in speed, and
	 * precision. As output it will produce a performance report
	 *
	 *
	 */
	void best_solve()
	{
		try_solve = true;
	}

	/*! \brief Set the Petsc solver
	 *
	 * \see KSPType in PETSC manual for a list of all PETSC solvers
	 *
	 */
	void setSolver(KSPType type)
	{
		PetscOptionsSetValue("-ksp_type",type);
	}

	/*! \brief Set the relative tolerance as stop criteria
	 *
	 * \see PETSC manual KSPSetTolerances for an explanation
	 *
	 * \param rtol Relative tolerance
	 *
	 */
	void setRelTol(PetscReal rtol_)
	{
		PetscOptionsSetValue("-ksp_rtol",std::to_string(rtol_).c_str());
	}

	/*! \brief Set the absolute tolerance as stop criteria
	 *
	 * \see PETSC manual KSPSetTolerances for an explanation
	 *
	 * \param abstol Absolute tolerance
	 *
	 */
	void setAbsTol(PetscReal abstol_)
	{
		PetscOptionsSetValue("-ksp_atol",std::to_string(abstol_).c_str());
	}

	/*! \brief Set the divergence tolerance
	 *
	 * \see PETSC manual KSPSetTolerances for an explanation
	 *
	 * \param divtol
	 *
	 */
	void setDivTol(PetscReal dtol_)
	{
		PetscOptionsSetValue("-ksp_dtol",std::to_string(dtol_).c_str());
	}

	/*! \brief Set the maximum number of iteration for Krylov solvers
	 *
	 *
	 */
	void setMaxIter(PetscInt n)
	{
		PetscOptionsSetValue("-ksp_max_it",std::to_string(n).c_str());
		maxits = n;
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
		PetscOptionsSetValue("-ksp_bcgsl_ell",std::to_string(l).c_str());
	}

	/*! \brief For GMRES based method,  the number of Krylov directions to orthogonalize against
	 *
	 * \param n number of directions
	 *
	 */
	void setRestart(PetscInt n)
	{
		PetscOptionsSetValue("-ksp_gmres_restart",std::to_string(n).c_str());
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

		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));
		PETSC_SAFE_CALL(KSPSetInitialGuessNonzero(ksp,PETSC_FALSE));

		Vector<double,PETSC_BASE> x(row,row_loc);
		Vec & x_ = x.getVec();

		if (try_solve == true)
			try_solve_simple(A_,b_,x_);
		else
			solve_simple(A_,b_,x_);

		x.update();
		return x;
	}

	/*! \brief It return the resiual error
	 *
	 *
	 */
	solError get_residual_error(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & x, const Vector<double,PETSC_BASE> & b)
	{
		return getSolNormError(A.getMat(),x.getVec(),b.getVec());
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
	 * \param x solution and initial guess
	 * \param initial_guess true if x has the initial guess
	 *
	 * \return true if succeed
	 *
	 */
	bool solve(SparseMatrix<double,int,PETSC_BASE> & A, Vector<double,PETSC_BASE> & x, const Vector<double,PETSC_BASE> & b)
	{
		Mat & A_ = A.getMat();
		const Vec & b_ = b.getVec();

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;

		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));
		PETSC_SAFE_CALL(KSPSetInitialGuessNonzero(ksp,PETSC_TRUE));


		Vec & x_ = x.getVec();

		if (try_solve == true)
			try_solve_simple(A_,b_,x_);
		else
			solve_simple(A_,b_,x_);

		x.update();
		return true;
	}
};

#endif

#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_ */
