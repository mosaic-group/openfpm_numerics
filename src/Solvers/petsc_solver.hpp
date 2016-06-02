/*
 * petsc_solver.hpp
 *
 *  Created on: Apr 26, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_

#include "Vector/Vector.hpp"
#include "Eigen/UmfPackSupport"
#include <Eigen/SparseLU>
#include <petscksp.h>
#include <petsctime.h>

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

template<>
class petsc_solver<double>
{
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

	// KSP Relative tolerance
	PetscReal rtol;

	// KSP Absolute tolerance
	PetscReal abstol;

	// KSP dtol tolerance to determine that the method diverge
	PetscReal dtol;

	// KSP Maximum number of iterations
	PetscInt maxits;

	// Parallel solver
	KSP ksp;

	// Solver used
	char s_type[32];

	// Tollerance set by the solver
	PetscScalar tol;

	// The full set of solvers
	openfpm::vector<std::string> solvs;

	// The full set of preconditionses for simple parallel solvers

	// PCType pre-conditioner type
	PCType pc;

	bool try_solve = false;

	// Residual behavior per iteration
	openfpm::vector<PetscReal> vres;

	/*! \brief procedure to collect the preconditioned residue at each iteration
	 *
	 *
	 */
	static PetscErrorCode monitor(KSP ksp,PetscInt it,PetscReal res,void* data)
	{
		openfpm::vector<PetscReal> * vres = static_cast<openfpm::vector<PetscReal> *>(data);

		vres->add(res);

		return 0;
	}

	/*! \brief try to solve the system with block jacobi pre-conditioner
	 *
	 *
	 *
	 */
	void try_solve_complex_bj(Mat & A_, const Vec & b_, Vec & x_)
	{
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol,abstol,dtol,5));
		PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));

		solve_complex(A_,b_,x_);
	}

	/*! \brief Try to solve the system using all the Krylov solver with simple preconditioners
	 *
	 *
	 *
	 */
	void try_solve_simple(Mat & A_, const Vec & b_, Vec & x_)
	{
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol,abstol,dtol,3000));
		for (size_t i = 0 ; i < solvs.size() ; i++)
		{
			setSolver(solvs.get(i).c_str());

			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));

			solve_simple(A_,b_,x_);
		}
	}

	/*! \brief solve complex use Krylov solver in combination with Block-Jacobi + Nested
	 *
	 *
	 * Solve the system using block-jacobi preconditioner + nested solver for each block
	 *
	 */
	void solve_complex(Mat & A_, const Vec & b_, Vec & x_)
	{
		// We get the size of the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;
		PetscInt * blks;
		PetscInt nlocal;
		PetscInt first;
		KSP *subksp;
		PC subpc;

		PETSC_SAFE_CALL(MatGetSize(A_,&row,&col));
		PETSC_SAFE_CALL(MatGetLocalSize(A_,&row_loc,&col_loc));

		// Set the preconditioner to Block-jacobi
		PC pc;
		KSPGetPC(ksp,&pc);
		PCSetType(pc,PCBJACOBI);
		PETSC_SAFE_CALL(KSPSetType(ksp,KSPGMRES));

		PetscInt m = 1;

		PetscMalloc1(m,&blks);
		for (size_t i = 0 ; i < m ; i++) blks[i] = row_loc;

		PCBJacobiSetLocalBlocks(pc,m,blks);
		PetscFree(blks);

		KSPSetUp(ksp);
		PCSetUp(pc);

		PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);

		for (size_t i = 0; i < nlocal; i++)
		{
			KSPGetPC(subksp[i],&subpc);

			PCSetType(subpc,PCLU);

//			PCSetType(subpc,PCJACOBI);
			KSPSetType(subksp[i],KSPPREONLY);
//			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedDefault,NULL,NULL));
//			KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

/*			if (!rank)
			{
				if (i%2)
				{
					PCSetType(subpc,PCILU);
				}
				else
				{
					PCSetType(subpc,PCNONE);
					KSPSetType(subksp[i],KSPBCGS);
					KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
				}
			}
			else
			{
					PCSetType(subpc,PCJACOBI);
					KSPSetType(subksp[i],KSPGMRES);
					KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
			}*/
		}


		KSPSolve(ksp,b_,x_);

		auto & v_cl = create_vcluster();
		if (try_solve == true)
		{
			// calculate error statistic about the solution
			solError err = statSolutionError(A_,b_,x_);

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "Method: " << s_type << " " << " pre-conditoner: " << PCJACOBI << "  iterations: " << err.its << std::endl;
				std::cout << "Norm of error: " << err.err_norm << "   Norm infinity: " << err.err_inf << std::endl;
			}
		}
	}

	void solve_ASM(Mat & A_, const Vec & b_, Vec & x_)
	{


#if 0
		PETSC_SAFE_CALL(KSPSetType(ksp,s_type));

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
		PETSC_SAFE_CALL(PCSetType(pc,PCASM));

		PCASMSetOverlap(pc,1);

		KSP       *subksp;        /* array of KSP contexts for local subblocks */
		PetscInt  nlocal,first;   /* number of local subblocks, first local subblock */
		PC        subpc;          /* PC context for subblock */
		PetscBool isasm;

		/*
		   Set runtime options
		*/
		KSPSetFromOptions(ksp);

		/*
		  Flag an error if PCTYPE is changed from the runtime options
		*/
		PetscObjectTypeCompare((PetscObject)pc,PCASM,&isasm);
		if (!isasm) SETERRQ(PETSC_COMM_WORLD,1,"Cannot Change the PCTYPE when manually changing the subdomain solver settings");

		/*
		  Call KSPSetUp() to set the block Jacobi data structures (including
		  creation of an internal KSP context for each block).

		  Note: KSPSetUp() MUST be called before PCASMGetSubKSP().
		*/
		KSPSetUp(ksp);

		/*
		   Extract the array of KSP contexts for the local blocks
		*/
		PCASMGetSubKSP(pc,&nlocal,&first,&subksp);

		/*
		   Loop over the local blocks, setting various KSP options
		   for each block.
		*/
		for (i=0; i<nlocal; i++)
		{
			KSPGetPC(subksp[i],&subpc);
			PCSetType(subpc,PCILU);
			KSPSetType(subksp[i],KSPGMRES);
			KSPSetTolerances(subksp[i],1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		}

		// if we are on on best solve set-up a monitor function

		if (try_solve == true)
		{
			// for bench-mark we are interested in non-preconditioned norm
			PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor,&vres,NULL));

			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
		}

		KSPGMRESSetRestart(ksp,200);

		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));

		auto & v_cl = create_vcluster();
//		if (try_solve == true)
//		{
			// calculate error statistic about the solution
			solError err = statSolutionError(A_,b_,x_);

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "Method: " << s_type << " " << " pre-conditoner: " << PCJACOBI << "  iterations: " << err.its << std::endl;
				std::cout << "Norm of error: " << err.err_norm << "   Norm infinity: " << err.err_inf << std::endl;
			}
//		}

#endif
	}

	void solve_AMG(Mat & A_, const Vec & b_, Vec & x_)
	{
		//		PETSC_SAFE_CALL(KSPSetType(ksp,s_type));
				PETSC_SAFE_CALL(KSPSetType(ksp,KSPRICHARDSON));

		//		-pc_hypre_boomeramg_tol <tol>

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

				// PETSC_SAFE_CALL(PCSetType(pc,PCJACOBI));

				PETSC_SAFE_CALL(PCSetType(pc,PCHYPRE));
		//		PCGAMGSetNSmooths(pc,0);
			    PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
			    PCFactorSetShiftAmount(pc, PETSC_DECIDE);
				PCHYPRESetType(pc, "boomeramg");
				MatSetBlockSize(A_,4);
				PetscOptionsSetValue("-pc_hypre_boomeramg_print_statistics","2");
				PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter","1000");
				PetscOptionsSetValue("-pc_hypre_boomeramg_nodal_coarsen","true");
				PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all","SOR/Jacobi");
				PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","Falgout");
				PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type","W");
				PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels","20");
				PetscOptionsSetValue("-pc_hypre_boomeramg_vec_interp_variant","0");
				KSPSetFromOptions(ksp);

				// if we are on on best solve set-up a monitor function

				if (try_solve == true)
				{
					// for bench-mark we are interested in non-preconditioned norm
					PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor,&vres,NULL));

					// Disable convergence check
					PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
				}

				// Solve the system
				PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));

				auto & v_cl = create_vcluster();
		//		if (try_solve == true)
		//		{
					// calculate error statistic about the solution
					solError err = statSolutionError(A_,b_,x_);

					if (v_cl.getProcessUnitID() == 0)
					{
						std::cout << "Method: " << s_type << " " << " pre-conditoner: " << PCJACOBI << "  iterations: " << err.its << std::endl;
						std::cout << "Norm of error: " << err.err_norm << "   Norm infinity: " << err.err_inf << std::endl;
					}
		//		}
	}

	/*! \brief solve simple use a Krylov solver + Simple selected Parallel Pre-conditioner
	 *
	 */
	void solve_Krylov_simple(Mat & A_, const Vec & b_, Vec & x_)
	{
		PETSC_SAFE_CALL(KSPSetType(ksp,s_type));

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

		// PETSC_SAFE_CALL(PCSetType(pc,PCJACOBI));

		PETSC_SAFE_CALL(PCSetType(pc,PCHYPRE));
//		PCGAMGSetNSmooths(pc,0);
	    PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
	    PCFactorSetShiftAmount(pc, PETSC_DECIDE);
		PCHYPRESetType(pc, "boomeramg");
		MatSetBlockSize(A_,4);
		PetscOptionsSetValue("-pc_hypre_boomeramg_print_statistics","2");
		PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter","1000");
		PetscOptionsSetValue("-pc_hypre_boomeramg_nodal_coarsen","true");
		PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all","SOR/Jacobi");
		PetscOptionsSetValue("-pc_hypre_boomeramg_coarsen_type","Falgout");
		PetscOptionsSetValue("-pc_hypre_boomeramg_cycle_type","W");
		PetscOptionsSetValue("-pc_hypre_boomeramg_max_levels","10");
		KSPSetFromOptions(ksp);
		// if we are on on best solve set-up a monitor function

		if (try_solve == true)
		{
			// for bench-mark we are interested in non-preconditioned norm
			PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor,&vres,NULL));

			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
		}

		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));

		auto & v_cl = create_vcluster();
//		if (try_solve == true)
//		{
			// calculate error statistic about the solution
			solError err = statSolutionError(A_,b_,x_);

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "Method: " << s_type << " " << " pre-conditoner: " << PCJACOBI << "  iterations: " << err.its << std::endl;
				std::cout << "Norm of error: " << err.err_norm << "   Norm infinity: " << err.err_inf << std::endl;
			}
//		}
	}

	/*! \brief solve simple use a Krylov solver + Simple selected Parallel Pre-conditioner
	 *
	 */
	void solve_simple(Mat & A_, const Vec & b_, Vec & x_)
	{
		PETSC_SAFE_CALL(KSPSetType(ksp,s_type));

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
			// for bench-mark we are interested in non-preconditioned norm
			PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor,&vres,NULL));

			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
		}

		KSPGMRESSetRestart(ksp,300);

		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));

		auto & v_cl = create_vcluster();
//		if (try_solve == true)
//		{
			// calculate error statistic about the solution
			solError err = statSolutionError(A_,b_,x_);

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "Method: " << s_type << " " << " pre-conditoner: " << PCJACOBI << "  iterations: " << err.its << std::endl;
				std::cout << "Norm of error: " << err.err_norm << "   Norm infinity: " << err.err_inf << std::endl;
			}
//		}
	}

	/*! \brief Calculate statistic on the error solution
	 *
	 * \brief A Matrix of the system
	 * \brief b Right hand side of the matrix
	 * \brief x Solution
	 *
	 */
	solError statSolutionError(Mat & A_, const Vec & b_, Vec & x_)
	{
		PetscScalar neg_one = -1.0;

		// We set the size of x according to the Matrix A
		PetscInt row;
		PetscInt col;
		PetscInt row_loc;
		PetscInt col_loc;
		PetscInt its;

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
		PETSC_SAFE_CALL(KSPGetIterationNumber(ksp,&its));

		solError err;
		err.err_norm = norm;
		err.err_inf = norm_inf;
		err.its = its;

		return err;
	}


public:

	petsc_solver()
	{
		strcpy(s_type,KSPBCGSL);
		PETSC_SAFE_CALL(KSPCreate(PETSC_COMM_WORLD,&ksp));
		PETSC_SAFE_CALL(KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits));

		// Add the solvers

//		solvs.add(std::string(KSPBCGS));
//		solvs.add(std::string(KSPIBCGS));
//		solvs.add(std::string(KSPFBCGS));
//		solvs.add(std::string(KSPFBCGSR));
//		solvs.add(std::string(KSPBCGSL));
		solvs.add(std::string(KSPGMRES));
//		solvs.add(std::string(KSPFGMRES));
//		solvs.add(std::string(KSPLGMRES));
//		solvs.add(std::string(KSPPGMRES));
//		solvs.add(std::string(KSPGCR));
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
		strcpy(s_type,type);
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
		rtol = rtol_;
		KSPSetTolerances(ksp,rtol,abstol,dtol,maxits);
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
		abstol = abstol_;
		KSPSetTolerances(ksp,0.0,abstol,30000.0,1000);
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
		dtol = dtol_;
		KSPSetTolerances(ksp,rtol,abstol,dtol,maxits);
	}

	/*! For the BiCGStab(L) it define the number of search directions
	 *
	 * BiCG Methods can fail for base break-down (or near break-down). BiCGStab(L) try
	 * to avoid this problem. Such method has a parameter L (2 by default) Bigger is L
	 * more the system will try to avoid the breakdown
	 *
	 * \size_t l Increasing L should reduce the probability of failure of the solver because of break-down of the base
	 *
	 */
	void searchDirections(PetscInt l)
	{
		KSPBCGSLSetEll(ksp, l);
	}

	/*! \brief Here we invert the matrix and solve the system
	 *
	 *  \warning umfpack is not a parallel solver, this function work only with one processor
	 *
	 *  \note if you want to use umfpack in a NON parallel, but on a distributed data, use solve with triplet
	 *
	 *	\tparam impl Implementation of the SparseMatrix
	 *
	 */
	Vector<double,PETSC_BASE> solve(SparseMatrix<double,int,PETSC_BASE> & A, const Vector<double,PETSC_BASE> & b)
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

		Vector<double,PETSC_BASE> x(row,row_loc);
		Vec & x_ = x.getVec();

		if (try_solve == true)
		{
			try_solve_simple(A_,b_,x_);
//			try_solve_complex_bj(A_,b_,x_);
		}
		else
		{
			solve_simple(A_,b_,x_);
		}

		x.update();
		return x;
	}
};


#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_HPP_ */
