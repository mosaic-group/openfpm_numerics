/*
 * petsc_solver_ASM_AMG_to_do.hpp
 *
 *  Created on: Jun 3, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_ASM_AMG_TO_DO_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_ASM_AMG_TO_DO_HPP_


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
	void try_solve_ASM(Mat & A_, const Vec & b_, Vec & x_)
	{
		PETSC_SAFE_CALL(KSPSetTolerances(ksp,rtol,abstol,dtol,5));
		for (size_t i = 0 ; i < solvs.size() ; i++)
		{
			setSolver(solvs.get(i).c_str());

			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));

//			solve_simple(A_,b_,x_);
			solve_ASM(A_,b_,x_);
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

		////////////////

		/////////////////////


//		PCASMSetLocalSubdomains(pc);

		PCASMSetOverlap(pc,5);

		KSP       *subksp;        /* array of KSP contexts for local subblocks */
		PetscInt  nlocal,first;   /* number of local subblocks, first local subblock */
		PC        subpc;          /* PC context for subblock */
		PetscBool isasm;

		/*
		   Set runtime options
		*/
//		KSPSetFromOptions(ksp);

		/*
		  Flag an error if PCTYPE is changed from the runtime options
		*/
//		PetscObjectTypeCompare((PetscObject)pc,PCASM,&isasm);
//		if (!isasm) SETERRQ(PETSC_COMM_WORLD,1,"Cannot Change the PCTYPE when manually changing the subdomain solver settings");

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
		for (size_t i = 0 ; i < nlocal ; i++)
		{
			KSPGetPC(subksp[i],&subpc);
//			PCFactorSetFill(subpc,30);
			PCFactorSetLevels(subpc,5);
			PCSetType(subpc,PCILU);
			KSPSetType(subksp[i],KSPRICHARDSON);
			KSPSetTolerances(subksp[i],1.e-3,0.1,PETSC_DEFAULT,PETSC_DEFAULT);
		}

		// if we are on on best solve set-up a monitor function

		if (try_solve == true)
		{
			// for bench-mark we are interested in non-preconditioned norm
			PETSC_SAFE_CALL(KSPMonitorSet(ksp,monitor,&vres,NULL));

			// Disable convergence check
			PETSC_SAFE_CALL(KSPSetConvergenceTest(ksp,KSPConvergedSkip,NULL,NULL));
		}

//		KSPGMRESSetRestart(ksp,100);

		// Solve the system
		PETSC_SAFE_CALL(KSPSolve(ksp,b_,x_));

		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp,&reason);

		std::cout << "Reason: " << reason << std::endl;

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

	/////////////////////

	IS             *is,*is_local;
	PetscInt Nsub;

	PCASMCreateSubdomains2D(128,128,4,4,1,1,&Nsub,&is,&is_local);


	if (create_vcluster().getProcessUnitID() == 1)
		ISView(is_local[1],PETSC_VIEWER_STDOUT_SELF);

	////////////////////

#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_PETSC_SOLVER_ASM_AMG_TO_DO_HPP_ */
