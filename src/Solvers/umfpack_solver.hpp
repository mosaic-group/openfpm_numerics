/*
 * Umfpack_solver.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_SOLVERS_UMFPACK_SOLVER_HPP_
#define OPENFPM_NUMERICS_SRC_SOLVERS_UMFPACK_SOLVER_HPP_

#include "Vector/Vector.hpp"
#include "Eigen/UmfPackSupport"
#include <Eigen/SparseLU>

#define UMFPACK_NONE 0

template<typename T>
class umfpack_solver
{
public:

	template<typename impl> static Vector<T> solve(const SparseMatrix<T,impl> & A, const Vector<T> & b)
	{
		std::cerr << "Error Umfpack only suppor double precision" << "/n";
	}
};

#define SOLVER_NOOPTION 0
#define SOLVER_PRINT_RESIDUAL_NORM_INFINITY 1
#define SOLVER_PRINT_DETERMINANT 2

template<>
class umfpack_solver<double>
{

public:

	/*! \brief Here we invert the matrix and solve the system
	 *
	 *  \warning umfpack is not a parallel solver, this function work only with one processor
	 *
	 *  \note if you want to use umfpack in a NON parallel, but on a distributed data, use solve with triplet
	 *
	 *	\tparam impl Implementation of the SparseMatrix
	 *
	 */
	template<typename impl> static Vector<double> solve(SparseMatrix<double,int,impl> & A, const Vector<double> & b, size_t opt = UMFPACK_NONE)
	{
		Vcluster & vcl = create_vcluster();

		Vector<double> x;

		// only master processor solve
		Eigen::UmfPackLU<Eigen::SparseMatrix<double,0,int> > solver;

		// Collect the matrix on master
		auto mat_ei = A.getMat();

		Eigen::Matrix<double, Eigen::Dynamic, 1> x_ei;

		// Collect the vector on master
		auto b_ei = b.getVec();

		// Copy b into x, this also copy the information on how to scatter back the information on x
		x = b;

		if (vcl.getProcessUnitID() == 0)
		{
			solver.compute(mat_ei);

			if(solver.info()!=Eigen::Success)
			{
				// decomposition failed
				std::cout << __FILE__ << ":" << __LINE__ << " solver failed" << "\n";
				return x;
			}

			x_ei = solver.solve(b_ei);

			if (opt & SOLVER_PRINT_RESIDUAL_NORM_INFINITY)
			{
				Eigen::Matrix<double, Eigen::Dynamic, 1> res;
				res = mat_ei * x_ei - b_ei;

				std::cout << "Infinity norm: " << res.lpNorm<Eigen::Infinity>() << "\n";
			}

			if (opt & SOLVER_PRINT_DETERMINANT)
			{
				std::cout << " Determinant: " << solver.determinant() << "\n";
			}

			x = x_ei;
		}

		// Vector is only on master, scatter back the information
		x.scatter();

		return x;
	}
};



#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_UMFPACK_SOLVER_HPP_ */
