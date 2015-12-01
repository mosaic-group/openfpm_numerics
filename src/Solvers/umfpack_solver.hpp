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

template<typename T>
class umfpack_solver
{
public:

	template<typename impl> static Vector<T> solve(const SparseMatrix<T,impl> & A, const Vector<T> & b)
	{
		std::cerr << "Error Umfpack only suppor double precision" << "/n";
	}
};

template<>
class umfpack_solver<double>
{
public:

	template<typename impl> static Vector<double> solve(const SparseMatrix<double,int,impl> & A, const Vector<double> & b)
	{
		Vector<double> x;

		Eigen::UmfPackLU<Eigen::SparseMatrix<double,0,int> > solver;

		solver.compute(A.getMat());

		if(solver.info()!=Eigen::Success)
		{
			// decomposition failed
			std::cout << __FILE__ << ":" << __LINE__ << " solver failed" << "\n";
			return x;
		}
		x.getVec() = solver.solve(b.getVec());

		return x;
	}
};



#endif /* OPENFPM_NUMERICS_SRC_SOLVERS_UMFPACK_SOLVER_HPP_ */
