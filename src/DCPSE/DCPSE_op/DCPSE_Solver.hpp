//
// Created by Abhinav Singh on 20.01.20.
//

#ifndef OPENFPM_PDATA_DCPSE_SOLVER_HPP
#define OPENFPM_PDATA_DCPSE_SOLVER_HPP
#ifdef HAVE_EIGEN


#include "DCPSE_op.hpp"
#include "Matrix/SparseMatrix.hpp"
#include "Vector/Vector.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "Vector/Vector_util.hpp"
#include "Vector/vector_dist.hpp"
#include "Solvers/umfpack_solver.hpp"
#include "Solvers/petsc_solver.hpp"
#include "util/eq_solve_common.hpp"

/*enum eq_struct
{
	VECTOR,
	SCALAR
};*/

//template<unsigned int prp_id> using prop_id = boost::mpl::int_<prp_id>;
/*! \brief Create a Matrix System for Ax=b
 *
 * This Class is for creating a placeholder for the matrix system
 *
 * Ax = b
 *
 *
 * \param Sys_eqs Equation Structure which has information about the system. Refer to EqnStruct.cpp for examples
 * \param parts Particle set
 *
 */
template<typename Sys_eqs, typename particles_type>
class DCPSE_scheme {

    //! type of the sparse matrix
    typename Sys_eqs::SparseMatrix_type A;

    //! Vector b
    typename Sys_eqs::Vector_type b;

    //! Vector x_ig (initial guess)
    typename Sys_eqs::Vector_type x_ig;

    //! Sparse matrix triplet type
    typedef typename Sys_eqs::SparseMatrix_type::triplet_type triplet;

    //! Distributed grid map
    typedef vector_dist<Sys_eqs::dims, typename Sys_eqs::stype, aggregate<size_t>> p_map_type;

    //! mapping grid
    p_map_type p_map;

    //! Grid points that has each processor
    openfpm::vector<size_t> pnt;

    //! Particles used to impose the system
    particles_type &parts;

    //! colums shift map
    //int col_sm[Sys_eqs::nvar];

    //! Each point in the grid has a global id, to decompose correctly the Matrix each processor contain a
    //! contiguos range of global id, example processor 0 can have from 0 to 234 and processor 1 from 235 to 512
    //! no processors can have holes in the sequence, this number indicate where the sequence start for this
    //! processor
    size_t s_pnt;

    //! row of the matrix
    size_t row;

    //! row on b
    size_t row_b;

    //! row on x_ig
    size_t row_x_ig;

    //! Total number of points
    size_t tot;

    //! solver options
    options_solver opt;

    size_t offset;


    /*! \brief Construct the gmap structure
 *
 */template<typename options>
    void construct_pmap(options opt = options_solver::STANDARD) {
        Vcluster<> &v_cl = create_vcluster();

        // Calculate the size of the local domain
        size_t sz = p_map.size_local();

        // Get the total size of the local grids on each processors
        v_cl.allGather(sz, pnt);
        v_cl.execute();
        s_pnt = 0;

        // calculate the starting point for this processor
        for (size_t i = 0; i < v_cl.getProcessUnitID(); i++)
            s_pnt += pnt.get(i);

        tot = sz;
        v_cl.sum(tot);
        v_cl.execute();

        // resize b if needed
        if (opt == options_solver::STANDARD) {
            b.resize(Sys_eqs::nvar * tot, Sys_eqs::nvar * sz);
            x_ig.resize(Sys_eqs::nvar * tot, Sys_eqs::nvar * sz);

        } else if (opt == options_solver::LAGRANGE_MULTIPLIER) {
            if (v_cl.rank() == v_cl.size() - 1) {
                b.resize(Sys_eqs::nvar * (tot + 1), Sys_eqs::nvar * (sz + 1));
                x_ig.resize(Sys_eqs::nvar * (tot + 1), Sys_eqs::nvar * (sz + 1));
            } else {
                b.resize(Sys_eqs::nvar * (tot + 1), Sys_eqs::nvar * sz);
                x_ig.resize(Sys_eqs::nvar * (tot + 1), Sys_eqs::nvar * sz);
            }
        }
            //Use Custom number of constraints using opt as an integer
        else {
            if (v_cl.rank() == v_cl.size() - 1) {
                b.resize(Sys_eqs::nvar * tot - offset, Sys_eqs::nvar * sz - offset);
                x_ig.resize(Sys_eqs::nvar * tot - offset, Sys_eqs::nvar * sz - offset);
            } else {
                b.resize(Sys_eqs::nvar * tot - offset, Sys_eqs::nvar * sz);
                x_ig.resize(Sys_eqs::nvar * tot - offset, Sys_eqs::nvar * sz);
            }
        }

        // Calculate the starting point

        // Counter
        size_t cnt = 0;

        // Create the re-mapping grid
        auto it = p_map.getDomainIterator();

        while (it.isNext()) {
            auto key = it.get();

            for (int i = 0; i < particles_type::dims; i++) {
                p_map.getPos(key)[i] = parts.getPos(key)[i];
            }

            p_map.template getProp<0>(key) = cnt + s_pnt;

            ++cnt;
            ++it;
        }

        // sync the ghost
        p_map.template ghost_get<0>();
    }

    //! Encapsulation of the b term as constant
    struct constant_b {
        //! scalar
        typename Sys_eqs::stype scal;

        /*! \brief Constrictor from a scalar
         *
         * \param scal scalar
         *
         */
        constant_b(typename Sys_eqs::stype scal) {
            this->scal = scal;
        }

        /*! \brief Get the b term on a grid point
         *
         * \note It does not matter the grid point it is a scalar
         *
         * \param  key grid position (unused because it is a constant)
         *
         * \return the scalar
         *
         */
        typename Sys_eqs::stype get(size_t key) {
            return scal;
        }
    };

    //! Encapsulation of the b term as constant
    template<unsigned int prp_id>
    struct variable_b {
        //! scalar
        typename Sys_eqs::stype scal;

        particles_type &parts;

        /*! \brief Constrictor from a scalar
         *
         * \param scal scalar
         *
         */
        variable_b(particles_type &parts)
                : parts(parts) {}

        /*! \brief Get the b term on a grid point
         *
         * \note It does not matter the grid point it is a scalar
         *
         * \param  key grid position (unused because it is a constant)
         *
         * \return the scalar
         *
         */
        inline typename Sys_eqs::stype get(size_t key) {
            return parts.template getProp<prp_id>(key);
        }
    };


    /*! \brief Check if the Matrix is consistent
 *
 */
    void consistency(options_solver opt)
    {
        openfpm::vector<triplet> &trpl = A.getMatrixTriplets();
        Vcluster<> &v_cl = create_vcluster();

        // A and B must have the same rows
        if (row != row_b) {
            std::cerr << "Error " << __FILE__ << ":" << __LINE__
                      << " the term B and the Matrix A for Ax=B must contain the same number of rows " << row
                      << "!=" << row_b << "\n";
            return;
        }

        if (row_b != p_map.size_local() * Sys_eqs::nvar) {
            std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " your system is underdetermined you set "
                      << row_b << " conditions " << " but i am expecting " << p_map.size_local() * Sys_eqs::nvar
                      << std::endl;
            return;
        }
        // Indicate all the non zero rows
        openfpm::vector<unsigned char> nz_rows;


        if (v_cl.rank() == v_cl.size()-1 && opt == options_solver::LAGRANGE_MULTIPLIER) {
            nz_rows.resize(row_b+Sys_eqs::nvar);
            }
        else{
            nz_rows.resize(row_b);
        };

        for (size_t i = 0; i < trpl.size(); i++) {
            if (trpl.get(i).row() - s_pnt * Sys_eqs::nvar >= nz_rows.size()) {
                std::cerr << "Error " << __FILE__ << ":" << __LINE__
                          << " It seems that you are setting colums that does not exist \n";
            }
            if (trpl.get(i).value() != 0) { nz_rows.get(trpl.get(i).row() - s_pnt * Sys_eqs::nvar) = true; }
        }

        // Indicate all the non zero colums
        // This check can be done only on single processor

        if (v_cl.getProcessingUnits() == 1) {
            openfpm::vector<unsigned> nz_cols;
            if (v_cl.rank() == v_cl.size()-1 && opt == options_solver::LAGRANGE_MULTIPLIER) {
                nz_cols.resize(row_b+Sys_eqs::nvar);
            }
            else{
                nz_cols.resize(row_b);
            };

            for (size_t i = 0; i < trpl.size(); i++) {
                if (trpl.get(i).value() != 0) { nz_cols.get(trpl.get(i).col()) = true; }
            }

            // all the rows must have a non zero element
            for (size_t i = 0; i < nz_rows.size(); i++) {
                if (nz_rows.get(i) == false) {
                    std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Ill posed matrix row " << i
                              << " is not filled " << " equation: " << "\n";
                }
            }

            // all the colums must have a non zero element
            for (size_t i = 0; i < nz_cols.size(); i++) {
                if (nz_cols.get(i) == false)
                    std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Ill posed matrix colum " << i
                              << " is not filled\n";
            }
        }
    }

    /*! \brief Solve an equation
    *
    *  \warning exp must be a scalar type
    *
    * \param exp where to store the result
    *
    */
    template<typename solType, typename expr_type>
    void copy_impl(solType & x, expr_type exp, unsigned int comp)
    {
        auto & parts = exp.getVector();

        auto it = parts.getDomainIterator();

        while (it.isNext()) {
            auto p = it.get();
            exp.value(p) = x(p.getKey() * Sys_eqs::nvar + comp + s_pnt * Sys_eqs::nvar);
            ++it;
        }
    }

    template<typename solType, typename exp1, typename ... othersExp>
    void copy_nested(solType &x, unsigned int &comp, exp1 exp, othersExp ... exps) {
        copy_impl(x, exp, comp);
        comp++;

        copy_nested(x, comp, exps ...);
    }


    template<typename solType, typename exp1>
    void copy_nested(solType &x, unsigned int &comp, exp1 exp) {
        copy_impl(x, exp, comp);
        comp++;
    }

public:

    /*! \brief Set the structure of the system of equation
     *
     * For example for stokes-flow where you are solving for V = velocity (Vector) and P = pressure (scalar)
     *
     * you should call this function with
     *
     * setEquationStructure({eq_struct::VECTOR,eq_struct::SCALAR})
     *
     */
/*    void setEquationStructure(std::initializer_list<eq_struct> l)
    {
    	int i = 0;
    	for (eq_struct e : l)
    	{
    		if (e == eq_struct::VECTOR)
    		{
    			for (int j = 0 ; j < Sys_eqs::dims ; j++)
    			{
    				col_sm[i+j] = i;
    			}
    			i += Sys_eqs::dims;
    		}
    		else
    		{
    			col_sm[i] = i;
    		}
    	}
    }*/


    /*! \brief Solve an equation
     *
     *  \warning exp must be a scalar type
     *
     * \param exp where to store the result
     *
     */
    template<typename ... expr_type>
    void solve(expr_type ... exps) {
        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
        typename Sys_eqs::solver_type solver;
//        umfpack_solver<double> solver;
        auto x = solver.solve(getA(opt), getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

    /*! \brief Solve an equation
     *
     *  \warning exp must be a scalar type
     *
     * \param Solver Manually created Solver instead from the Equation structure
     * \param exp where to store the result
     *
     */
    template<typename SolverType, typename ... expr_type>
    void solve_with_solver(SolverType &solver, expr_type ... exps) {
#ifdef SE_CLASS1

        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
#endif
        auto x = solver.solve(getA(opt), getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

    /*! \brief Solve an equation
 *
 *  \warning exp must be a scalar type
 *
 * \param Solver Manually created Solver instead from the Equation structure
 * \param exp where to store the result
 *
 */
    template<typename SolverType, typename ... expr_type>
    void solve_with_solver_ig(SolverType &solver,expr_type ... exps) {
#ifdef SE_CLASS1

        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
#endif
        auto x = solver.solve(getA(opt),get_x_ig(opt),getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

        /*! \brief Successive Solve an equation
 *
 *  \warning exp must be a scalar type
 *
 * \param Solver Manually created Solver instead from the Equation structure
 * \param exp where to store the result
 *
 */
    template<typename SolverType, typename ... expr_type>
    void solve_with_solver_successive(SolverType &solver,expr_type ... exps) {
#ifdef SE_CLASS1

        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
#endif
        auto x = solver.solve_successive(getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

    /*! \brief Successive Solve an equation with inital guess
 *
 *  \warning exp must be a scalar type
 *
 * \param Solver Manually created Solver instead from the Equation structure
 * \param exp where to store the result
 *
 */
    template<typename SolverType, typename ... expr_type>
    void solve_with_solver_ig_successive(SolverType &solver,expr_type ... exps) {
#ifdef SE_CLASS1

        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
#endif
        auto x = solver.solve_successive(get_x_ig(opt),getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

    /*! \brief Solve an equation with a given Nullspace
     *
     *  \warning exp must be a scalar type
     *
     * \param Solver Manually created Solver instead from the Equation structure
     * \param exp where to store the result
     *
     */
/*    template<typename NullspaceType, typename SolverType, typename ... expr_type>
    void solve_with_nullspace_solver(NullspaceType &nullspace,SolverType &solver, expr_type ... exps) {
#ifdef SE_CLASS1

        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
#endif
        auto x = solver.nullspace_solve(nullspace,getA(opt), getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }*/

    /*! \brief Solve an equation with a constant Nullspace from PETSC Hack
     *
     *  \warning exp must be a scalar type
     *
     * \param Solver Manually created Solver instead from the Equation structure
     * \param exp where to store the result
     *
     */
    template<typename SolverType, typename ... expr_type>
    void solve_with_default_nullspace_solver(SolverType &solver, expr_type ... exps) {
#ifdef SE_CLASS1

        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };
#endif
        auto x = solver.with_nullspace_solve(getA(opt), getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

    /*! \brief Solve an equation
     *
     *  \warning exp must be a scalar type
     *
     * \param exp where to store the result
     *
     */
    template<typename SolverType, typename ... expr_type>
    void try_solve_with_solver(SolverType &solver, expr_type ... exps) {
        if (sizeof...(exps) != Sys_eqs::nvar) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the number of properties you gave does not match the solution in\
    													dimensionality, I am expecting " << Sys_eqs::nvar <<
                      " properties " << std::endl;
        };

        auto x = solver.try_solve(getA(opt), getB(opt));

        unsigned int comp = 0;
        copy_nested(x, comp, exps ...);
    }

    void reset_b()
    {
    	row_b = 0;
    }
    void reset_x_ig()
    {
        row_x_ig = 0;
    }

    void reset(particles_type &part, options_solver opt = options_solver::STANDARD)
    {
    	row = 0;
    	row_b = 0;
        row_x_ig = 0;


        p_map.clear();
    	p_map.resize(part.size_local());

    	A.getMatrixTriplets().clear();

    	construct_pmap(opt);
    }

    void reset_nodec()
    {
    	row = 0;
    	row_b = 0;
        row_x_ig = 0;

    	A.getMatrixTriplets().clear();
    }

    /*! \brief Constructor for the solver
     *
     *
     * \param parts Particle set
     * \param option_solver opt=options_solver::LAGRANGE_MULTIPLIER can be used for purely Neumann system
     *
     */
    DCPSE_scheme(particles_type &part, options_solver opt = options_solver::STANDARD)
            : parts(part), p_map(part.getDecomposition(), 0), row(0), row_b(0), opt(opt) {
        p_map.resize(part.size_local());

        construct_pmap(opt);
    }

    /*DCPSE_scheme(particles_type &part, int option_num)
            : parts(part), p_map(part.getDecomposition(), 0), row(0), row_b(0),opt(options_solver::CUSTOM),offset(option_num) {
        p_map.resize(part.size_local());
        construct_pmap(option_num);
    }*/


    /*! \brief Impose an operator in the Matrix System
    *
    * This function impose an operator on a particular particle region to produce the system
    *
    * Ax = b
    *
    *
    * \param op Operator to impose (A term)
    * \param subset Vector with indices of particles where the operator has to be imposed
    * \param prp_id<>() Property number in the aggregate (Scalar only) for imposing on the RHS b.
    * \param id Equation id in the system that we are imposing given by ed_id type
    *
    */
    template<typename T, typename index_type, unsigned int prp_id>
    void impose(const T &op, openfpm::vector<index_type> &subset,
                const prop_id<prp_id> &num,
                eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        variable_b<prp_id> vb(parts);

        impose_git(op, vb, id.getId(), itd);
    }

    /*! \brief Impose b part only in the Matrix System Ax=b
    *
    * This function impose RHS of an existing Ax=b system.
    *
    *
    * \param subset Vector with indices of particles where the operator has to be imposed
    * \param num right hand side of the term (b term) Constant in this case
    * \param id Equation id in the system that we are imposing given by ed_id type
    *
    */
    template<typename index_type, unsigned int prp_id>
    void impose_b(openfpm::vector<index_type> &subset,
                const prop_id<prp_id> &num,
                eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        variable_b<prp_id> vb(parts);

        impose_git_b(vb, id.getId(), itd);
    }

    /*! \brief Impose x as initial guess for the Matrix System Ax=b
    *
    * This function impose an initial guess for the matrix solver Ax=b
    *
    *
    * \param subset Vector with indices of particles where the operator has to be imposed
    * \param the constant guess num.
    * \param id Equation id in the system that we are imposing given by ed_id type
    *
    */
    template<typename index_type, unsigned int prp_id>
    void impose_x_ig(openfpm::vector<index_type> &subset,
                  const prop_id<prp_id> &num,
                  eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        variable_b<prp_id> vx(parts);

        impose_git_x(vx, id.getId(), itd);
    }

    /*! \brief Impose an operator in the Matrix System
     *
     * This function impose an operator on a particular particle region to produce the system
     *
     * Ax = b
     *
     *
     * \param op Operator to impose (A term)
     * \param subset Vector with indices of particles where the operator has to be imposed
     * \param RHS Expression of the Vector to be imposed
     * \param id Equation id in the system that we are imposing given by ed_id type
     *
     */
    template<typename T, typename index_type, typename RHS_type, typename sfinae = typename std::enable_if<!std::is_fundamental<RHS_type>::type::value>::type>
    void impose(const T &op, openfpm::vector<index_type> &subset,
                const RHS_type &rhs,
                eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        impose_git(op, rhs, id.getId(), itd);
    }
    /*! \brief Impose b part only in the Matrix System Ax=b
    *
    * This function impose RHS of an existing Ax=b system.
    *
    *
    *
    * \param subset Vector with indices of particles where the operator has to be imposed
    * \param Expression of the Vector to be imposed
    * \param id Equation id in the system that we are imposing given by ed_id type
    *
    */
    template<typename index_type, typename RHS_type, typename sfinae = typename std::enable_if<!std::is_fundamental<RHS_type>::type::value>::type>
    void impose_b(openfpm::vector<index_type> &subset,
                const RHS_type &rhs,
                eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();
        impose_git_b(rhs, id.getId(), itd);
    }
    /*! \brief Impose initial guess x in the Matrix System Ax=b
    *
    * This function impose initial guess x of an existing Ax=b system.
    *
    * \param subset Vector with indices of particles where the operator has to be imposed as a guess
    * \param num Constant for all the particles
    * \param id Equation id in the system that we are imposing given by ed_id type
    *
    */
    template<typename index_type, typename RHS_type, typename sfinae = typename std::enable_if<!std::is_fundamental<RHS_type>::type::value>::type>
    void impose_x_ig(openfpm::vector<index_type> &subset,
                  const RHS_type &rhs,
                  eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();
        impose_git_x(rhs, id.getId(), itd);
    }

    /*! \brief Impose an operator in the Matrix System
     *
     * This function impose an operator on a particular particle region to produce the system
     *
     * Ax = b
     *
     *
     * \param op Operator to impose (A term)
     * \param subset Vector with indices of particles where the operator has to be imposed
     * \param num Constant for all the particles
     * \param id Equation id in the system that we are imposing given by ed_id type
     *
     */
    template<typename T, typename index_type>
    void impose(const T &op,
                openfpm::vector<index_type> &subset,
                const typename Sys_eqs::stype num,
                eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        constant_b b(num);

        impose_git(op, b, id.getId(), itd);
    }
    /*! \brief Impose b part only in the Matrix System Ax=b
    *
    * This function impose RHS of an existing Ax=b system.
    *
    *
    *
    * \param subset Vector with indices of particles where the operator has to be imposed
    * \param num Constant for all the particles
    * \param id Equation id in the system that we are imposing given by ed_id type
    *
    */
    template< typename index_type>
    void impose_b(openfpm::vector<index_type> &subset,
                const typename Sys_eqs::stype num,
                eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        constant_b b(num);

        impose_git_b(b, id.getId(), itd);
    }

    /*! \brief Impose initial guess x in the Matrix System Ax=b
*
* This function impose RHS of an existing Ax=b system.
*
*
*
* \param subset Vector with indices of particles where the operator has to be imposed as a guess
* \param num Constant for all the particles
* \param id Equation id in the system that we are imposing given by ed_id type
*
*/
    template< typename index_type>
    void impose_x_ig(openfpm::vector<index_type> &subset,
                  const typename Sys_eqs::stype num,
                  eq_id id = eq_id()) {
        auto itd = subset.template getIteratorElements<0>();

        constant_b x_ig(num);

        impose_git_x(x_ig, id.getId(), itd);
    }

    /*! \brief produce the Matrix
 *
 *  \return the Sparse matrix produced
 *
 */
    template<typename options>
    typename Sys_eqs::SparseMatrix_type &getA(options opt) {
        if (A.isMatrixFilled()) return A;
        if (opt == options_solver::STANDARD) {
            A.resize(tot * Sys_eqs::nvar, tot * Sys_eqs::nvar,
                     p_map.size_local() * Sys_eqs::nvar,
                     p_map.size_local() * Sys_eqs::nvar);
        }
        else if (opt == options_solver::LAGRANGE_MULTIPLIER) {
            auto &v_cl = create_vcluster();
            openfpm::vector<triplet> &trpl = A.getMatrixTriplets();

            if (v_cl.rank() == v_cl.size() - 1) {
                A.resize(Sys_eqs::nvar * (tot + 1), Sys_eqs::nvar * (tot + 1),
                         Sys_eqs::nvar * (p_map.size_local() + 1),
                         Sys_eqs::nvar * (p_map.size_local() + 1));

                for (int j = 0; j < Sys_eqs::nvar; j++) {
                    for (int i = 0; i < tot; i++) {
                        triplet t1;
                        t1.row() = tot * Sys_eqs::nvar + j;
                        t1.col() = i * Sys_eqs::nvar + j;
                        t1.value() = 1;
                        trpl.add(t1);
                    }
                    for (int i = 0; i < p_map.size_local(); i++) {
                        triplet t2;
                        t2.row() = s_pnt * Sys_eqs::nvar + i * Sys_eqs::nvar + j;
                        t2.col() = tot * Sys_eqs::nvar + j;
                        t2.value() = 1;
                        trpl.add(t2);
                    }
                    triplet t3;
                    t3.col() = tot * Sys_eqs::nvar + j;
                    t3.row() = tot * Sys_eqs::nvar + j;
                    t3.value() = 0;
                    trpl.add(t3);
                }
            } else {
                A.resize(Sys_eqs::nvar * (tot + 1), Sys_eqs::nvar * (tot + 1),
                         p_map.size_local() * Sys_eqs::nvar,
                         p_map.size_local() * Sys_eqs::nvar);
                for (int j = 0; j < Sys_eqs::nvar; j++) {
                    for (int i = 0; i < p_map.size_local(); i++) {
                        triplet t2;
                        t2.row() = s_pnt * Sys_eqs::nvar + i * Sys_eqs::nvar + j;
                        t2.col() = tot * Sys_eqs::nvar + j;
                        t2.value() = 1;
                        trpl.add(t2);
                    }
                }
            }
        }
        else{
            auto &v_cl = create_vcluster();
            if (v_cl.rank() == v_cl.size() - 1) {
                A.resize(tot * Sys_eqs::nvar - offset, tot * Sys_eqs::nvar - offset,
                         p_map.size_local() * Sys_eqs::nvar - offset,
                         p_map.size_local() * Sys_eqs::nvar - offset);
            }
            else {
                A.resize(tot * Sys_eqs::nvar - offset, tot * Sys_eqs::nvar - offset,
                         p_map.size_local() * Sys_eqs::nvar,
                         p_map.size_local() * Sys_eqs::nvar);
            }
        }
#ifdef SE_CLASS1
        consistency(opt);
#endif
        return A;

    }

    /*! \brief produce the B vector
     *
     *  \return the vector produced
     *
     */
    typename Sys_eqs::Vector_type &getB(options_solver opt = options_solver::STANDARD) {
/*#ifdef SE_CLASS1
        consistency(opt);
#endif*/
        if (opt == options_solver::LAGRANGE_MULTIPLIER) {
            auto &v_cl = create_vcluster();
            if (v_cl.rank() == v_cl.size() - 1) {
                for(int j=0;j<Sys_eqs::nvar;j++)
                {b(tot * Sys_eqs::nvar+j) = 0;}
            }
        }
        return b;
    }

    /*! \brief produce the B vector
     *
     *  \return the vector produced
     *
     */
    typename Sys_eqs::Vector_type &get_x_ig(options_solver opt = options_solver::STANDARD) {
/*#ifdef SE_CLASS1
        consistency(opt);
#endif*/
        if (opt == options_solver::LAGRANGE_MULTIPLIER) {
            auto &v_cl = create_vcluster();
            if (v_cl.rank() == v_cl.size() - 1) {
                for(int j=0;j<Sys_eqs::nvar;j++)
                    {x_ig(tot * Sys_eqs::nvar+j) = 0;}
            }
        }
        return x_ig;
    }


    template<typename bop, typename iterator>
    void impose_git_b(bop num,
                      long int id,
                      const iterator &it_d) {
        auto it = it_d;
        // iterate all particles points
        while (it.isNext()) {
            // get the particle
            auto key = it.get();
            // Calculate the non-zero colums
            b(p_map.template getProp<0>(key) * Sys_eqs::nvar + id) = num.get(key);
//       std::cout << "b=(" << p_map.template getProp<0>(key)*Sys_eqs::nvar + id << "," << num.get(key)<<")" <<"\n";

            // if SE_CLASS1 is defined check the position
#ifdef SE_CLASS1
            //			T::position(key,gs,s_pos);
#endif
            ++row_b;
            ++it;
        }
    }

    template<typename xop, typename iterator>
    void impose_git_x(xop num,
                      long int id,
                      const iterator &it_d) {
        auto it = it_d;
        // iterate all particles points
        while (it.isNext()) {
            // get the particle
            auto key = it.get();
            // Calculate the non-zero colums
            x_ig(p_map.template getProp<0>(key) * Sys_eqs::nvar + id) = num.get(key);
//       std::cout << "b=(" << p_map.template getProp<0>(key)*Sys_eqs::nvar + id << "," << num.get(key)<<")" <<"\n";

            // if SE_CLASS1 is defined check the position
#ifdef SE_CLASS1
            //			T::position(key,gs,s_pos);
#endif
            ++row_x_ig;
            ++it;
        }
    }


    /*! \brief Impose an operator
     *
     * This function impose an operator on a particular grid region to produce the system
     *
     * Ax = b
     *
     * ## Stokes equation 2D, lid driven cavity with one splipping wall
     * \snippet eq_unit_test.hpp Copy the solution to grid
     *
     * \param op Operator to impose (A term)
     * \param num right hand side of the term (b term)
     * \param id Equation id in the system that we are imposing
     * \param it_d iterator that define where you want to impose
     *
     */
    template<typename T, typename bop, typename iterator>
    void impose_git(const T &op,
                    bop num,
                    long int id,
                    const iterator &it_d) {
        openfpm::vector<triplet> &trpl = A.getMatrixTriplets();

        auto it = it_d;

        //std::unordered_map<long int, typename particles_type::stype> cols;

        tsl::hopscotch_map<long int, typename particles_type::stype> cols;

        // iterate all particles points
        while (it.isNext()) {
            // get the particle
            auto key = it.get();

            // Calculate the non-zero colums
            typename Sys_eqs::stype coeff = 1.0;
            op.template value_nz<Sys_eqs>(p_map, key, cols, coeff, 0);

            // indicate if the diagonal has been set
            bool is_diag = false;

            // create the triplet
            for (auto it2 = cols.begin(); it2 != cols.end(); ++it2) {
                trpl.add();
                trpl.last().row() = p_map.template getProp<0>(key) * Sys_eqs::nvar + id;
                trpl.last().col() = it2->first;
                trpl.last().value() = it2->second;
                if (trpl.last().row() == trpl.last().col())
                {is_diag = true;}
            }

            // If does not have a diagonal entry put it to zero
            if (is_diag == false)
            {
                trpl.add();
                trpl.last().row() = p_map.template getProp<0>(key) * Sys_eqs::nvar + id;
                trpl.last().col() = p_map.template getProp<0>(key) * Sys_eqs::nvar + id;
                trpl.last().value() = 0.0;
            }
            b(p_map.template getProp<0>(key) * Sys_eqs::nvar + id) = num.get(key);
            cols.clear();

            // if SE_CLASS1 is defined check the position
#ifdef SE_CLASS1
            //			T::position(key,gs,s_pos);
#endif

            ++row;
            ++row_b;
            ++row_x_ig;
            ++it;
        }
    }

};

#include "DCPSE/DCPSE_op/EqnsStruct.hpp"
#endif //Eigen
#endif //OPENFPM_PDATA_DCPSE_SOLVER_HPP
