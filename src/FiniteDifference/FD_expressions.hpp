/*
 * FD_expressions.hpp
 *
 *  Created on: May 7, 2020
 *      Author: i-bird
 */

#ifndef FD_EXPRESSIONS_HPP_
#define FD_EXPRESSIONS_HPP_

namespace FD
{

	template <typename exp1,typename exp2, typename impl>
	class grid_dist_expression_op
	{
	};

	/*! \brief Main class that encapsulate a vector properties operand to be used for expressions construction
	 *
	 * \tparam prp property involved
	 * \tparam vector involved
	 *
	 */
	template<unsigned int prp, typename grid>
	class grid_dist_expression
	{
		//! The vector
		grid & g;

	public:

		//! The type of the internal vector
		typedef grid gtype;

		//! Property id of the point
		static const unsigned int prop = prp;

		int var_id = 0;

		void setVarId(int var_id)
		{
			this->var_id = var_id;
		}

		//! constructor for an external vector
		grid_dist_expression(grid & g)
		:g(g)
		{}

		/*! \brief Return the vector on which is acting
		 *
		 * It return the vector used in getVExpr, to get this object
		 *
		 * \return the vector
		 *
		 */
		grid & getGrid()
		{
			return g;
		}

		/*! \brief Return the vector on which is acting
		*
		* It return the vector used in getVExpr, to get this object
		*
		* \return the vector
		*
		*/
		const grid & getGrid() const
		{
			return g;
		}

		/*! \brief This function must be called before value
		 *
		 * it initialize the expression if needed
		 *
		 */
		inline void init() const
		{}

		/*! \brief Evaluate the expression
		 *
		 * \param k where to evaluate the expression
		 *
		 * \return the result of the expression
		 *
		 */
		inline auto value(const grid_dist_key_dx<grid::dims> & k) const -> decltype(g.template getProp<prp>(k))
		{
			return g.template getProp<prp>(k);
		}

		/*! \brief Fill the vector property with the evaluated expression
		 *
		 * \param v_exp expression to evaluate
		 *
		 * \return itself
		 *
		 */
		template<unsigned int prp2> grid & operator=(const grid_dist_expression<prp2,grid> & g_exp)
		{
			g_exp.init();

			auto it = g.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = g_exp.value(key);

				++it;
			}

			return g;
		}

		/*! \brief Fill the vector property with the evaluated expression
		 *
		 * \param v_exp expression to evaluate
		 *
		 * \return itself
		 *
		 */
		template<typename exp1, typename exp2, typename op> grid & operator=(const grid_dist_expression_op<exp1,exp2,op> & g_exp)
		{
			g_exp.init();

			auto it = g.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = g_exp.value(key);

				++it;
			}

			return g;
		}

		/*! \brief Fill the vector property with the double
		 *
		 * \param d value to fill
		 *
		 * \return the internal vector
		 *
		 */
		grid & operator=(double d)
		{
			auto it = g.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = d;

				++it;
			}

			return g;
		}

		template<typename Sys_eqs, typename gmap_type, typename unordered_map_type, typename coeff_type>
		inline void value_nz(const gmap_type & g_map,
				grid_dist_key_dx<Sys_eqs::dims> & key,
				 const grid_sm<Sys_eqs::dims,void> & gs,
				 typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
				 coeff_type & cols,
				 typename Sys_eqs::stype coeff,
				 unsigned int comp) const
		{
			cols[g_map.template getProp<0>(key)*Sys_eqs::nvar + var_id + comp] += coeff;
		}

	/*    inline vector_dist_expression_op<vector_dist_expression<prp,vector>,boost::mpl::int_<1>,VECT_COMP> operator[](int comp)
		{
			int comp_n[1];

			comp_n[0] = comp;

			vector_dist_expression_op<vector_dist_expression<prp,vector>,boost::mpl::int_<1>,VECT_COMP> v_exp(*this,comp_n,var_id);

			return v_exp;
		}*/
	};

	/*! \Create an expression from a vector property
	 *
	 * \tpatam prp property
	 * \param v
	 *
	 */
	template <unsigned int prp,typename grid> inline grid_dist_expression<prp,grid > getV(grid & g)
	{
		grid_dist_expression<prp,grid > exp_g(g);

		return exp_g;
	}

};


#endif /* FD_EXPRESSIONS_HPP_ */
