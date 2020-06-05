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

	/*! \brief Main class that encapsulate a grid properties operand to be used for expressions construction
	 *
	 * \tparam prp property involved
	 * \tparam grid involved
	 *
	 */
	template<unsigned int prp, typename grid>
	class grid_dist_expression
	{
		//! The grid
		grid & g;

	public:

		//! The type of the internal grid
		typedef grid gtype;

		//! Property id of the point
		static const unsigned int prop = prp;

		int var_id = 0;

		void setVarId(int var_id)
		{
			this->var_id = var_id;
		}

		//! constructor for an external grid
		grid_dist_expression(grid & g)
		:g(g)
		{}

		/*! \brief Return the grid on which is acting
		 *
		 * It return the grid used in getVExpr, to get this object
		 *
		 * \return the grid
		 *
		 */
		grid & getGrid()
		{
			return g;
		}

		/*! \brief Return the grid on which is acting
		*
		* It return the grid used in getVExpr, to get this object
		*
		* \return the grid
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

		/*! \brief Fill the grid property with the evaluated expression
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

		/*! \brief Fill the grid property with the evaluated expression
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

		/*! \brief Fill the grid property with the double
		 *
		 * \param d value to fill
		 *
		 * \return the internal grid
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

		template<typename Sys_eqs, typename gmap_type, typename unordered_map_type>
		inline void value_nz(const gmap_type & g_map,
				grid_dist_key_dx<Sys_eqs::dims> & key,
				 const grid_sm<Sys_eqs::dims,void> & gs,
				 typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
				 unordered_map_type & cols,
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


	/*! \brief Main class that encapsulate a double constant
	 *
	 * \param prp no meaning
	 *
	 */
	template<unsigned int dim>
	class grid_dist_expression<dim,double>
	{
		//! constant parameter
		double d;

	public:

		//! constructor from a constant expression
		inline grid_dist_expression(const double & d)
		:d(d)
		{}

		/*! \brief This function must be called before value
		 *
		 * it initialize the expression if needed
		 *
		 */
		inline void init() const
		{}

		/*! \brief Evaluate the expression
		 *
		 * \param k ignored position in the grid
		 *
		 * It just return the value set in the constructor
		 *
		 * \return the constant value
		 *
		 */
		inline double value(const grid_dist_key_dx<dim> & k) const
		{
			return d;
		}


/*
	    template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	    inline void value_nz(pmap_type & p_map, const vect_dist_key_dx & key, unordered_map_type & cols, coeff_type & coeff, unsigned int comp) const
	    {
	        cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += coeff;
	    }*/
	};


	/*! \brief Main class that encapsulate a float constant
	 *
	 * \param prp no meaning
	 *
	 */
	template<unsigned int dim>
	class grid_dist_expression<dim,float>
	{
		//! constant value
		float d;

	public:

		//! type of object the structure return then evaluated
		typedef float vtype;

		//! constrictor from constant value
		inline grid_dist_expression(const float & d)
		:d(d)
		{}

		/*! \brief This function must be called before value
		 *
		 * it initialize the expression if needed
		 *
		 */
		inline void init() const
		{}

		/*! \brief Evaluate the expression
		 *
		 * \param k ignored position in the grid
		 *
		 * It just return the value set in the constructor
		 *
		 * \return the constant value set in the constructor
		 *
		 */
		inline float value(const grid_dist_key_dx<dim> & k) const
		{
			return d;
		}
	};

	struct sum {};

	template <typename exp1, typename exp2>
	class grid_dist_expression_op<exp1,exp2,sum>
	{
		//! expression 1
		const exp1 o1;

		//! expression 1
		const exp2 o2;

	public:

		typedef typename exp1::gtype gtype;

		//! Costruct a FD expression out of two expressions
		inline grid_dist_expression_op(const exp1 & o1,const exp2 & o2)
				:o1(o1),o2(o2)
		{}

		/*! \brief This function must be called before value
		*
		* it initialize the expression if needed
		*
		 */
		inline void init() const
		{
			o1.init();
			o2.init();
		}

		/*! \brief Evaluate the expression
		 *
		 * \param key where to evaluate the expression
		 *
		 * \return the result of the expression
		 *
		 */
		inline auto value(grid_dist_key_dx<gtype::dims> & key) const -> typename std::remove_reference<decltype(o1.value(key))>::type
		{
			typename std::remove_reference<decltype(o1.value(key))>::type val;

			return o1.value(key) + o2.value(key);
		}

		/*! \brief Return the grid on which is acting
		 *
		 * It return the grid used in getVExpr, to get this object
		 *
		 * \return the grid
		 *
		 */
		gtype & getGrid()
		{
			return o1.getGrid();
		}

		/*! \brief Return the grid on which is acting
		*
		* It return the grid used in getVExpr, to get this object
		*
		* \return the grid
		*
		*/
		const gtype & getGrid() const
		{
			return o1.getGrid();
		}
	};


	struct mul {};

	template <typename exp1, typename exp2>
	class grid_dist_expression_op<exp1,exp2,mul>
	{
		//! expression 1
		const exp1 o1;

		//! expression 1
		const exp2 o2;

	public:

		typedef typename exp2::gtype gtype;

		//! Costruct a FD expression out of two expressions
		inline grid_dist_expression_op(const exp1 & o1,const exp2 & o2)
				:o1(o1),o2(o2)
		{}

		/*! \brief This function must be called before value
		*
		* it initialize the expression if needed
		*
		 */
		inline void init() const
		{
			o1.init();
			o2.init();
		}

		/*! \brief Evaluate the expression
		 *
		 * \param key where to evaluate the expression
		 *
		 * \return the result of the expression
		 *
		 */
		inline auto value(grid_dist_key_dx<gtype::dims> & key) const -> typename std::remove_reference<decltype(o1.value(key))>::type
		{
			typename std::remove_reference<decltype(o1.value(key))>::type val;

			return o1.value(key) * o2.value(key);
		}

		/*! \brief Return the grid on which is acting
		 *
		 * It return the grid used in getVExpr, to get this object
		 *
		 * \return the grid
		 *
		 */
		gtype & getGrid()
		{
			return o1.getGrid();
		}

		/*! \brief Return the grid on which is acting
		*
		* It return the grid used in getVExpr, to get this object
		*
		* \return the grid
		*
		*/
		const gtype & getGrid() const
		{
			return o1.getGrid();
		}

		template<typename Sys_eqs, typename gmap_type, typename unordered_map_type>
		inline void value_nz(const gmap_type & g_map,
				grid_dist_key_dx<Sys_eqs::dims> & key,
				 const grid_sm<Sys_eqs::dims,void> & gs,
				 typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
				 unordered_map_type & cols,
				 typename Sys_eqs::stype coeff,
				 unsigned int comp) const
		{
			typename Sys_eqs::stype coeff_tmp = o1.value(key) * coeff;

			o2.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,coeff_tmp,comp);
		}
	};

	/*! \Create an expression from a grid property
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

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, unsigned int p2, typename g1, typename g2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,g1>,FD::grid_dist_expression<p2,g2>,FD::sum>
operator+(const FD::grid_dist_expression<p1,g1> & ga, const FD::grid_dist_expression<p2,g2> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,g1>,FD::grid_dist_expression<p2,g2>,FD::sum> exp_sum(ga,gb);

	return exp_sum;
}

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, typename op1, unsigned int prp1, typename g1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<prp1,g1>,FD::sum>
operator+(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, const FD::grid_dist_expression<prp1,g1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<prp1,g1>,FD::sum> exp_sum(ga,gb);

	return exp_sum;
}

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, typename op1, unsigned int prp1, typename g1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::sum>
operator+(const FD::grid_dist_expression<prp1,g1> & ga, const FD::grid_dist_expression_op<exp1,exp2,op1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::sum> exp_sum(ga,gb);

	return exp_sum;
}

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, typename op1, typename exp3 , typename exp4, typename op2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression_op<exp3,exp4,op2>,FD::sum>
operator+(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, const FD::grid_dist_expression_op<exp3,exp4,op2> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression_op<exp3,exp4,op2>,FD::sum> exp_sum(ga,gb);

	return exp_sum;
}

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename g1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1>,FD::grid_dist_expression<g1::dims,double>,FD::sum>
operator+(const FD::grid_dist_expression<prp1,g1> & ga, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1>,FD::grid_dist_expression<g1::dims,double>,FD::sum> exp_sum(ga,FD::grid_dist_expression<g1::dims,double>(d));

	return exp_sum;
}

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int prp1 , typename g1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<g1::dims,double>,FD::grid_dist_expression<prp1,g1>,FD::sum>
operator+(double d, const FD::grid_dist_expression<prp1,g1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<g1::dims,double>,FD::grid_dist_expression<prp1,g1>,FD::sum> exp_sum(FD::grid_dist_expression<g1::dims,double>(d),gb);

	return exp_sum;
}

/* \brief sum two distributed grid expression
 *
 * \param ga grid expression one
 * \param gb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, typename op1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double>,FD::sum>
operator+(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double>,FD::sum> exp_sum(ga,FD::grid_dist_expression<exp1::gtype::dims,double>(d));

	return exp_sum;
}



/* \brief Multiply two distributed grid expression
 *
 * \param va grid expression one
 * \param vb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p2, typename g2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<g2::dims,double>,FD::grid_dist_expression<p2,g2>,FD::mul>
operator*(double d, const FD::grid_dist_expression<p2,g2> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<g2::dims,double>,FD::grid_dist_expression<p2,g2>,FD::mul> exp_sum(FD::grid_dist_expression<g2::dims,double>(d),vb);

	return exp_sum;
}

/* \brief Multiply two distributed grid expression
 *
 * \param va grid expression one
 * \param vb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p2, typename g2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p2,g2>,FD::grid_dist_expression<g2::dims,double>,FD::mul>
operator*(const FD::grid_dist_expression<p2,g2> & va, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p2,g2>,FD::grid_dist_expression<g2::dims,double>,FD::mul> exp_sum(va,FD::grid_dist_expression<g2::dims,double>(d));

	return exp_sum;
}

/* \brief Multiply two distributed grid expression
 *
 * \param va grid expression one
 * \param vb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1,unsigned int p2, typename v2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1>,FD::grid_dist_expression<p2,v2>,FD::mul>
operator*(const FD::grid_dist_expression<p1,v1> & va, const FD::grid_dist_expression<p2,v2> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1>,FD::grid_dist_expression<p2,v2>,FD::mul> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two distributed grid expression
 *
 * \param va grid expression one
 * \param vb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1, typename exp1, typename exp2, typename op1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul>
operator*(const FD::grid_dist_expression<p1,v1> & va, const FD::grid_dist_expression_op<exp1,exp2,op1> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two distributed grid expression
 *
 * \param va grid expression one
 * \param vb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int p1, typename v1, typename exp1, typename exp2, typename op1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<p1,v1>,FD::mul>
operator*(const FD::grid_dist_expression_op<exp1,exp2,op1> & va, const FD::grid_dist_expression<p1,v1> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<p1,v1>,FD::mul> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two distributed grid expression
 *
 * \param va grid expression one
 * \param vb grid expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1, typename exp2, typename op1, typename exp3 , typename exp4, typename op2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression_op<exp3,exp4,op2>,FD::mul>
operator*(const FD::grid_dist_expression_op<exp1,exp2,op1> & va, const FD::grid_dist_expression_op<exp3,exp4,op2> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression_op<exp3,exp4,op2>,FD::mul> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply a distributed grid expression by a number
 *
 * \param va grid expression
 * \param d number
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, typename op1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double>,FD::mul>
operator*(const FD::grid_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double>,FD::mul> exp_sum(va,FD::grid_dist_expression<exp1::gtype::dims,double>(d));

	return exp_sum;
}

/* \brief Multiply a distributed grid expression by a number
 *
 * \param d number
 * \param vb grid expression
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2, typename op1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<exp1::gtype::dims,double>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul>
operator*(double d, const FD::grid_dist_expression_op<exp1,exp2,op1> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<exp1::gtype::dims,double>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul> exp_sum(FD::grid_dist_expression<exp1::gtype::dims,double>(d),vb);

	return exp_sum;
}

#endif /* FD_EXPRESSIONS_HPP_ */
