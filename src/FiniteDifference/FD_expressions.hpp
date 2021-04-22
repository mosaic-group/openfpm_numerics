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

	constexpr int NORM_EXPRESSION = 0;
	constexpr int STAG_EXPRESSION = 1;
	constexpr int GRID_COMP = 2;

	template <typename exp1,typename exp2, typename impl>
	class grid_dist_expression_op
	{
	};

	struct g_comp {};

	template<unsigned int prp, typename grid, unsigned int impl>
	class grid_dist_expression
	{};

	/*! \brief Main class that encapsulate a grid properties operand to be used for expressions construction
	 *
	 * \tparam prp property involved
	 * \tparam grid involved
	 *
	 */
	template<unsigned int prp, typename grid>
	class grid_dist_expression<prp,grid,NORM_EXPRESSION>
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
		inline auto value(const grid_dist_key_dx<grid::dims> & k, comb<grid::dims> & c_where, int comp = 0) const -> decltype(g.template getProp<prp>(k))
		{
			return g.template getProp<prp>(k);
		}

		/*! \brief Evaluate the expression
		 *
		 * \param k where to evaluate the expression
		 *
		 * \return the result of the expression
		 *
		 */
		inline auto value_ref(const grid_dist_key_dx<grid::dims> & k, comb<grid::dims> & c_where) const -> decltype(g.template getProp<prp>(k))
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
		template<unsigned int prp2> grid & operator=(const grid_dist_expression<prp2,grid,NORM_EXPRESSION> & g_exp)
		{
			g_exp.init();

			comb<grid::dims> s_pos;
			s_pos.zero;

			auto it = g.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = g_exp.value(key,s_pos);

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

			comb<grid::dims> s_pos;
			s_pos.mone();

			auto it = g.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = g_exp.value(key,s_pos);

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

			//comb<grid::dims> s_pos;
			//s_pos.mone();

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
				 unsigned int comp,
				 comb<Sys_eqs::dims> & c_where) const
		{
			cols[g_map.template getProp<0>(key)*Sys_eqs::nvar + var_id + comp] += coeff;
		}

	    inline grid_dist_expression_op<grid_dist_expression<prp,grid,NORM_EXPRESSION>,boost::mpl::int_<1>,g_comp> operator[](int comp)
		{
			int comp_n[1];

			comp_n[0] = comp;

			grid_dist_expression_op<grid_dist_expression<prp,grid,NORM_EXPRESSION>,boost::mpl::int_<1>,g_comp> v_exp(*this,comp_n,var_id);

			return v_exp;
		}

        //Need more treatment for staggered (c_where based on exp)
        inline typename gtype::stype get(grid_dist_key_dx<gtype::dims> & key)
        {
            comb<gtype::dims> c_where;
            c_where.zero();
            return this->value(key,c_where);
        }

		int isConstant(){
		    return false;
		}
	};

	template<unsigned int prp, typename grid>
	class grid_dist_expression<prp,grid,STAG_EXPRESSION>
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
		inline auto value_ref(const grid_dist_key_dx<grid::dims> & k, comb<grid::dims> & c_where) const -> decltype(g.template getProp<prp>(k))
		{
			return g.template getProp<prp>(k);
		}

		/*! \brief Evaluate the expression
		 *
		 * \param k where to evaluate the expression
		 *
		 * \return the result of the expression
		 *
		 */
		inline auto value(grid_dist_key_dx<grid::dims> & k, comb<grid::dims> & c_where, int comp = 0) const -> typename std::remove_reference<decltype(g.template getProp<prp>(k))>::type
		{
			comb<grid::dims> c_o1 = g. getStagPositions()[prp].get(comp);

            // x0, dx are defined in proper dir є(x, y, z)

            int c = 0;

            typename std::remove_reference<decltype(g.template getProp<prp>(k))>::type inte = 0;
        	for (int i = 0 ; i < grid::dims ; i++)
        	{
        		if (c_where[i] != c_o1[i])
        		{
        			int sign = (c_where[i] > c_o1[i])?1:-1;

        			inte += g.template getProp<prp>(k);

        			long int x0 = k.getKeyRef().get(i);

        			k.getKeyRef().set_d(i, x0 + sign);
        			inte += g.template getProp<prp>(k);
        			k.getKeyRef().set_d(i, x0);
        			c += 2;
        		}
        	}

        	if (c != 0)
        	{inte /= c;}
        	else
        	{inte = g.template getProp<prp>(k);}

			return inte;
		}

		/*! \brief Fill the grid property with the evaluated expression
		 *
		 * \param v_exp expression to evaluate
		 *
		 * \return itself
		 *
		 */
		template<unsigned int prp2> grid & operator=(const grid_dist_expression<prp2,grid,STAG_EXPRESSION> & g_exp)
		{
			g_exp.init();

			auto it = g.getDomainIterator();

			comb<grid::dims> s_pos = g.getStagPosition()[prp].get(0);

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = g_exp.value(key,s_pos);

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

			comb<grid::dims> s_pos = g.getStagPositions()[prp].get(0);

			while (it.isNext())
			{
				auto key = it.get();

				g.template getProp<prp>(key) = g_exp.value(key,s_pos);

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
				 unsigned int comp,
				 comb<Sys_eqs::dims> & c_where) const
		{
			comb<grid::dims> c_o1 = g.getStagPositions()[prp].get(comp);

            // x0, dx are defined in proper dir є(x, y, z)

            int c = 0;
        	for (int i = 0 ; i < grid::dims ; i++)
        	{
        		if (c_where[i] != c_o1[i])
        		{c += 2;}
        	}

        	if (c != 0)
        	{
            	for (int i = 0 ; i < grid::dims ; i++)
            	{
            		if (c_where[i] != c_o1[i])
            		{
            			int sign = (c_where[i] > c_o1[i])?1:-1;

            			cols[g_map.template getProp<0>(key)*Sys_eqs::nvar + var_id + comp] += coeff / c;

            			long int x0 = key.getKeyRef().get(i);

            			key.getKeyRef().set_d(i, x0 + sign);
            			cols[g_map.template getProp<0>(key)*Sys_eqs::nvar + var_id + comp] += coeff / c;
            			key.getKeyRef().set_d(i, x0);
            		}
            	}
        	}
        	else
        	{cols[g_map.template getProp<0>(key)*Sys_eqs::nvar + var_id + comp] += coeff;}
		}

	    inline grid_dist_expression_op<grid_dist_expression<prp,grid,STAG_EXPRESSION>,boost::mpl::int_<1>,g_comp> operator[](int comp)
		{
			int comp_n[1];

			comp_n[0] = comp;

			grid_dist_expression_op<grid_dist_expression<prp,grid,STAG_EXPRESSION>,boost::mpl::int_<1>,g_comp> v_exp(*this,comp_n,var_id);

			return v_exp;
		}

	};

	/*! \brief Main class that encapsulate a double constant
	 *
	 * \param prp no meaning
	 *
	 */
	template<unsigned int dim>
	class grid_dist_expression<dim,double,NORM_EXPRESSION>
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
		inline double value(const grid_dist_key_dx<dim> & k, comb<dim> & c_where) const
		{
			return d;
		}
	};

	/*! \brief Main class that encapsulate a double constant
	 *
	 * \param prp no meaning
	 *
	 */
	template<unsigned int dim>
	class grid_dist_expression<dim,double,STAG_EXPRESSION>
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
		inline double value(const grid_dist_key_dx<dim> & k, comb<dim> & c_where) const
		{
			return d;
		}
	};

	/*! \brief Main class that encapsulate a float constant
	 *
	 * \param prp no meaning
	 *
	 */
	template<unsigned int dim, unsigned int impl>
	class grid_dist_expression<dim,float,impl>
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
		inline auto value(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> typename std::remove_reference<decltype(o1.value(key,c_where))>::type
		{
			typename std::remove_reference<decltype(o1.value(key,c_where))>::type val;

			return o1.value(key,c_where) + o2.value(key,c_where);
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
                             unsigned int comp,
                             comb<Sys_eqs::dims> & c_where) const
        {
            o1.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,coeff,comp,c_where);
            o2.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,coeff,comp,c_where);
        }
	};

	struct sub {};


    template <typename exp1, typename exp2>
	class grid_dist_expression_op<exp1,exp2,sub>
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
		inline auto value(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> typename std::remove_reference<decltype(o1.value(key,c_where))>::type
		{
			typename std::remove_reference<decltype(o1.value(key,c_where))>::type val;

			return o1.value(key,c_where) - o2.value(key,c_where);
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
                             unsigned int comp,
                             comb<Sys_eqs::dims> & c_where) const
        {
            o1.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,coeff,comp,c_where);
            o2.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,-coeff,comp,c_where);
        }
	};

    struct subuni{};

    template <typename exp1>
    class grid_dist_expression_op<exp1,void,subuni>
    {
        //! expression 1
        const exp1 o1;

    public:

        typedef typename exp1::gtype gtype;

        //! Costruct a FD expression out of two expressions
        inline grid_dist_expression_op(const exp1 & o1)
                :o1(o1)
        {}

        /*! \brief This function must be called before value
        *
        * it initialize the expression if needed
        *
         */
        inline void init() const
        {
            o1.init();
        }

        /*! \brief Evaluate the expression
         *
         * \param key where to evaluate the expression
         *
         * \return the result of the expression
         *
         */
        inline auto value(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> typename std::remove_reference<decltype(o1.value(key,c_where))>::type
        {
            typename std::remove_reference<decltype(o1.value(key,c_where))>::type val;

            return -o1.value(key,c_where);
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
                             unsigned int comp,
                             comb<Sys_eqs::dims> & c_where) const
        {
            o1.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,-coeff,comp,c_where);
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
		inline auto value(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> typename std::remove_reference<decltype(o1.value(key,c_where))>::type
		{
			typename std::remove_reference<decltype(o1.value(key,c_where))>::type val;

			return o1.value(key,c_where) * o2.value(key,c_where);
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
				 unsigned int comp,
				 comb<Sys_eqs::dims> & c_where) const
		{
			typename Sys_eqs::stype coeff_tmp = o1.value(key,c_where) * coeff;

			o2.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,coeff_tmp,comp,c_where);
		}
	};

	/*! \brief selector for position or properties left side expression
	 *
	 * \tparam vector type of the original vector
	 *
	 * \tparam prp property id
	 *
	 */
	template <typename grid_type, unsigned int prp>
	struct pos_or_propL
	{
		//! return the value (position or property) of the particle k in the vector v
		static inline auto value(grid_type & v, const grid_dist_key_dx<grid_type::dims> & k) -> decltype(v.template getProp<prp>(k))
		{
			return v.template getProp<prp>(k);
		}
	};

	template<unsigned int, bool is_valid>
	struct get_grid_dist_expression_op
	{
		template<typename exp_type>
		inline static auto get(exp_type & o1, grid_dist_key_dx<exp_type::gtype::dims> & key, comb<exp_type::gtype::dims> & c_where) -> decltype(o1.value(key,c_where))
		{
			return o1.value(key,c_where);
		}

		template<typename exp_type>
		inline static auto get_ref(exp_type & o1, grid_dist_key_dx<exp_type::gtype::dims> & key, comb<exp_type::gtype::dims> & c_where) -> decltype(o1.value_ref(key,c_where))
		{
			return o1.value_ref(key,c_where);
		}

		template<unsigned int prop, typename exp_type, typename grid_type>
		inline static void assign(exp_type & o1, grid_type & g, const grid_dist_key_dx<exp_type::gtype::dims> & key)
		{
			pos_or_propL<grid_type,exp_type::prop>::value(g,key) = o1.value(key);
		}

		template<unsigned int prop, typename grid_type>
		inline static void assign_double(double d, grid_type & g, const grid_dist_key_dx<grid_type::dims> & key)
		{
			pos_or_propL<grid_type,prop>::value(g,key) = d;
		}
	};

	template<>
	struct get_grid_dist_expression_op<1,false>
	{
		template<typename exp_type>
		static int get(exp_type & o1, grid_dist_key_dx<exp_type::gtype::dims> & key, comb<exp_type::gtype::dims> & c_where, const int (& comp)[1])
		{
			return 0;
		}

		template<unsigned int prop, typename exp_type, typename grid_type>
		inline static void assign(exp_type & o1, grid_type & g, const grid_dist_key_dx<exp_type::gtype::dims> & key)
		{
		}

		template<unsigned int prop, typename grid_type>
		inline static void assign_double(double d, grid_type & g, const grid_dist_key_dx<grid_type::dims> & key)
		{
		}
	};

	template<>
	struct get_grid_dist_expression_op<1,true>
	{
		template<typename exp_type>
		static auto get(exp_type & o1, grid_dist_key_dx<exp_type::gtype::dims> & key, comb<exp_type::gtype::dims> & c_where, const int (& comp)[1]) -> decltype(o1.value(key,c_where)[0])
		{
			return o1.value(key,c_where,comp[0])[comp[0]];
		}

		template<typename exp_type>
		static auto get_ref(exp_type & o1, grid_dist_key_dx<exp_type::gtype::dims> & key, comb<exp_type::gtype::dims> & c_where, const int (& comp)[1]) -> decltype(o1.value_ref(key,c_where)[0])
		{
			return o1.value_ref(key,c_where)[comp[0]];
		}

		template<unsigned int prop,typename exp_type, typename grid_type>
		inline static void assign(exp_type & o1, grid_type & g, grid_dist_key_dx<exp_type::gtype::dims> & key,comb<exp_type::gtype::dims> & c_where, const int (& comp)[1])
		{
			pos_or_propL<grid_type,prop>::value(g,key)[comp[0]] = o1.value(key,c_where);
		}

		template<unsigned int prop, typename grid_type>
		inline static void assign_double(double d, grid_type & g, const grid_dist_key_dx<grid_type::dims> & key, const int (& comp)[1])
		{
			pos_or_propL<grid_type,prop>::value(g,key)[comp[0]] = d;
		}
	};

	template<>
	struct get_grid_dist_expression_op<2,true>
	{
		template<typename exp_type>
		static auto get(exp_type & o1, grid_dist_key_dx<exp_type::gtype::dims> & key, comb<exp_type::gtype::dims> & c_where, const int (& comp)[2]) -> decltype(o1.value(key,c_where)[0][0])
		{
			return o1.value(key)[comp[0]][comp[1]];
		}

		template<unsigned int prop,typename exp_type, typename grid_type>
		inline static void assign(exp_type & o1, grid_type & g, const grid_dist_key_dx<grid_type::dims> & key, const int (& comp)[2])
		{
			pos_or_propL<grid_type,prop>::value(g,key)[comp[0]][comp[1]] = o1.value(key);
		}

		template<unsigned int prop, typename grid_type>
		inline static void assign_double(double d, grid_type & g, const grid_dist_key_dx<grid_type::dims> & key, const int (& comp)[2])
		{
			pos_or_propL<grid_type,prop>::value(g,key)[comp[0]][comp[1]] = d;
		}
	};

	/*! \brief it take an expression and create the negatove of this expression
	 *
	 *
	 */
	template <typename exp1,int n>
	class grid_dist_expression_op<exp1,boost::mpl::int_<n>,g_comp>
	{
		//! expression 1
		exp1 o1;

		//! component
		int comp[n];

		int var_id = 0;
	    void setVarId(int var_id)
	    {
	        this->var_id = var_id;
	    }

		typedef grid_dist_expression_op<exp1,boost::mpl::int_<n>,g_comp> myself;

	public:

	        typedef std::false_type is_ker;

		typedef typename exp1::gtype gtype;

		//! Property id of the point
		static const unsigned int prop = exp1::prop;

		//! constructor from an expresssion

		grid_dist_expression_op(const exp1 & o1, int (& comp)[n], int var_id)
		:o1(o1),var_id(var_id)
		{
			for (int i = 0 ; i < n ; i++)
			{this->comp[i] = comp[i];}
		}

	    /*! \brief Return the vector on which is acting
	    *
	    * It return the vector used in getVExpr, to get this object
	    *
	    * \return the vector
	    *
	    */
	    const gtype & getGrid() const
	    {
	        return o1.getGrid();
	    }

	    /*! \brief Return the vector on which is acting
	    *
	    * It return the vector used in getVExpr, to get this object
	    *
	    * \return the vector
	    *
	    */
	    gtype & getGrid()
	    {
	        return o1.getGrid();
	    }

		//! initialize the expression tree
		inline void init() const
		{
			o1.init();
		}

		//! property on which this view is acting
		typedef typename boost::mpl::at<typename gtype::value_type::type,boost::mpl::int_<exp1::prop>>::type property_act;

		/*! \brief Return the result of the expression
		 *
		 * \note this function must be deactivated on transitional objects. Suppose we are slicing a tensor of rank 2
		 *            an object of rank 1 is implicitly created for such object we have to deactivate this function
		 *            because ill-formed
		 *
		 * \param key point where to evaluate
		 *
		 *
		 */
		inline auto value(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> decltype(get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,key,c_where,comp))
		{
			return get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,key,c_where,comp);
		}

		/*! \brief Return the result of the expression
		 *
		 * \note this function must be deactivated on transitional objects. Suppose we are slicing a tensor of rank 2
		 *            an object of rank 1 is implicitly created for such object we have to deactivate this function
		 *            because ill-formed
		 *
		 * \param key point where to evaluate
		 *
		 *
		 */
		inline auto value_ref(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> decltype(get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get(o1,key,c_where,comp))
		{
			return get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::get_ref(o1,key,c_where,comp);
		}

		/*! \brief Return the result of the expression
		 *
		 * \note this function must be deactivated on transitional objects. Suppose we are slicing a tensor of rank 2
		 *            an object of rank 1 is implicitly created for such object we have to deactivate this function
		 *            because ill-formed
		 *
		 * \param key point where to evaluate
		 *
		 *
		 */
		inline auto get(grid_dist_key_dx<gtype::dims> & key, comb<gtype::dims> & c_where) const -> decltype(value(key,c_where))
		{
			return this->value(key,c_where);
		}

		template<typename Sys_eqs, typename gmap_type, typename unordered_map_type>
		inline void value_nz(const gmap_type & g_map,
				grid_dist_key_dx<Sys_eqs::dims> & key,
				 const grid_sm<Sys_eqs::dims,void> & gs,
				 typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
				 unordered_map_type & cols,
				 typename Sys_eqs::stype coeff,
				 unsigned int comp_,
				 comb<Sys_eqs::dims> & c_where) const
	    {
	#ifdef SE_CLASS1

	    	if (n != 1)
	    	{
	    		std::cout << __FILE__ << ":" << __LINE__ << " Error it only work for tensore of rank 1 ... like vectors " << std::endl;
	    	}

	#endif

	        o1.template value_nz<Sys_eqs>(g_map,key,gs,spacing,cols,coeff,comp_ + var_id + comp[0],c_where);
	    }

	    inline grid_dist_expression_op<exp1,boost::mpl::int_<2>,g_comp> operator[](int comp_)
	    {
	    	int comp_n[n+1];

	    	for (int i = 0 ; i < n ; i++)
	    	{comp_n[i] = comp[i];}
	    	comp_n[n] = comp_;

	    	grid_dist_expression_op<exp1,boost::mpl::int_<2>,g_comp> v_exp(o1,comp_n,var_id);

	    	return v_exp;
	    }

	    //Need more treatment for staggered (c_where based on exp)
        inline typename gtype::stype get(grid_dist_key_dx<gtype::dims> & key)
        {
		    comb<gtype::dims> c_where;
		    c_where.zero();
            return this->value(key,c_where);
        }


        /*! \brief Fill the vector property with the evaluated expression
		 *
		 * \param v_exp expression to evaluate
		 *
		 * \return itself
		 *
		 */
		template<unsigned int prp2, unsigned int impl> gtype & operator=(const grid_dist_expression<prp2,gtype,impl> & v_exp)
		{
			v_exp.init();

			auto & g = getGrid();

			auto it = g.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<exp1::prop>(v_exp,g,key,comp);

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
		template<typename exp1_, typename exp2_, typename op> gtype & operator=(const grid_dist_expression_op<exp1_,exp2_,op> & v_exp)
		{
			v_exp.init();

			auto & g = getGrid();

			auto it = g.getDomainIterator();

            comb<gtype::dims> c_where;
            c_where.zero();

			while (it.isNext())
			{
				auto key = it.get();

				get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign<exp1::prop>(v_exp,g,key,c_where,comp);

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
		gtype & operator=(double d)
		{
			auto & v = getGrid();

			auto it = v.getDomainIterator();

			while (it.isNext())
			{
				auto key = it.get();

				//pos_or_propL<vtype,exp1::prp>::value(v,key) = d;
				get_grid_dist_expression_op<n,n == rank_gen<property_act>::type::value>::template assign_double<exp1::prop>(d,v,key,comp);


				++it;
			}

			return v;
		}

		int isConstant(){
		    return false;
		}
	};


	/*! \Create an expression from a grid property
	 *
	 * \tpatam prp property
	 * \param v
	 *
	 */
	template <unsigned int prp,typename grid> inline grid_dist_expression<prp,grid,NORM_EXPRESSION> getV(grid & g)
	{
		grid_dist_expression<prp,grid,NORM_EXPRESSION> exp_g(g);

		return exp_g;
	}

	/*! \Create an expression from a grid property
	 *
	 * \tpatam prp property
	 * \param v
	 *
	 */
	template <unsigned int prp, typename grid> inline grid_dist_expression<prp,grid,STAG_EXPRESSION> getV_stag(grid & g)
	{
		grid_dist_expression<prp,grid,STAG_EXPRESSION> exp_g(g);

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
template<unsigned int p1, unsigned int p2, typename g1, typename g2, unsigned int impl_p1, unsigned int impl_p2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,g1,impl_p1>,FD::grid_dist_expression<p2,g2,impl_p2>,FD::sum>
operator+(const FD::grid_dist_expression<p1,g1,impl_p1> & ga, const FD::grid_dist_expression<p2,g2,impl_p2> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,g1,impl_p1>,FD::grid_dist_expression<p2,g2,impl_p2>,FD::sum> exp_sum(ga,gb);

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
template<typename exp1 , typename exp2, typename op1, unsigned int prp1, typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sum>
operator+(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, const FD::grid_dist_expression<prp1,g1,impl_p1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sum> exp_sum(ga,gb);

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
template<typename exp1 , typename exp2, typename op1, unsigned int prp1, typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::sum>
operator+(const FD::grid_dist_expression<prp1,g1,impl_p1> & ga, const FD::grid_dist_expression_op<exp1,exp2,op1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::sum> exp_sum(ga,gb);

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
template<unsigned int prp1 , typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::sum>
operator+(const FD::grid_dist_expression<prp1,g1,impl_p1> & ga, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::sum> exp_sum(ga,FD::grid_dist_expression<g1::dims,double,impl_p1>(d));

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
template<unsigned int prp1 , typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sum>
operator+(double d, const FD::grid_dist_expression<prp1,g1,impl_p1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sum> exp_sum(FD::grid_dist_expression<g1::dims,double,impl_p1>(d),gb);

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
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::sum>
operator+(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::sum> exp_sum(ga,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>(d));

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
template<unsigned int p1, unsigned int p2, typename g1, typename g2, unsigned int impl_p1, unsigned int impl_p2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,g1,impl_p1>,FD::grid_dist_expression<p2,g2,impl_p2>,FD::sub>
operator-(const FD::grid_dist_expression<p1,g1,impl_p1> & ga, const FD::grid_dist_expression<p2,g2,impl_p2> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,g1,impl_p1>,FD::grid_dist_expression<p2,g2,impl_p2>,FD::sub> exp_sum(ga,gb);

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
template<typename exp1 , typename exp2, typename op1, unsigned int prp1, typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sub>
operator-(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, const FD::grid_dist_expression<prp1,g1,impl_p1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sub> exp_sum(ga,gb);

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
template<typename exp1 , typename exp2, typename op1, unsigned int prp1, typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::sub>
operator-(const FD::grid_dist_expression<prp1,g1,impl_p1> & ga, const FD::grid_dist_expression_op<exp1,exp2,op1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::sub> exp_sum(ga,gb);

	return exp_sum;
}

/* \brief minus of a distributed grid expression
 *
 * \param ga grid expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename exp1 , typename exp2_, typename op1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2_,op1>,void,FD::subuni>
operator-(const FD::grid_dist_expression_op<exp1,exp2_,op1> & ga)
{
    FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2_,op1>,void,FD::subuni> exp_sum(ga);

    return exp_sum;
}

template<unsigned int prp1 , typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,void,FD::subuni>
operator-(const FD::grid_dist_expression<prp1,g1,impl_p1> & ga)
{
    FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,void,FD::subuni> exp_sum(ga);

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
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression_op<exp3,exp4,op2>,FD::sub>
operator-(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, const FD::grid_dist_expression_op<exp3,exp4,op2> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression_op<exp3,exp4,op2>,FD::sub> exp_sum(ga,gb);

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
template<unsigned int prp1 , typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::sub>
operator-(const FD::grid_dist_expression<prp1,g1,impl_p1> & ga, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<prp1,g1,impl_p1>,FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::sub> exp_sum(ga,FD::grid_dist_expression<g1::dims,double,impl_p1>(d));

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
template<unsigned int prp1 , typename g1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sub>
operator-(double d, const FD::grid_dist_expression<prp1,g1,impl_p1> & gb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<g1::dims,double,impl_p1>,FD::grid_dist_expression<prp1,g1,impl_p1>,FD::sub> exp_sum(FD::grid_dist_expression<g1::dims,double,impl_p1>(d),gb);

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
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::sub>
operator-(const FD::grid_dist_expression_op<exp1,exp2,op1> & ga, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::sub> exp_sum(ga,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>(d));

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
template<unsigned int p2, typename g2, unsigned int impl_p2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<g2::dims,double,impl_p2>,FD::grid_dist_expression<p2,g2,impl_p2>,FD::mul>
operator*(double d, const FD::grid_dist_expression<p2,g2,impl_p2> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<g2::dims,double,impl_p2>,FD::grid_dist_expression<p2,g2,impl_p2>,FD::mul> exp_sum(FD::grid_dist_expression<g2::dims,double,impl_p2>(d),vb);

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
template<unsigned int p2, typename g2, unsigned int impl_p2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p2,g2,impl_p2>,FD::grid_dist_expression<g2::dims,double,impl_p2>,FD::mul>
operator*(const FD::grid_dist_expression<p2,g2,impl_p2> & va, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p2,g2,impl_p2>,FD::grid_dist_expression<g2::dims,double,impl_p2>,FD::mul> exp_sum(va,FD::grid_dist_expression<g2::dims,double,impl_p2>(d));

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
template<unsigned int p1, typename v1,unsigned int p2, typename v2, unsigned int impl_p1, unsigned int impl_p2>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1,impl_p1>,FD::grid_dist_expression<p2,v2,impl_p2>,FD::mul>
operator*(const FD::grid_dist_expression<p1,v1, impl_p1> & va, const FD::grid_dist_expression<p2,v2,impl_p2> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1,impl_p1>,FD::grid_dist_expression<p2,v2,impl_p2>,FD::mul> exp_sum(va,vb);

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
template<unsigned int p1, typename v1, typename exp1, typename exp2, typename op1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1,impl_p1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul>
operator*(const FD::grid_dist_expression<p1,v1,impl_p1> & va, const FD::grid_dist_expression_op<exp1,exp2,op1> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<p1,v1,impl_p1>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul> exp_sum(va,vb);

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
template<unsigned int p1, typename v1, typename exp1, typename exp2, typename op1, unsigned int impl_p1>
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<p1,v1,impl_p1>,FD::mul>
operator*(const FD::grid_dist_expression_op<exp1,exp2,op1> & va, const FD::grid_dist_expression<p1,v1,impl_p1> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<p1,v1,impl_p1>,FD::mul> exp_sum(va,vb);

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
inline FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::mul>
operator*(const FD::grid_dist_expression_op<exp1,exp2,op1> & va, double d)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression_op<exp1,exp2,op1>,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::mul> exp_sum(va,FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>(d));

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
inline FD::grid_dist_expression_op<FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul>
operator*(double d, const FD::grid_dist_expression_op<exp1,exp2,op1> & vb)
{
	FD::grid_dist_expression_op<FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>,FD::grid_dist_expression_op<exp1,exp2,op1>,FD::mul> exp_sum(FD::grid_dist_expression<exp1::gtype::dims,double,FD::NORM_EXPRESSION>(d),vb);

	return exp_sum;
}

#endif /* FD_EXPRESSIONS_HPP_ */
