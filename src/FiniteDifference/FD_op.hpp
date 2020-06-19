/*
 * FD_op.hpp
 *
 *  Created on: May 1, 2020
 *      Author: Abhinav Singh
 */

#ifndef FD_OP_HPP_
#define FD_OP_HPP_
#include "util/common.hpp"
#include "FD_expressions.hpp"


namespace FD
{
	constexpr int CENTRAL = 0;
	constexpr int CENTRAL_ONE_SIDE_FORWARD = 1;
	constexpr int CENTRAL_ONE_SIDE_BACKWARD = 2;


	template<unsigned int dir, unsigned int ord, unsigned int impl>
	struct Derivative_impl
	{
		template<typename rtype,typename expr_type>
		static rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Error derivative order: " << ord << " have not been implemented yet" << std::endl;
		}
	};

	template<unsigned int dir>
	struct Derivative_impl<dir,2,CENTRAL>
	{
		template<typename rtype,typename expr_type>
		static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key)
		{
			rtype ret;

			long int old_val = key.getKeyRef().get(dir);
			key.getKeyRef().set_d(dir, old_val + 1);
			ret = o1.value(key)*1.0/o1.getGrid().spacing(dir)/2.0;
			key.getKeyRef().set_d(dir,old_val);


			old_val = key.getKeyRef().get(dir);
			key.getKeyRef().set_d(dir, old_val - 1);
			ret += o1.value(key)*(-1.0/o1.getGrid().spacing(dir)/2.0 );
			key.getKeyRef().set_d(dir,old_val);

			return ret;
		}
	};

	enum one_side_direction
	{
		OS_CENTRAL,
		OS_FORWARD,
		OS_BACKWARD,
	};

	template<typename gtype, typename key_type>
	one_side_direction use_one_side(gtype & g, unsigned int dir, key_type & key)
	{
		if (g.getDecomposition().periodicity()[dir] == true)
		{
			return one_side_direction::OS_CENTRAL;
		}

		auto keyg = g.getGKey(key);
		if (keyg.get(dir) == 0)
		{
			return one_side_direction::OS_FORWARD;
		}
		else if (keyg.get(dir) == g.getGridInfoVoid().size(dir) - 1)
		{
			return one_side_direction::OS_BACKWARD;
		}

		return one_side_direction::OS_CENTRAL;
	}

	template<unsigned int dir>
	struct Derivative_impl<dir,2,CENTRAL_ONE_SIDE_FORWARD>
	{
		template<typename rtype,typename expr_type>
		static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key)
		{
			rtype ret;

			ret = o1.value(key)*(-3.0/o1.getGrid().spacing(dir)/2.0);

			long int old_val = key.getKeyRef().get(dir);
			key.getKeyRef().set_d(dir, old_val + 1);
			ret += o1.value(key)*2.0/o1.getGrid().spacing(dir);
			key.getKeyRef().set_d(dir,old_val);


			old_val = key.getKeyRef().get(dir);
			key.getKeyRef().set_d(dir, old_val + 2);
			ret += o1.value(key)*(-0.5/o1.getGrid().spacing(dir) );
			key.getKeyRef().set_d(dir,old_val);

			return ret;
		}
	};

	template<unsigned int dir>
	struct Derivative_impl<dir,2,CENTRAL_ONE_SIDE_BACKWARD>
	{
		template<typename rtype,typename expr_type>
		static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key)
		{
			rtype ret;

			ret = o1.value(key)*(3.0/o1.getGrid().spacing(dir)/2.0);

			long int old_val = key.getKeyRef().get(dir);
			key.getKeyRef().set_d(dir, old_val - 1);
			ret += o1.value(key)*(-2.0/o1.getGrid().spacing(dir));
			key.getKeyRef().set_d(dir,old_val);


			old_val = key.getKeyRef().get(dir);
			key.getKeyRef().set_d(dir, old_val - 2);
			ret += o1.value(key)*0.5/o1.getGrid().spacing(dir);
			key.getKeyRef().set_d(dir,old_val);

			return ret;
		}
	};

	template<unsigned dir_, unsigned int ord_, unsigned int impl_>
	struct GRID_DERIVATIVE
	{
		typedef boost::mpl::int_<dir_> dir;
		typedef boost::mpl::int_<ord_> ord;
		typedef boost::mpl::int_<impl_> impl;
	};

	template <typename exp1, unsigned int dir, unsigned int ord, unsigned int impl>
	class grid_dist_expression_op<exp1,void,GRID_DERIVATIVE<dir,ord,impl> >
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
		inline auto value(grid_dist_key_dx<gtype::dims> & key) const -> typename std::remove_reference<decltype(o1.value(key))>::type
		{
			typedef typename std::remove_reference<decltype(o1.value(key))>::type r_type;

			r_type val;

			one_side_direction os = use_one_side(getGrid(),dir,key);

			if (os == one_side_direction::OS_CENTRAL)
			{
				val = Derivative_impl<dir,ord,impl>::template calculate<r_type>(o1,key);
			}
			else if (os == one_side_direction::OS_FORWARD)
			{
				val = Derivative_impl<dir,ord,impl+1>::template calculate<r_type>(o1,key);
			}
			else
			{
				val = Derivative_impl<dir,ord,impl+2>::template calculate<r_type>(o1,key);
			}

			return val;
		}

        template<typename Sys_eqs, typename gmap_type, typename unordered_map_type>
        inline void value_nz(const gmap_type & g_map,
                             grid_dist_key_dx<Sys_eqs::dims> & kmap,
                             const grid_sm<Sys_eqs::dims,void> & gs,
                             typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                             unordered_map_type & cols,
                             typename Sys_eqs::stype coeff,
                             unsigned int comp) const
        {

        /*    if (is_grid_staggered<Sys_eqs>::value())
            {
                D<d,arg,Sys_eqs,BACKWARD>::value(g_map,kmap,gs,spacing,cols,coeff);
                return;
            }*/

            //o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,coeff,comp);

            long int old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,coeff/spacing[dir]/2.0,comp);
            kmap.getKeyRef().set_d(dir,old_val);

            old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-coeff/spacing[dir]/2.0,comp);
            kmap.getKeyRef().set_d(dir,old_val);


           // cols[g_map.template getProp<0>(kmap)*Sys_eqs::nvar + FD::var_id + comp] += coeff;
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
	};


	template<unsigned int dir, unsigned int ord, unsigned int impl>
	class Derivative
	{

	public:

		Derivative()
		{
		}

		template<typename operand_type>
		grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<dir,ord,impl>> operator()(operand_type arg)
		{
			return grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<dir,ord,impl>>(arg);
		}
	};

	typedef Derivative<0,2,CENTRAL> Derivative_x;
    typedef Derivative<1,2,CENTRAL> Derivative_y;

};


#endif
