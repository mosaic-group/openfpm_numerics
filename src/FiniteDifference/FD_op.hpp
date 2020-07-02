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


    // ord_d means order of the derivative and not order of convergence
    template<unsigned int dir,unsigned int ord_d,unsigned int ord, unsigned int impl>
    struct Derivative_impl
    {
        template<typename rtype,typename expr_type>
        static rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key)
        {
            std::cerr << __FILE__ << ":" << __LINE__ << " Error the give derivative scheme with convergence order " << ord << " has not been implemented yet" << std::endl;
        }
    };


    template<unsigned int dir>
    struct Derivative_impl<dir,1,2,CENTRAL>
    {
        template<typename rtype,typename expr_type>
        static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key,comb<expr_type::gtype::dims> & c_where)
        {
            // x0, dx are defined in proper dir є(x, y, z)
            auto dx = o1.getGrid().spacing(dir);
            long int x0 = key.getKeyRef().get(dir);

            key.getKeyRef().set_d(dir, x0 + 1);
            rtype ret = o1.value(key,c_where)*1.0/2.0;

            key.getKeyRef().set_d(dir, x0 - 1);
            ret += o1.value(key,c_where)*(-1.0/2.0);

            key.getKeyRef().set_d(dir, x0);
            return ret/dx;
        }

        template<typename Sys_eqs, typename o1_type, typename gmap_type, typename unordered_map_type>
        inline static void calculate_nz(o1_type o1,
                                 const gmap_type & g_map,
                                 grid_dist_key_dx<Sys_eqs::dims> & kmap,
                                 const grid_sm<Sys_eqs::dims,void> & gs,
                                 typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                                 unordered_map_type & cols,
                                 typename Sys_eqs::stype coeff,
                                 unsigned int comp)
        {
            long int old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,coeff/spacing[dir]/2.0,comp);
            kmap.getKeyRef().set_d(dir,old_val);

            old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-coeff/spacing[dir]/2.0,comp);
            kmap.getKeyRef().set_d(dir,old_val);
        }
    };

    template<unsigned int dir>
    struct Derivative_impl<dir,2,2,CENTRAL>
    {
        template<typename rtype,typename expr_type>
        static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key, comb<expr_type::gtype::dims> & c_where)
        {
            // x0, dx are defined in proper dir є(x, y, z)
            auto dx = o1.getGrid().spacing(dir);
            long int x0 = key.getKeyRef().get(dir);

            rtype ret = o1.value(key,c_where)*(-2.0);

            key.getKeyRef().set_d(dir, x0 + 1);
            ret += o1.value(key,c_where)*(1.0);

            key.getKeyRef().set_d(dir, x0 - 1);
            ret += o1.value(key,c_where)*(1.0);

            key.getKeyRef().set_d(dir, x0);
            return ret/(dx*dx);
        }

        template<typename Sys_eqs, typename o1_type, typename gmap_type, typename unordered_map_type>
        inline static void calculate_nz(o1_type o1,
                                        const gmap_type & g_map,
                                        grid_dist_key_dx<Sys_eqs::dims> & kmap,
                                        const grid_sm<Sys_eqs::dims,void> & gs,
                                        typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                                        unordered_map_type & cols,
                                        typename Sys_eqs::stype coeff,
                                        unsigned int comp)
        {   long int old_val = kmap.getKeyRef().get(dir);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-2.0*coeff/(spacing[dir]*spacing[dir]),comp);

            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);

            old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
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
    struct Derivative_impl<dir,1,2,CENTRAL_ONE_SIDE_FORWARD>
    {
        template<typename rtype,typename expr_type>
        static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key, comb<expr_type::gtype::dims> & c_where)
        {
            // x0, dx are defined in proper dir є(x, y, z)
            auto dx = o1.getGrid().spacing(dir);
            long int x0 = key.getKeyRef().get(dir);

            rtype ret = o1.value(key,c_where)*(-3.0/2.0);

            key.getKeyRef().set_d(dir, x0 + 1);
            ret += o1.value(key,c_where)*2.0;

            key.getKeyRef().set_d(dir, x0 + 2);
            ret += o1.value(key,c_where)*(-0.5);

            key.getKeyRef().set_d(dir, x0);
            return ret/dx;
        }

        template<typename Sys_eqs, typename o1_type, typename gmap_type, typename unordered_map_type>
        inline static void calculate_nz(o1_type o1,
                                        const gmap_type & g_map,
                                        grid_dist_key_dx<Sys_eqs::dims> & kmap,
                                        const grid_sm<Sys_eqs::dims,void> & gs,
                                        typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                                        unordered_map_type & cols,
                                        typename Sys_eqs::stype coeff,
                                        unsigned int comp)
        {   long int old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir));
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-3.0*coeff/spacing[dir]/2.0,comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,2.0*coeff/spacing[dir],comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 2);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-0.5*coeff/spacing[dir],comp);
            kmap.getKeyRef().set_d(dir,old_val);
        }
    };

    template<unsigned int dir>
    struct Derivative_impl<dir,2,2,CENTRAL_ONE_SIDE_FORWARD>
    {
        template<typename rtype,typename expr_type>
        static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key, comb<expr_type::gtype::dims> & c_where)
        {
            // x0, dx are defined in proper dir є(x, y, z)
            auto dx = o1.getGrid().spacing(dir);
            long int x0 = key.getKeyRef().get(dir);

            rtype ret = o1.value(key,c_where)*(2.0);

            key.getKeyRef().set_d(dir, x0 + 1);
            ret += o1.value(key,c_where)*(-5.0);

            key.getKeyRef().set_d(dir, x0 + 2);
            ret += o1.value(key,c_where)*(4.0);

            key.getKeyRef().set_d(dir, x0 + 3);
            ret += o1.value(key,c_where)*(-1.0);

            key.getKeyRef().set_d(dir, x0);
            return ret/(dx*dx);
        }

        template<typename Sys_eqs, typename o1_type, typename gmap_type, typename unordered_map_type>
        inline static void calculate_nz(o1_type o1,
                                        const gmap_type & g_map,
                                        grid_dist_key_dx<Sys_eqs::dims> & kmap,
                                        const grid_sm<Sys_eqs::dims,void> & gs,
                                        typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                                        unordered_map_type & cols,
                                        typename Sys_eqs::stype coeff,
                                        unsigned int comp)
        {
            long int old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir));
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,2.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-5.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 2);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,4.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) + 3);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-1.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
        }
    };

    template<unsigned int dir>
    struct Derivative_impl<dir,1,2,CENTRAL_ONE_SIDE_BACKWARD>
    {
        template<typename rtype,typename expr_type>
        static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key, comb<expr_type::gtype::dims> & c_where)
        {
            // x0, dx are defined in proper dir є(x, y, z)
            auto dx = o1.getGrid().spacing(dir);
            long int x0 = key.getKeyRef().get(dir);

            rtype ret = o1.value(key,c_where)*(3.0/2.0);

            key.getKeyRef().set_d(dir, x0 - 1);
            ret += o1.value(key,c_where)*(-2.0);

            key.getKeyRef().set_d(dir, x0 - 2);
            ret += o1.value(key,c_where)*0.5;

            key.getKeyRef().set_d(dir, x0);
            return ret/dx;
        }

        template<typename Sys_eqs, typename o1_type, typename gmap_type, typename unordered_map_type>
        inline static void calculate_nz(o1_type o1,
                                        const gmap_type & g_map,
                                        grid_dist_key_dx<Sys_eqs::dims> & kmap,
                                        const grid_sm<Sys_eqs::dims,void> & gs,
                                        typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                                        unordered_map_type & cols,
                                        typename Sys_eqs::stype coeff,
                                        unsigned int comp)
        {   long int old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir));
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,3.0*coeff/spacing[dir]/2.0,comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-2.0*coeff/spacing[dir],comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 2);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,0.5*coeff/spacing[dir],comp);
            kmap.getKeyRef().set_d(dir,old_val);
        }
    };

    template<unsigned int dir>
    struct Derivative_impl<dir,2,2,CENTRAL_ONE_SIDE_BACKWARD>
    {
        template<typename rtype,typename expr_type>
        static inline rtype calculate(expr_type & o1, grid_dist_key_dx<expr_type::gtype::dims> & key, comb<expr_type::gtype::dims> & c_where)
        {
            // x0, dx are defined in proper dir є(x, y, z)
            auto dx = o1.getGrid().spacing(dir);
            long int x0 = key.getKeyRef().get(dir);

            rtype ret = o1.value(key,c_where)*(2.0);

            key.getKeyRef().set_d(dir, x0 - 1);
            ret += o1.value(key,c_where)*(-5.0);

            key.getKeyRef().set_d(dir, x0 - 2);
            ret += o1.value(key,c_where)*(4.0);

            key.getKeyRef().set_d(dir, x0 - 3);
            ret += o1.value(key,c_where)*(-1.0);

            key.getKeyRef().set_d(dir, x0);
            return ret/(dx*dx);
        }

        template<typename Sys_eqs, typename o1_type, typename gmap_type, typename unordered_map_type>
        inline static void calculate_nz(o1_type o1,
                                        const gmap_type & g_map,
                                        grid_dist_key_dx<Sys_eqs::dims> & kmap,
                                        const grid_sm<Sys_eqs::dims,void> & gs,
                                        typename Sys_eqs::stype (& spacing )[Sys_eqs::dims],
                                        unordered_map_type & cols,
                                        typename Sys_eqs::stype coeff,
                                        unsigned int comp)
        {               long int old_val = kmap.getKeyRef().get(dir);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir));
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,2.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 1);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-5.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 2);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,4.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
            kmap.getKeyRef().set_d(dir, kmap.getKeyRef().get(dir) - 3);
            o1.template value_nz<Sys_eqs>(g_map,kmap,gs,spacing,cols,-1.0*coeff/(spacing[dir]*spacing[dir]),comp);
            kmap.getKeyRef().set_d(dir,old_val);
        }
    };

    template<unsigned dir_, unsigned int ord_d_, unsigned int ord_, unsigned int impl_>
    struct GRID_DERIVATIVE
    {
        typedef boost::mpl::int_<dir_> dir;
        typedef boost::mpl::int_<ord_> ord;
        typedef boost::mpl::int_<ord_d_> ord_d;
        typedef boost::mpl::int_<impl_> impl;
    };

    template <typename exp1, unsigned int dir, unsigned int  ord_d,unsigned int ord, unsigned int impl>
    class grid_dist_expression_op<exp1,void,GRID_DERIVATIVE<dir,ord_d,ord,impl> >
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
            typedef typename std::remove_reference<decltype(o1.value(key,c_where))>::type r_type;

            r_type val;

            one_side_direction os = use_one_side(getGrid(),dir,key);

            if (os == one_side_direction::OS_CENTRAL)
            {
                val = Derivative_impl<dir,ord_d,ord,impl>::template calculate<r_type>(o1,key,c_where);
            }
            else if (os == one_side_direction::OS_FORWARD)
            {
                val = Derivative_impl<dir,ord_d,ord,impl+1>::template calculate<r_type>(o1,key,c_where);
            }
            else
            {
                val = Derivative_impl<dir,ord_d,ord,impl+2>::template calculate<r_type>(o1,key,c_where);
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
            one_side_direction os = use_one_side(getGrid(),dir,kmap);
            if (os == one_side_direction::OS_CENTRAL)
            {
                Derivative_impl<dir,ord_d,ord,impl>::template calculate_nz<Sys_eqs>(o1,g_map,kmap,gs,spacing,cols,coeff,comp);
            }
            else if (os == one_side_direction::OS_FORWARD)
            {
                Derivative_impl<dir,ord_d,ord,impl+1>::template calculate_nz<Sys_eqs>(o1,g_map,kmap,gs,spacing,cols,coeff,comp);
            }
            else
            {
                Derivative_impl<dir,ord_d,ord,impl+2>::template calculate_nz<Sys_eqs>(o1,g_map,kmap,gs,spacing,cols,coeff,comp);
            }

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




    template<unsigned int dir,unsigned int ord_d ,unsigned int ord, unsigned int impl>
    class Derivative
    {

    public:

        Derivative()
        {
        }

        template<typename operand_type>
        grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<dir,ord_d,ord,impl>> operator()(operand_type arg)
        {
            return grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<dir,ord_d,ord,impl>>(arg);
        }
    };

    template<unsigned int dim, unsigned int impl>
    class Laplacian
    {

    public:

        Laplacian()
        {
        }

        template<typename operand_type>
        FD::grid_dist_expression_op<
                FD::grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<0,2,2,impl>>,
                FD::grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<1,2,2,impl>>,
                FD::sum>
        operator()(operand_type arg)
        {
            grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<0,2,2,impl>> d_xx(arg);
            grid_dist_expression_op<operand_type,void,GRID_DERIVATIVE<1,2,2,impl>> d_yy(arg);
            return d_xx + d_yy;
        }
    };

	class L2Error
	{
	public:
		L2Error() {}

		template<unsigned int p1, unsigned int p2, typename g1, typename g2, unsigned int impl_p1, unsigned int impl_p2>
		auto operator()(const FD::grid_dist_expression<p1,g1,impl_p1> ga, const FD::grid_dist_expression<p2,g2,impl_p2> gb) const ->
		typename std::remove_reference<decltype(ga.getGrid().template getProp<p1>(ga.getGrid().getDomainIterator().get()))>::type
		{
			typename std::remove_reference<decltype(ga.getGrid().template getProp<p1>(ga.getGrid().getDomainIterator().get()))>::type result = 0.0;
			comb<FD::grid_dist_expression<p1,g1,impl_p1>::gtype::dims> s_pos;
			auto diff_expr = ga + -1*gb;

			auto & v_cl = create_vcluster();

			auto it = ga.getGrid().getDomainIterator();
			while (it.isNext())
			{
				auto key = it.get();
				auto diff = diff_expr.value(key,s_pos);
				result += diff*diff;
				++it;
			}

			v_cl.sum(result);
			v_cl.execute();

			return result;
		}
	};


    typedef Derivative<0,1,2,CENTRAL> Derivative_x;
    typedef Derivative<1,1,2,CENTRAL> Derivative_y;
    typedef Derivative<0,2,2,CENTRAL> Derivative_xx;
    typedef Derivative<1,2,2,CENTRAL> Derivative_yy;
    typedef Laplacian<2,CENTRAL> Laplacian_xy;
};


#endif
