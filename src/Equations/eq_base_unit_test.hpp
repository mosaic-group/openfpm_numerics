/*
 * eq_base_unit_test.hpp
 *
 *  Created on: Sep 18, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_BASE_UNIT_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_BASE_UNIT_TEST_HPP_

// quadratic and cubic term in p


typedef mul<p[x],p[x]> px2;
typedef mul<p[y],p[y]> py2;
typedef mul<p[x],p[y]> px_py;
typedef add<px2,py2> px2_p_py2;
typedef add<px2,py2> px2_s_py2;
typedef mul<px2,px> px3;
typedef mul<py2,py> py3;
typedef mul<gm,nu> gm_nu;

// quadratic term in p

typedef mul<px3,py> px3_m_py;
typedef mul<px2,py2> px2_m_py2;
typedef mul<px,py3> px_m_py3;

typedef mul<px2,px2_s_py2> px2_r;
typedef mul<px_py,px2_s_py2> px_py_r;
typedef mul<py2,px2_s_py2> py2_r;

// Division of the quadratic term by px2_s_py2

typedef div<px3_m_py,px2_p_py2> p_c1;
typedef div<px2_m_py2,px2_p_py2> p_c2;
typedef div<px_m_py3,px2_p_py2> p_c3;

typedef div<px2_r,px2_p_py2> p_c4;
typedef div<px_py_r,px2_p_py2> p_c5;
typedef div<py2_r,px2_p_py2> p_c6;

// Create u_alpha_beta equation

template <unsigned int i, unsigned int j> using u = div<2,add<D<i,v[j]>,D<j,v[i]>>>;

// terms

typedef mul<gm,p_c1,u[x][x]> t1;
typedef mul<gm,p_c2,u[x][y]> t2;
typedef mul<gm,p_c3,u[y][y]> t3;

typedef mul<gm,p_c4,u[x][x]> t4;
typedef mul<gm,p_c5,u[x][y]> t5;
typedef mul<gm,p_c6,u[y][y]> t6;

#endif /* OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_BASE_UNIT_TEST_HPP_ */
