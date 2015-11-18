/*
 * eq_2_ag_unit_test.hpp
 *
 *  Created on: Sep 18, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_2_AG_UNIT_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_2_AG_UNIT_TEST_HPP_

// Left part equation 2

namespace eq2
{
	typedef mul<etad2,D<y,t1>> t1_l;
	typedef mul<eta,D<y,t2>> t2_l;
	typedef mul<etad2,D<y,t3>> t3_l;

	typedef mul<eta,D<x,t4>> t4_l;
	typedef mul<eta2,D<x,t5>> t5_l;
	typedef mul<eta,D<x,t6>> t6_l;

	typedef mul<eta,Lap<v>> meLv;
	typedef D<y,P> DyP;

	// left part equation 2

	typedef add<meLv,DyP,Dt1,Dt2,Dt3,Dt4,Dt6> left;

	typedef div<-2,D<x,hp>> t1_r;
	typedef mul<zeta,D<y,mul<del_nu,py2>>> t2_r;
	typedef mul<zeta,D<x,mul<del_nu,px_py>>> t3_r;
	typedef mul<zeta,D<y,mul<del_nu,dev<2,px2_p_py2>>>> t4_r;
	typedef mul<zeta,D<y,mul<del_nu,dev<2,px2_p_py2>>>> t5_r;
	typedef mul<nud2,D<x,mul<hp,px2_m_py2>>> t6_r;
	typedef D<x,sigma[x][y]> t7_r;
	typedef D<y,sigma[y][x]> t8_r;
	typedef mul<nud2,D<y,mul<gm,lambda,del_nu,px2_m_py2>> t9_r

	typedef add<t1_r,t2_r,t3_r,t4_r,t5_r,t6_r,t7_r,t8_r,t9_r> right;
}



#endif /* OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_2_AG_UNIT_TEST_HPP_ */
