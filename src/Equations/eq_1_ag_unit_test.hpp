/*
 * eq_1_ag_unit_test.hpp
 *
 * This file model the Eq1 in 2D of Active Gel
 *
 * R. Ramaswamy, G. Bourantas, F. Jülicher, and I. F. Sbalzarini.
 * A hybrid particle-mesh method for incompressible active polar
 * viscous gels. J. Comput. Phys., 291:334–361, 2015.
 *
 *  Created on: Sep 18, 2015
 *      Author: Pietro Incardona
 */

#ifndef OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_1_AG_UNIT_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_1_AG_UNIT_TEST_HPP_


////////////////////////////////////// Left part  eq 1 //////////////////////////////

// Creating square parentesis term left part equation 1

// Derivation and multiplication of t1 ... t6

namespace eq1
{
	typedef mul<etad2,D<x,t1>> t1_l;
	typedef mul<eta,D<x,t2>> t2_l;
	typedef mul<etad2,D<x,t3>> t3_l;

	typedef mul<eta,D<y,t4>> t4_l;
	typedef mul<eta2,D<y,t5>> t5_l;
	typedef mul<eta,D<y,t6>> t6_l;

	typedef mul<eta,Lap<v>> meLv;
	typedef D<x,P> DxP;

	// Assemble the left part
	typedef add<meLV,minus<DxP>,t1_l,t2_l,t3_l,t4_l,t5_l,t6_l> left;

	////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////// Right part  eq 1 /////////////////////////////

	typedef div<2,D<x,hp>> t1_r;
	typedef mul<zeta,D<x,mul<del_nu,px2>>> t2_r;
	typedef mul<zeta,D<y,mul<del_nu,px_py>>> t3_r;
	typedef mul<zeta,D<x,mul<del_nu,dev<2,px2_p_py2>>>> t4_r;
	typedef mul<nud2,D<x,mul<-2,hp,px_py>> t5_r;
	typedef mul<etad2,D<y,mul<hp,px2_m_py2>>> t6_r;
	typedef D<x,sigma[x][x]> t7_r;
	typedef D<y,sigma[x][y]> t8_r;

	typedef add<t1_r,t2_r,t3_r,t4_r,t5_r,t6_r,t7_r,t8_r,g_ext> eq1r;

	typedef Eq<eq1l,eq1r> equation;
}

/////////////////////////////////////////////////////////////////////////////////////

#endif /* OPENFPM_NUMERICS_SRC_EQUATIONS_EQ_1_AG_UNIT_TEST_HPP_ */
