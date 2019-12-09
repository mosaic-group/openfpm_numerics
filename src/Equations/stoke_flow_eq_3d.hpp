/*
 * stoke_flow_eq_3d.hpp
 *
 *  Created on: May 28, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_EQUATIONS_STOKE_FLOW_EQ_3D_HPP_
#define OPENFPM_NUMERICS_SRC_EQUATIONS_STOKE_FLOW_EQ_3D_HPP_

constexpr unsigned int x = 0;
constexpr unsigned int y = 1;
constexpr unsigned int z = 2;

// Model the equations

constexpr unsigned int v[] = {0,1,2};
constexpr unsigned int P = 3;

typedef Field<v[x],lid_nn_3d> v_x;
typedef Field<v[y],lid_nn_3d> v_y;
typedef Field<v[z],lid_nn_3d> v_z;
typedef Field<P,lid_nn_3d> Prs;

// Eq1 V_x

typedef mul<eta<lid_nn_3d>,Lap<v_x>> eta_lap_vx;
typedef D<x,Prs> p_x;
typedef minus<p_x> m_p_x;
typedef sum<eta_lap_vx,m_p_x> vx_eq;

// Eq2 V_y

typedef mul<eta<lid_nn_3d>,Lap<v_y>> eta_lap_vy;
typedef D<y,Prs> p_y;
typedef minus<p_y> m_p_y;
typedef sum<eta_lap_vy,m_p_y> vy_eq;

// Eq3 V_z

typedef mul<eta<lid_nn_3d>,Lap<v_z>> eta_lap_vz;
typedef D<z,Prs> p_z;
typedef minus<p_z> m_p_z;
typedef sum<eta_lap_vz,m_p_z> vz_eq;

// Eq4 Incompressibility

typedef D<x,v_x,FORWARD> dx_vx;
typedef D<y,v_y,FORWARD> dy_vy;
typedef D<z,v_z,FORWARD> dz_vz;
typedef sum<sum<dx_vx,dy_vy>,dz_vz> ic_eq;


// Directional Avg
typedef Avg<x,v_y> avg_x_vy;
typedef Avg<z,v_y> avg_z_vy;

typedef Avg<y,v_x> avg_y_vx;
typedef Avg<z,v_x> avg_z_vx;

typedef Avg<y,v_z> avg_y_vz;
typedef Avg<x,v_z> avg_x_vz;

// Directional Avg

typedef Avg<x,v_y,FORWARD> avg_x_vy_f;
typedef Avg<z,v_y,FORWARD> avg_z_vy_f;

typedef Avg<y,v_x,FORWARD> avg_y_vx_f;
typedef Avg<z,v_x,FORWARD> avg_z_vx_f;

typedef Avg<y,v_z,FORWARD> avg_y_vz_f;
typedef Avg<x,v_z,FORWARD> avg_x_vz_f;

#define EQ_1 0
#define EQ_2 1
#define EQ_3 2
#define EQ_4 3


#endif /* OPENFPM_NUMERICS_SRC_EQUATIONS_STOKE_FLOW_EQ_3D_HPP_ */
