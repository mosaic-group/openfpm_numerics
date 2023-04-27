//
// Created by foggia on 19th Jan 2022
// It's a modification of Abhinav's test, adapted for grids
//

#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <boost/test/unit_test.hpp>

#include "config.h"
#include "Grid/grid_dist_id.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "util/util_debug.hpp"
#include "util/common.hpp"

const double a = 2.8e-4;
const double b = 5e-3;
const double tau = .1;
const double k = .005;
const int dim = 2;

void *gridGlobal;
typedef grid_dist_id<2,double,aggregate<double,double,double,double,double,double>> grid_type;

// State types for systems with different number of ODEs
typedef state_type_ofpm_impl<1,texp_g<dim,double>> state_type_1ode;
typedef state_type_ofpm_impl<2,texp_g<dim,double>> state_type_2ode;
typedef state_type_ofpm_impl<3,texp_g<dim,double>> state_type_3ode;

template<typename DX, typename DY>
struct Fitz {

  DX & ddx;
  DY & ddy;
  
  //Constructor
  Fitz(DX & m_ddx, DY & m_ddy) :
    ddx(m_ddx),
    ddy(m_ddy) {}
  
  void operator()(const state_type_2ode & x,
		  state_type_2ode & dxdt,
		  const double t) const {

    grid_type & temp = *(grid_type *) gridGlobal;
    auto u{FD::getV<4>(temp)};
    auto v{FD::getV<5>(temp)};
    u = x.data.get<0>();
    v = x.data.get<1>();
    
    temp.ghost_get<4,5>();
    dxdt.data.get<0>() = ddx(u) + ddy(u) + (1.0);
    dxdt.data.get<1>() = ddx(v) + ddy(v) + (2.0);

    // One point stay fixed
    
    auto key2 = temp.getDomainIterator().get();

    if (create_vcluster().rank() == 0) {
      dxdt.data.get<0>().value_ref(key2) = 0.0;
      dxdt.data.get<1>().value_ref(key2) = 0.0;
    }

    double u_max{0.0};
    double v_max{0.0};

    auto it = temp.getDomainIterator();

    while (it.isNext())
    {
      auto key = it.get();

      if (u_max < dxdt.data.get<0>().value(key))
        u_max = dxdt.data.get<0>().value(key);
          
      if (v_max < dxdt.data.get<1>().value(key))
        v_max = dxdt.data.get<1>().value(key);

      ++it;
    }
  }
};

void Exponential_struct_ofp2(const state_type_3ode & x,
			     state_type_3ode & dxdt,
			     const double t) {

  // sytem: dx1/dt = x1 --> solution: x1(t) = exp(t)
  // sytem: dx2/dt = 2*x2 --> solution: x2(t) = exp(2t)
  dxdt.data.get<0>() = x.data.get<0>();
  dxdt.data.get<1>() = 2.0 * x.data.get<1>();
  dxdt.data.get<2>() = x.data.get<0>();
}


void Exponential(const state_type_1ode & x,
		 state_type_1ode & dxdt,
		 const double t) {

  // sytem: dx/dt = x --> solution: x(t) = exp(t)
  dxdt = x;
}

// void sigmoid(const state_type_1ode & x,
// 	     state_type_1ode & dxdt,
// 	     const double t) {
//   dxdt = x * (1.0 - x);
// }

BOOST_AUTO_TEST_SUITE(odeInt_grid_tests)

BOOST_AUTO_TEST_CASE(odeint_grid_test_exponential) {

  size_t edgeSemiSize{40};
  const size_t sz[dim] = {edgeSemiSize,edgeSemiSize};
  Box<dim,double> box{{0.0,0.0}, {1.0,1.0}};
  periodicity<dim> bc{{NON_PERIODIC,NON_PERIODIC}};
  double spacing[dim];
  spacing[0] = 1.0 / (sz[0] - 1);
  spacing[1] = 1.0 / (sz[1] - 1);
  Ghost<dim,long int> ghost{2};
  BOOST_TEST_MESSAGE("Test: exponential");
  BOOST_TEST_MESSAGE("Init grid_dist_id ...");

  grid_dist_id<dim,double,aggregate<double,double,double>> grid{sz,box,ghost,bc};

  auto it{grid.getDomainIterator()};
  while (it.isNext()) {
    auto key = it.get();
    grid.template get<0>(key) = std::exp(0);   // Initial state
    grid.template get<1>(key) = std::exp(0.4); // Analytical solution
    ++it;
  }
  grid.ghost_get<0>();

  auto Init{FD::getV<0>(grid)};   // Initial state
  auto Sol{FD::getV<1>(grid)};    // Analytical solution
  auto OdeSol{FD::getV<2>(grid)}; // Numerical solution

  state_type_1ode x0;
  x0.data.get<0>() = Init;

  double t{0.0};
  double tf{0.4};
  const double dt{0.1};
  
  boost::numeric::odeint::runge_kutta4<state_type_1ode,double,
  				       state_type_1ode,double,
  				       boost::numeric::odeint::vector_space_algebra_ofp> rk4; // Time integrator
  size_t steps{boost::numeric::odeint::integrate_const(rk4,Exponential,x0,0.0,tf,dt)};

  OdeSol = x0.data.get<0>(); // Numerical solution

  // Error
  auto it2{grid.getDomainIterator()};
  double worst{0.0};
  while (it2.isNext()) {
    auto p{it2.get()};
    if (std::fabs(grid.template get<1>(p) - grid.template get<2>(p)) > worst) {
      worst = std::fabs(grid.template get<1>(p) - grid.template get<2>(p));
    }
    ++it2;
  }

  std::cout << worst << std::endl;
  BOOST_REQUIRE(worst < 1e-6);

  // Another way
  x0.data.get<0>() = Init;
  while (t < tf) {
    rk4.do_step(Exponential,x0,t,dt);
    OdeSol = x0.data.get<0>();
    t += dt;
  }

  OdeSol = x0.data.get<0>();

  // Error
  auto it3{grid.getDomainIterator()};
  double worst2{0.0};
  while (it3.isNext()) {
    auto p{it3.get()};
    if (std::fabs(grid.template get<1>(p) - grid.template get<2>(p)) > worst2) {
      worst2 = fabs(grid.template get<1>(p) - grid.template get<2>(p));
    }
    ++it3;
  }
  std::cout << worst2 << std::endl;
  BOOST_REQUIRE(worst2 < 1e-6);
  BOOST_REQUIRE_EQUAL(worst,worst2);
}



BOOST_AUTO_TEST_CASE(odeint_grid_test_STRUCT_exponential) {

  size_t edgeSemiSize{40};
  const size_t sz[dim] = {edgeSemiSize,edgeSemiSize};
  Box<dim, double> box{{0.0,0.0},{1.0,1.0}};
  periodicity<dim> bc{{NON_PERIODIC,NON_PERIODIC}};
  double spacing[dim];
  spacing[0] = 1.0 / (sz[0] - 1);
  spacing[1] = 1.0 / (sz[1] - 1);
  Ghost<dim,long int> ghost{2};
  BOOST_TEST_MESSAGE("Test: exponential");
  BOOST_TEST_MESSAGE("Init grid_dist_id ...");

  grid_dist_id<dim,double,aggregate<double,double,double,double,double,double>> grid{sz,box,ghost,bc};

  auto it{grid.getDomainIterator()};
  while (it.isNext()) {
    auto key = it.get();
    
    grid.get<0>(key) = std::exp(0);
    grid.template get<0>(key) = std::exp(0.0); // Initial state 1
    grid.template get<1>(key) = std::exp(0.4); // Analytical solution 1
    grid.template get<2>(key) = std::exp(0.0); // Initial state 2
    grid.template get<3>(key) = std::exp(0.8); // Analytical solution 2
    ++it;
  }
  grid.ghost_get<0>();

  auto Init1{FD::getV<0>(grid)};   // Initial state 1
  auto Sol1{FD::getV<1>(grid)};    // Analytical solution 1
  
  auto Init2{FD::getV<2>(grid)};   // Initial state 2
  auto Sol2{FD::getV<3>(grid)};    // Analytical solution 2
  
  auto OdeSol1{FD::getV<4>(grid)}; // Numerical solution 1
  auto OdeSol2{FD::getV<5>(grid)}; // Numerical solution 2

  state_type_3ode x0;
  x0.data.get<0>() = Init1;
  x0.data.get<1>() = Init2;
  x0.data.get<2>() = Init1;

  double t{0};
  double tf{0.4};
  const double dt{0.1};

  // size_t steps{boost::numeric::odeint::integrate(Exponential_struct,x0,0.0,tf,dt)};
  // size_t steps{boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta4<state_type_3ode,double,state_type_3ode,double,boost::numeric::odeint::vector_space_algebra_ofp>(),
  // 						       Exponential_struct_ofp,x0,0.0,tf,dt)};
  
  typedef boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<state_type_3ode,double,state_type_3ode,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;
  integrate_adaptive(stepper_type(),Exponential_struct_ofp2,x0,t,tf,dt);
  
  OdeSol1 = x0.data.get<0>();
  OdeSol2 = x0.data.get<1>();

  // Error
  auto it2{grid.getDomainIterator()};
  double worst{0.0};
  double worst2{0.0};
  while (it2.isNext()) {
    auto p{it2.get()};
    if (std::fabs(grid.getProp<1>(p) - grid.getProp<4>(p)) > worst)
      worst = std::fabs(grid.getProp<1>(p) - grid.getProp<4>(p));
    if (std::fabs(grid.getProp<3>(p) - grid.getProp<5>(p)) > worst2)
      worst2 = std::fabs(grid.getProp<3>(p) - grid.getProp<5>(p));
    ++it2;
  }

  std::cout << worst << " " << worst2 << std::endl;
  BOOST_REQUIRE(worst < 1e-6);
  BOOST_REQUIRE(worst2 < 1e-6);


  // A different way
  x0.data.get<0>() = Init1;
  x0.data.get<1>() = Init2;
  x0.data.get<2>() = Init1;
  
  boost::numeric::odeint::runge_kutta4<state_type_3ode,double,state_type_3ode,double,boost::numeric::odeint::vector_space_algebra_ofp> rk4;
  while (t < tf) {
    rk4.do_step(Exponential_struct_ofp2,x0,t,dt);
    t+=dt;
  }

  OdeSol1 = x0.data.get<0>();
  OdeSol2 = x0.data.get<1>();

  // Error
  auto it3{grid.getDomainIterator()};
  double worst3{0.0};
  double worst4{0.0};
  while (it3.isNext()) {
    auto p{it3.get()};
    if (std::fabs(grid.getProp<1>(p) - grid.getProp<4>(p)) > worst3)
      worst3 = std::fabs(grid.getProp<1>(p) - grid.getProp<4>(p));
    if (std::fabs(grid.getProp<3>(p) - grid.getProp<5>(p)) > worst4)
      worst4 = std::fabs(grid.getProp<3>(p) - grid.getProp<5>(p));
    ++it3;
  }

  std::cout << worst3 << " " << worst4 << std::endl;
  BOOST_REQUIRE(worst3 < 1e-6);
  BOOST_REQUIRE(worst4 < 5e-5);
}




BOOST_AUTO_TEST_CASE(odeint_grid_test2_exponential) {

  size_t edgeSemiSize{40};
  const size_t sz[dim] = {edgeSemiSize,edgeSemiSize};
  Box<dim,double> box{{0.0,0.0}, {1.0,1.0}};
  periodicity<dim> bc{{NON_PERIODIC,NON_PERIODIC}};
  double spacing[dim];
  spacing[0] = 1.0 / (sz[0] - 1);
  spacing[1] = 1.0 / (sz[1] - 1);
  Ghost<dim,long int> ghost{2};
  BOOST_TEST_MESSAGE("Test: exponential");
  BOOST_TEST_MESSAGE("Init grid_dist_id ...");

  grid_dist_id<dim,double,aggregate<double,double,double>> grid{sz,box,ghost,bc};

  double t{0.0};
  double tf{0.5};
  const double dt{0.1};

  auto it{grid.getDomainIterator()};
  while (it.isNext()) {
    auto key = it.get();
    grid.template get<0>(key) = std::exp(t);  // Initial state
    grid.template get<1>(key) = std::exp(tf); // Analytical solution
    ++it;
  }
  grid.ghost_get<0>();

  auto Init{FD::getV<0>(grid)};   // Initial state
  auto Sol{FD::getV<1>(grid)};    // Analytical solution
  auto OdeSol{FD::getV<2>(grid)}; // Numerical solution

  state_type_1ode x0;
  x0.data.get<0>() = Init;

  typedef boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54<state_type_1ode,double,state_type_1ode,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;

  
  integrate_adaptive(stepper_type(),Exponential,x0,t,tf,dt);
  OdeSol = x0.data.get<0>(); // Numerical solution
  
  // Error
  auto it2{grid.getDomainIterator()};
  double worst{0.0};
  while (it2.isNext()) {
    auto p{it2.get()};
    if (std::fabs(grid.template get<1>(p) - grid.template get<2>(p)) > worst) {
      worst = std::fabs(grid.template get<1>(p) - grid.template get<2>(p));
    }
    ++it2;
  }
  
  std::cout << worst << std::endl;
  BOOST_REQUIRE(worst < 1e-6);
  
  // Another way
  boost::numeric::odeint::runge_kutta4<state_type_1ode,double,state_type_1ode,double,boost::numeric::odeint::vector_space_algebra_ofp> rk4;
  x0.data.get<0>() = Init;
  for (size_t i = 0; i < int(tf/dt); ++i, t += dt) {
    rk4.do_step(Exponential,x0,t,dt);
    t += dt;
  }
  OdeSol = x0.data.get<0>();

  // Error
  auto it3{grid.getDomainIterator()};
  double worst2{0.0};
  while (it3.isNext()) {
    auto p{it3.get()};
    if (std::fabs(grid.template get<1>(p) - grid.template get<2>(p)) > worst2) {
      worst2 = fabs(grid.template get<1>(p) - grid.template get<2>(p));
    }
    ++it3;
  }
  std::cout << worst2 << std::endl;
  BOOST_REQUIRE(worst2 < 1e-6);

  // Yet another way
  // x0.data.get<0>() = Init;
  // integrate(rk4,Exponential,x0,t,tf,dt);

  // OdeSol = x0.data.get<0>();

  // // Error
  // auto it4{grid.getDomainIterator()};
  // double worst3{0.0};
  // while (it4.isNext()) {
  //   auto p{it4.get()};
  //   if (std::fabs(grid.template get<1>(p) - grid.template get<2>(p)) > worst3) {
  //     worst3 = fabs(grid.template get<1>(p) - grid.template get<2>(p));
  //   }
  //   ++it4;
  // }
  // std::cout << worst3 << std::endl;
  // BOOST_REQUIRE(worst3 < 1e-6);
  
  // BOOST_REQUIRE_EQUAL(worst,worst2);
  // BOOST_REQUIRE_EQUAL(worst2,worst3);
}



// BOOST_AUTO_TEST_CASE(odeint_base_test3) 
// {
//     size_t edgeSemiSize = 40;
//     const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
//     Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
//     size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
//     double spacing[2];
//     spacing[0] = 1.0 / (sz[0] - 1);
//     spacing[1] = 1.0 / (sz[1] - 1);
//     double rCut = 3.9 * spacing[0];
//     Ghost<2, double> ghost(rCut);
//     BOOST_TEST_MESSAGE("Init vector_dist...");

//     vector_dist<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

//     double t=0.0,tf=0.5;
//     const double dt=0.1;

//     auto it = Particles.getGridIterator(sz);
//     while (it.isNext())
//     {
//         Particles.add();
//         auto key = it.get();
//         mem_id k0 = key.get(0);
//         double xp0 = k0 * spacing[0];
//         Particles.getLastPos()[0] = xp0;
//         mem_id k1 = key.get(1);
//         double yp0 = k1 * spacing[1];
//         Particles.getLastPos()[1] = yp0;
//         Particles.getLastProp<0>() = 1.0/(1.0+exp(-t)); // Carefull in putting the constant, f = A*sigmoid does not respect f' = f*(1.0-f) but f*(1.0-f/A), for simplicity I remove the constant
//         Particles.getLastProp<1>() = 1.0/(1.0+exp(-tf)); // Carefull in putting the constant, f = A*sigmoid does not respect f' = f*(1.0-f)  but f*(1.0-f/A), for simplicity I remove the constant
//         ++it;
//     }
//     Particles.map();
//     Particles.ghost_get<0>();
//     auto Init = getV<0>(Particles);
//     auto Sol = getV<1>(Particles);
//     auto OdeSol = getV<2>(Particles);

//     state_type x0;
//     x0=Init;
//     // The rhs of x' = f(x)
//     //size_t steps=boost::numeric::odeint::integrate(sigmoid,x0,0.0,tf,dt);
//     //typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type > > stepper_type;
//     //integrate_adaptive( stepper_type() , sigmoid , x0 , t , tf , dt);
//     size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< state_type >(),sigmoid,x0,0.0,tf,dt);

//     OdeSol=x0;
//     auto it2 = Particles.getDomainIterator();
//     double worst = 0.0;
//     while (it2.isNext()) {
//         auto p = it2.get();
//         if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst) {
//             worst = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
//         }
//         ++it2;
//     }

//     BOOST_REQUIRE(worst < 1e-8);

//     x0=Init;
//     boost::numeric::odeint::runge_kutta4< state_type > rk4;
//     for( size_t i=0 ; i<int(tf/dt) ; ++i,t+=dt )
//     {
//         rk4.do_step(sigmoid,x0,t,dt);
//         t+=dt;
//     }

//     OdeSol=x0;
//     auto it3 = Particles.getDomainIterator();
//     double worst2 = 0.0;
//     while (it3.isNext()) {
//         auto p = it3.get();
//         if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst2) {
//             worst2 = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
//         }
//         ++it3;
//     }
//     //std::cout<<worst2<<std::endl;
//     BOOST_REQUIRE(worst < 1e-6);
//     BOOST_REQUIRE_EQUAL(worst,worst2);
// }



#ifdef HAVE_EIGEN

BOOST_AUTO_TEST_CASE(dcpse_op_react_diff_test) {

  size_t edgeSemiSize{5};
  const size_t sz[dim] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
  Box<dim,double> box{{0.0, 0.0},{1.0, 1.0}};
  periodicity<dim> bc{{PERIODIC,PERIODIC}};
  double spacing[dim];
  spacing[0] = 1.0 / (sz[0]);
  spacing[1] = 1.0 / (sz[1]);
  Ghost<dim,double> ghost{spacing[0] * 3};

  BOOST_TEST_MESSAGE("Test: reaction diffusion");
  BOOST_TEST_MESSAGE("Init grid_dist_id ...");
  double sigma2 = spacing[0] * spacing[1] / (2 * 4);

  // properties: u, v, du, dv
  grid_dist_id<dim,double,aggregate<double,double,double,double,double,double>> domain{sz,box,ghost,bc};

  auto it{domain.getDomainIterator()};
  while (it.isNext()) {
    auto key{it.get()};
    domain.get<0>(key) = 0.0; // u
    domain.get<1>(key) = 0.0; // v
    domain.get<2>(key) = 0.0; // du/dt
    domain.get<3>(key) = 0.0; // dv/dt

    auto gkey = it.getGKey(key);

    if (gkey.get(0)==sz[0] / 2 && gkey.get(1) == sz[1]/2)
    {
      domain.get<0>(key) = 1.0;
      domain.get<1>(key) = 1.0;
    }
    ++it;
  }
  domain.ghost_get<0>();
  
  FD::Derivative<0,2,2,FD::CENTRAL> ddx;
  FD::Derivative<1,2,2,FD::CENTRAL> ddy;

  gridGlobal=(void *) & domain;

  auto u{FD::getV<0>(domain)};
  auto v{FD::getV<1>(domain)};
  auto fu{FD::getV<2>(domain)};
  auto fv{FD::getV<3>(domain)};
  
  Fitz<decltype(ddx),decltype(ddy)> system(ddx,ddy);
  state_type_2ode x0;
  
  x0.data.get<0>() = u;
  x0.data.get<1>() = v;

  double dt{0.001};
  double t{0.0};
  double tf{10.5};

  //typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;
  typedef boost::numeric::odeint::runge_kutta4<state_type_2ode,double,state_type_2ode,double,boost::numeric::odeint::vector_space_algebra_ofp> stepper_type;

  integrate_adaptive(stepper_type(),system,x0,t,tf,dt);
  fu = x0.data.get<0>();
  fv = x0.data.get<1>();

  domain.ghost_get<2,3>();
  u = ddx(fu) + ddy(fu);
  v = ddx(fv) + ddy(fv);

  auto it2{domain.getDomainIterator()};
  
  if (create_vcluster().rank() == 0)
    ++it2;
  
  while (it2.isNext()) {
    auto p{it2.get()};
    BOOST_REQUIRE_CLOSE(domain.get<0>(p),-1.0,1);   
    ++it2;
  }
}
#endif

BOOST_AUTO_TEST_SUITE_END()