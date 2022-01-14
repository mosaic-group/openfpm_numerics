
#include "config.h"
#include <type_traits>
#include <cstring>
#include "util/common.hpp"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "Grid/grid_dist_id.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include "OdeIntegrators/boost_vector_algebra_ofp.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "util/PathsAndFiles.hpp"

// Include own files

// Global variables
constexpr int DIM = 2;
constexpr int CONC = 0;
constexpr int RHS = 1;
constexpr int ANALYT_CONC = 2;
constexpr int ABS_ERROR = 3;

typedef aggregate<double,double,double,double> prop;
typedef grid_dist_id<DIM,double,prop> grid_type;
void *gridGlobal; // Pointer to grid, so one can use it everywhere

// Class that contains the system to solve
template <typename DX, typename DY>
class diff {

  double coeff;
  DX ddx;
  DY ddy;

public:
  diff(double m_coeff, DX & m_ddx, DY & m_ddy)
    : coeff(m_coeff),
      ddx(m_ddx),
      ddy(m_ddy) {}

  void operator() (const state_type_ofpm_impl<1,texp_g<2,double>> & x,
		   state_type_ofpm_impl<1,texp_g<2,double>> & dxdt,
		   const double t) {

    // Get global pointer for grid
    grid_type & grid = *(grid_type*) gridGlobal;
    
    // Expressions
    auto conc{FD::getV<CONC>(grid)};
    auto rhs{FD::getV<RHS>(grid)};

    // Initialize the concentration from the state
    conc = x.data.get<0>();

    grid.template ghost_get<CONC>();

    // System
    rhs = coeff * (ddx(conc) + ddy(conc));

    // Copy back the rhs into the state
    dxdt.data.get<0>() = rhs;
    
  }
};

// Struct that has the observer in which the state will be written in a file
// and the error will be computed
struct boundary_observer {
  
public:
  boundary_observer() {}
  void operator() (const state_type_ofpm_impl<1,texp_g<2,double>> & x,
		   double t) {

    // Get global pointer for grid
    grid_type & grid = *(grid_type*) gridGlobal;
    
    // Impose BC if needed

    // Compute Linf error
    double maxError{0.0};
    auto it{grid.getDomainIterator()};
    while (it.isNext()) {
      auto key{it.get()};
      grid.template get<ABS_ERROR>(key) = std::fabs(grid.template get<ANALYT_CONC>(key) - grid.template get<CONC>(key));
      if (grid.template get<ABS_ERROR>(key) > maxError)
	maxError = grid.template get<ABS_ERROR>(key);
      ++it;
    }
    
    Vcluster<> & v_cl = create_vcluster();
    v_cl.max(maxError);
    v_cl.execute();

    std::cout << "t = " << t << " Linf: " << maxError << std::endl;
  }
  
};

BOOST_AUTO_TEST_SUITE(odeInt_grid_tests)

BOOST_AUTO_TEST_CASE(odeint_for_grid) 
{
    state_type_ofpm_impl<3,texp_g<3,double>> state_;

    bool test = std::is_same<decltype(state_.data),aggregate<texp_g<3,double>,texp_g<3,double>,texp_g<3,double>> >::value;
    BOOST_REQUIRE_EQUAL(test,true);

  Vcluster<> & v_cl = create_vcluster();

  int n;
  double finalT, dt;

    n = 100;
    finalT = 0.1;
    dt = 0.00001;

  std::string cwd{boost::filesystem::current_path().string()};
  std::string output_dir{cwd + "/output_odeint"};
  create_directory_if_not_exist(output_dir);
  std::string g_output{output_dir + "/test_diff2D"};

  // Grid creation
  Box<DIM,double> domain{{-1.0,-1.0},{1.0,1.0}};
  Ghost<DIM,long int> ghost{3};
  // periodicity<DIM> bc{{NON_PERIODIC,NON_PERIODIC}};
  periodicity<DIM> bc{{PERIODIC,PERIODIC}};
  size_t sz[DIM];

  sz[0] = n;
  sz[1] = n;
  grid_type grid{sz,domain,ghost,bc};
  grid.setPropNames({"num_conc","analyt_conc","abs_error"});

  // Simulation parameters
  double dx{grid.getSpacing()[0]};
  int nT{int(std::floor(finalT/dt)+1)};
  double diff_coeff{1.0};
  
  // Analytical solution
  {
    auto it{grid.getDomainIterator()};
    while (it.isNext()) {
      auto key{it.get()};
      double x{grid.getPos(key)[0]};
      grid.template get<ANALYT_CONC>(key) = -std::sin(x)*std::exp(-finalT*diff_coeff);
      ++it;
    }
  }
  
  // Initial condition
  {
    auto it{grid.getDomainIterator()};
    while (it.isNext()) {
      auto key{it.get()};
      double x{grid.getPos(key)[0]};
      grid.template get<CONC>(key) = -std::sin(x);
      ++it;
    }
  }

  // Initialize the global pointer for grid
  gridGlobal = (void *) & grid;
  
  // Define FD objects needed
  FD::Derivative<0,2,2,FD::CENTRAL> ddx;
  FD::Derivative<1,2,2,FD::CENTRAL> ddy;

  // Objects
  diff<decltype(ddx),decltype(ddy)> system{diff_coeff,ddx,ddy};                                                                             // System
  boundary_observer observer;                                                                                                               // Observer
  state_type_ofpm_impl<1,texp_g<2,double>> state;                                                                                                                  // State
  boost::numeric::odeint::euler< state_type_ofpm_impl<1,texp_g<2,double>> ,double, state_type_ofpm_impl<1,texp_g<2,double>> ,double,boost::numeric::odeint::vector_space_algebra_ofp> euler; // Time integrator

  // State and integration
  auto conc{FD::getV<CONC>(grid)};
  state.data.get<0>() = conc;
  size_t steps{boost::numeric::odeint::integrate_const(euler,system,state,0.0,finalT,dt,observer)};
  
  grid.write_frame(g_output,0,VTK_WRITER);
}

BOOST_AUTO_TEST_SUITE_END()
