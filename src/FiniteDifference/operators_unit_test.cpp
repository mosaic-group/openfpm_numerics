/*
 * operators_unit_tests.cpp
 *
 *  Created on: Dec 09, 2019
 *      Author: amfoggia
 */

//#define SE_CLASS1
#ifndef OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_UNIT_TEST_HPP_
#define OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_UNIT_TEST_HPP_

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "FiniteDifference/Derivative.hpp"
#include "FiniteDifference/Laplacian.hpp"
#include "FiniteDifference/Average.hpp"
#include "FiniteDifference/sum.hpp"
#include "FiniteDifference/mul.hpp"
#include "FiniteDifference/FDScheme.hpp"
#include "eq.hpp"
#include "FiniteDifference/operators.hpp"
#include "Vector/Vector.hpp"


struct  op_sys_nn
{
  // dimensionaly of the equation (2D problem 3D problem ...)
  static const unsigned int dims = 2;
  
  // number of fields in the system v_x, v_y, P so a total of 3
  static const unsigned int nvar = 1;
  
  // boundary conditions PERIODIC OR NON_PERIODIC
  static const bool boundary[];
  
  // type of space float, double, ...
  typedef float stype;
  
  // type of base grid, it is the distributed grid that will store the result
  // Note the first property is a 2D vector (velocity), the second is a scalar (Pressure)
  typedef grid_dist_id<2,float,aggregate<float>,CartDecomposition<2,float>> b_grid;
  
  // type of SparseMatrix, for the linear system, this parameter is bounded by the solver
  // that you are using, in case of umfpack using <double,int> it is the only possible choice
  // By default SparseMatrix is EIGEN based
  typedef SparseMatrix<double,int> SparseMatrix_type;
  
  // type of Vector for the linear system, this parameter is bounded by the solver
  // that you are using, in case of umfpack using <double> it is the only possible choice
  typedef Vector<double> Vector_type;
  
  // Define that the underline grid where we discretize the system of equation is staggered
  static const int grid_type = STAGGERED_GRID;
};


const bool op_sys_nn::boundary[] = {NON_PERIODIC,NON_PERIODIC};

constexpr unsigned int x = 0;
constexpr unsigned int y = 1;
constexpr unsigned int z = 2;
constexpr unsigned int V = 0;

template<typename T>
struct debug;

BOOST_AUTO_TEST_SUITE( operators_test )

BOOST_AUTO_TEST_CASE( operator_plus )
{
  Field<x,op_sys_nn> f1({-3,2});
  Field<y,op_sys_nn> f2;
  Field<z,op_sys_nn> f3;

  std::cout << "f1 pos: " << f1.def_pos.to_string() << std::endl;
  std::cout << f2.def_pos.to_string() << std::endl;

  auto sum2 = f1+f2;
  auto sum3 = f1+f2+f3;
  
  auto mul2 = f1*f2;
  auto mul3 = f2*f1*f3;

  auto test1 = f2*(f3+f1);

  Der<x,op_sys_nn,CENTRAL> dx;
  auto der1 = dx(f1);
  auto der2 = dx(sum2);

  coeff<double,op_sys_nn> c{3.0};
  auto coeff_test = c*f1;

  // debug<decltype(f1*f1)> dbg;
  
  // ----------------------------------------------------------------------
  
  Box<2,float> domain({0.0,0.0},{1.0,1.0});
  Ghost<2,float> g(0.01);
  size_t szu[] = {6,6};
  long int sz[] = {6,6};
  periodicity<op_sys_nn::dims> periodicity = {NON_PERIODIC,NON_PERIODIC};
  
  grid_dist_id<2,float,aggregate<float>> g_dist(szu,domain,g);
  Padding<2> pd({1,1},{0,0});
  Ghost<2,long int> stencil_max(1);
  FDScheme<op_sys_nn> fd(pd,stencil_max,domain, g_dist);
 
  auto it = g_dist.getDomainGhostIterator();

  // while (it.isNext()) {
    
  //   auto kmap = it.get();
  //   auto orig = kmap;
  //   std::cout << kmap.getKey().to_string() << " ----- ";

  //   long int old_val = kmap.getKeyRef().get(0);
  //   kmap.getKeyRef().set_d(0, kmap.getKeyRef().get(0) + 1);
  //   std::cout << kmap.getKey().to_string() << " ----- ";
    
  //   kmap.getKeyRef().set_d(1, kmap.getKeyRef().get(1) + 1);
  //   std::cout << kmap.getKey().to_string() << " -----> ";
    
  //   kmap.getKeyRef() = orig.getKeyRef();
  //   std::cout << kmap.getKey().to_string() << "\n";
    
  //   ++it;
  // }

  // ----------------------------------------------------------------------

  std::initializer_list<char> cc = {0,0};
  std::initializer_list<char> ll = {0,-1};
  std::initializer_list<char> bl = {-1,-1};
  std::initializer_list<char> bb = {-1,0};
  
  Field<y,op_sys_nn> vy{bb};
  grid_sm<op_sys_nn::dims,void> gs = g_dist.getGridInfoVoid();
  float spacing[2];
  spacing[0] = 0.1;
  spacing[1] = 0.1;
  std::unordered_map<long int,float> cols;
  float coeff = 1.0;
  // auto & test_grid = fd.getMap();
  Avg<x,decltype(vy),FORWARD> avg_vy{vy};

  std::cout << "------------------------\n";

  //fd.impose<y>(avg_vy,0.0,{-1,0},{-1,sz[1]-1},bb);

  fd.impose<y>(vy,0.0,{0,0},{0,sz[1]-1},bl);
  typedef typename op_sys_nn::SparseMatrix_type::triplet_type triplet;
  openfpm::vector<triplet> & trpl = fd.getA().getMatrixTriplets();

  for (int i = 0; i < trpl.size(); ++i)
    std::cout << "row: " << trpl.get(i).row() << " "
	      << "col: " << trpl.get(i).col() << " "
	      << "val: " << trpl.get(i).value() << "\n-----------------------" << std::endl;
  

  // auto it2 = g_dist.getDomainIterator();

  // while(it2.isNext()) {
  //   auto kmap = it2.get();
  //   std::cout << "key: " << kmap.getKey().to_string() << " index: " << test_grid.template get<0>(kmap)*op_sys_nn::nvar + 0 << std::endl;
  //   ++it2;
  // }
  

  // auto it2 = test_grid.getDomainIterator();
  // auto kmap = it2.get();
  // vx.value(test_grid,kmap,gs,spacing,cols,coeff,{-1,-1});

  // std::cout << "----------\n";
  // for (auto cols_it = cols.cbegin(); cols_it != cols.cend(); ++cols_it)
  //   std::cout << "Key:[" << cols_it->first << "] Value:[" << cols_it->second << "]\n";
  // std::cout << "----------\n";
  
  // ++it2;
  // ++it2;
  // ++it2;

  // kmap = it2.get();
  // vx.value(test_grid,kmap,gs,spacing,cols,coeff,{-1,-1});
  // std::cout << it2.get().getKey().to_string() << std::endl;

  // std::cout << "----------\n";  
  // for (auto cols_it = cols.cbegin(); cols_it != cols.cend(); ++cols_it)
  //   std::cout << "Key:[" << cols_it->first << "] Value:[" << cols_it->second << "]\n";
  // std::cout << "----------\n";  
  
  
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* OPENFPM_NUMERICS_SRC_FINITEDIFFERENCE_OPERATORS_UNIT_TEST_HPP_ */

