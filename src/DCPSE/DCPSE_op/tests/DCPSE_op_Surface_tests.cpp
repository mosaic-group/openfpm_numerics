//
// Created by Abhinav Singh on 15.11.21.
//
//#define SE_CLASS1

#include "config.h"
#ifdef HAVE_EIGEN
#ifdef HAVE_PETSC

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../DCPSE_surface_op.hpp"
#include "../DCPSE_Solver.hpp"
#include "../../DcpseInterpolation.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include <iostream>
#include "util/SphericalHarmonics.hpp"

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)
    BOOST_AUTO_TEST_CASE(dcpse_surface_simple) {
        double boxP1{-1.5}, boxP2{1.5};
        double boxSize{boxP2 - boxP1};
        size_t n=256;
        size_t sz[2] = {n,n};
        double grid_spacing{boxSize/(sz[0]-1)};
        double rCut{3.9 * grid_spacing};

        Box<2,double> domain{{boxP1,boxP1},{boxP2,boxP2}};
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost{rCut + grid_spacing/8.0};
        auto &v_cl=create_vcluster();

        vector_dist_ws<2, double, aggregate<double,double,double[2],double,double[2]>> Sparticles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // 1. Particles on a line
        if (v_cl.rank() == 0) {
            for (int i = 0; i < n; ++i) {
                double xp = -1.5+i*grid_spacing;
                Sparticles.add();
                Sparticles.getLastPos()[0] = xp;
                Sparticles.getLastPos()[1] = 0;
                Sparticles.getLastProp<3>() = std::sin(xp);
                Sparticles.getLastProp<2>()[0] = 0;
                Sparticles.getLastProp<2>()[1] = 1.0;
                Sparticles.getLastProp<1>() = -std::sin(xp);//sin(theta)*exp(-finalT/(radius*radius));
                Sparticles.getLastSubset(0);
            }
        }
        Sparticles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Sparticles.ghost_get<0,3>();
        //Sparticles.write("Sparticles");
        //Here template parameters are Normal property no.
        SurfaceDerivative_xx<2> SDxx(Sparticles, 2, rCut,grid_spacing);
        SurfaceDerivative_yy<2> SDyy(Sparticles, 2, rCut,grid_spacing);
        //SurfaceDerivative_x<2> SDx(Sparticles, 4, rCut,grid_spacing);
        //SurfaceDerivative_y<2> SDy(Sparticles, 4, rCut,grid_spacing);
        auto INICONC = getV<3>(Sparticles);
        auto CONC = getV<0>(Sparticles);
        auto TEMP = getV<4>(Sparticles);
        auto normal = getV<2>(Sparticles);
        //auto ANASOL = getV<1>(domain);
        CONC=SDxx(INICONC)+SDyy(INICONC);
        //TEMP[0]=(-normal[0]*normal[0]+1.0) * SDx(INICONC) - normal[0]*normal[1] * SDy(INICONC);
        //TEMP[1]=(-normal[1]*normal[1]+1.0) * SDy(INICONC) - normal[0]*normal[1] * SDx(INICONC);
        //Sparticles.ghost_get<4>();
        //CONC=SDxx(TEMP[0]) + SDyy(TEMP[1]);
        auto it2 = Sparticles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            }

            ++it2;
        }
        Sparticles.deleteGhost();
        //std::cout<<v_cl.rank()<<":WORST:"<<worst<<std::endl;
        //Sparticles.write("Sparticles");
        BOOST_REQUIRE(worst < 0.03);
}
    BOOST_AUTO_TEST_CASE(dcpse_surface_circle) {

        double boxP1{-1.5}, boxP2{1.5};
        double boxSize{boxP2 - boxP1};
        size_t n=512;
        auto &v_cl=create_vcluster();
        //std::cout<<v_cl.rank()<<":Enter res: "<<std::endl;
        //std::cin>>n;
        size_t sz[2] = {n,n};
        double grid_spacing{boxSize/(sz[0]-1)};
        double rCut{5.1 * grid_spacing};

        Box<2,double> domain{{boxP1,boxP1},{boxP2,boxP2}};
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost{rCut + grid_spacing/8.0};
        vector_dist_ws<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // Surface prameters
        const double radius{1.0};
        std::array<double,2> center{0.0,0.0};
        Point<2,double> coord;
        const double pi{3.14159265358979323846};

        // 1. Particles on surface
        double theta{0.0};
        double dtheta{2*pi/double(n)};
        if (v_cl.rank() == 0) {
            for (int i = 0; i < n; ++i) {
                coord[0] = center[0] + radius * std::cos(theta);
                coord[1] = center[1] + radius * std::sin(theta);
                Sparticles.add();
                Sparticles.getLastPos()[0] = coord[0];
                Sparticles.getLastPos()[1] = coord[1];
                Sparticles.getLastProp<3>() = std::sin(theta);
                Sparticles.getLastProp<2>()[0] = std::cos(theta);
                Sparticles.getLastProp<2>()[1] = std::sin(theta);
                Sparticles.getLastProp<1>() = -std::sin(theta);;//sin(theta)*exp(-finalT/(radius*radius));
                Sparticles.getLastSubset(0);
                theta += dtheta;
            }
        }
        Sparticles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Sparticles.ghost_get<0,3>();
        //Sparticles.write("Sparticles");
        //Here template parameters are Normal property no.
        SurfaceDerivative_xx<2> SDxx(Sparticles, 2, rCut,grid_spacing);
        SurfaceDerivative_yy<2> SDyy(Sparticles, 2, rCut,grid_spacing);
        //SurfaceDerivative_xy<2> SDxy(Sparticles, 3, rCut,grid_spacing);
        //SurfaceDerivative_x<2> SDx(Sparticles, 3, rCut,grid_spacing);
        //SurfaceDerivative_y<2> SDy(Sparticles, 3, rCut,grid_spacing);
        auto INICONC = getV<3>(Sparticles);
        auto CONC = getV<0>(Sparticles);
        auto TEMP = getV<4>(Sparticles);
        auto normal = getV<2>(Sparticles);

        CONC=SDxx(INICONC)+SDyy(INICONC);
        //TEMP[0]=(-normal[0]*normal[0]+1.0) * SDx(INICONC) - normal[0]*normal[1] * SDy(INICONC);
        //TEMP[1]=(-normal[1]*normal[1]+1.0) * SDy(INICONC) - normal[0]*normal[1] * SDx(INICONC);
        //TEMP[0]=(-normal[0]*normal[0]+1.0);
        //TEMP[1]=normal[0]*normal[1];
        //Sparticles.ghost_get<2,4>();
        //CONC=SDx(TEMP[0]) + SDy(TEMP[1]);
        //CONC= (SDx(TEMP[0])*SDx(INICONC)+TEMP[0]*SDxx(INICONC)-(SDx(TEMP[1])*SDy(INICONC)+TEMP[1]*SDxy(INICONC)))+
        //        (SDy((-normal[1]*normal[1]+1.0))*SDy(INICONC)+(-normal[1]*normal[1]+1.0)*SDyy(INICONC)-(SDy(TEMP[1])*SDx(INICONC)+TEMP[1]*SDxy(INICONC)));
        auto it2 = Sparticles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            }
            ++it2;
        }
        Sparticles.deleteGhost();
        //Sparticles.write("Sparticles");
        //std::cout<<worst;
        BOOST_REQUIRE(worst < 0.03);
}
    BOOST_AUTO_TEST_CASE(dcpse_surface_solver_circle) {
        double boxP1{-1.5}, boxP2{1.5};
        double boxSize{boxP2 - boxP1};
        size_t n=512,k=2;
        auto &v_cl=create_vcluster();
        /*if(v_cl.rank()==0)
        std::cout<<v_cl.rank()<<":Enter res: "<<std::endl;
        std::cin>>n;
        if(v_cl.rank()==0)
        std::cout<<v_cl.rank()<<":Enter Freq: "<<std::endl;
        std::cin>>k;*/
        size_t sz[2] = {n,n};
        double grid_spacing{boxSize/(sz[0]-1)};
        double rCut{3.9 * grid_spacing};

        Box<2,double> domain{{boxP1,boxP1},{boxP2,boxP2}};
        size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
        Ghost<2,double> ghost{rCut + grid_spacing/8.0};
        vector_dist_ws<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles(0, domain,bc,ghost);
        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        // Surface prameters
        const double radius{1.0};
        std::array<double,2> center{0.0,0.0};
        Point<2,double> coord;
        const double pi{3.14159265358979323846};

        // 1. Particles on surface
        double theta{0.0};
        double dtheta{2*pi/double(n)};
        if (v_cl.rank() == 0) {
            for (int i = 0; i < n; ++i) {
                coord[0] = center[0] + radius * std::cos(theta);
                coord[1] = center[1] + radius * std::sin(theta);
                Sparticles.add();
                Sparticles.getLastPos()[0] = coord[0];
                Sparticles.getLastPos()[1] = coord[1];
                Sparticles.getLastProp<3>() = -openfpm::math::intpowlog(k,2)*std::sin(k*theta);
                Sparticles.getLastProp<2>()[0] = std::cos(theta);
                Sparticles.getLastProp<2>()[1] = std::sin(theta);
                Sparticles.getLastProp<1>() = std::sin(k*theta);;//sin(theta)*exp(-finalT/(radius*radius));
                Sparticles.getLastSubset(0);
                if(coord[0]==1. && coord[1]==0.)
                {Sparticles.getLastSubset(1);}
                theta += dtheta;
            }
        }
        Sparticles.map();

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        Sparticles.ghost_get<0>();
        //Sparticles.write("Sparticles");
        vector_dist_subset<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles_bulk(Sparticles,0);
        vector_dist_subset<2, double, aggregate<double,double,double[2],double,double[2],double>> Sparticles_boundary(Sparticles,1);
        auto & bulk=Sparticles_bulk.getIds();
        auto & boundary=Sparticles_boundary.getIds();
        //Here template parameters are Normal property no.
        SurfaceDerivative_xx<2> SDxx(Sparticles, 2, rCut,grid_spacing);
        SurfaceDerivative_yy<2> SDyy(Sparticles, 2, rCut,grid_spacing);
        auto INICONC = getV<3>(Sparticles);
        auto CONC = getV<0>(Sparticles);
        auto TEMP = getV<4>(Sparticles);
        auto normal = getV<2>(Sparticles);
        auto ANASOL = getV<1>(Sparticles);
        DCPSE_scheme<equations2d1,decltype(Sparticles)> Solver(Sparticles);
        Solver.impose(SDxx(CONC)+SDyy(CONC), bulk, INICONC);
        Solver.impose(CONC, boundary, ANASOL);
        Solver.solve(CONC);
        auto it2 = Sparticles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
            }
            ++it2;
        }
        Sparticles.deleteGhost();
        //Sparticles.write("Sparticles");
        //std::cout<<worst;
        BOOST_REQUIRE(worst < 0.03);
}
BOOST_AUTO_TEST_CASE(dcpse_surface_sphere_copy) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n=10000;
  size_t n_sp=n;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n,n,n};
  double grid_spacing{boxSize/(sz[0]-1)};
  double grid_spacing_surf=std::sqrt(4.0*M_PI/n);//grid_spacing*30;
  double rCut{2.5 * grid_spacing_surf};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};

  constexpr int K = 1;
  // particles
  vector_dist_ws<3, double, aggregate<double,double,double[3],double,double[3],double>> Sparticles(0, domain,bc,ghost);
  // 1. particles on the Spherical surface
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    //std::vector<Vector3f> data;
    //GenerateSphere(1,data);
    for(int i=1;i<n_sp;i++)
        {
            double y = 1.0 - (i /double(n_sp - 1.0)) * 2.0;
            double radius = sqrt(1 - y * y);
            double Golden_theta = Golden_angle * i;
            double x = cos(Golden_theta) * radius;
            double z = sin(Golden_theta) * radius;
            Sparticles.add();
            Sparticles.getLastPos()[0] = x;
            Sparticles.getLastPos()[1] = y;
            Sparticles.getLastPos()[2] = z;
            double rm=sqrt(x*x+y*y+z*z);
            Sparticles.getLastProp<2>()[0] = x/rm;
            Sparticles.getLastProp<2>()[1] = y/rm;
            Sparticles.getLastProp<2>()[2] = z/rm;
            Sparticles.getLastProp<4>()[0] = 1.0 ;
            Sparticles.getLastProp<4>()[1] = std::atan2(sqrt(x*x+y*y),z);
            Sparticles.getLastProp<4>()[2] = std::atan2(y,x);
            if(i<=2*(K)+1)
            {Sparticles.getLastSubset(1);}
            else
            {Sparticles.getLastSubset(0);}
        }
    //std::cout << "n: " << n << " - grid spacing: " << grid_spacing << " - rCut: " << rCut << "Surf Normal spacing" << grid_spacing<<std::endl;
  }

  Sparticles.map();
  Sparticles.ghost_get<3>();

  vector_dist_subset<3,double,aggregate<double,double,double[3],double,double[3],double>> Sparticles_bulk(Sparticles,0);
  vector_dist_subset<3,double,aggregate<double,double,double[3],double,double[3],double>> Sparticles_boundary(Sparticles,1);
  auto &bulkIds=Sparticles_bulk.getIds();
  auto &bdrIds=Sparticles_boundary.getIds();
  std::unordered_map<const lm,double,key_hash,key_equal> Alm;
  //Setting max mode l_max
  //Setting amplitudes to 1
  for(int l=0;l<=K;l++){
      for(int m=-l;m<=l;m++){
          Alm[std::make_tuple(l,m)]=0;
      }
  }
  Alm[std::make_tuple(1,0)]=1;
  auto it2 = Sparticles.getDomainIterator();
  while (it2.isNext()) {
      auto p = it2.get();
      Point<3, double> xP = Sparticles.getProp<4>(p);
      /*double Sum=0;
      for(int m=-spL;m<=spL;++m)
      {
        Sum+=openfpm::math::Y(spL,m,xP[1],xP[2]);
      }*/
      //Sparticles.getProp<ANADF>(p) = Sum;//openfpm::math::Y(K,K,xP[1],xP[2]);openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);;
      Sparticles.getProp<3>(p)=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<1>(p)=2.0;//-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      ++it2;
  }
  auto f=getV<2>(Sparticles);
  auto Df=getV<0>(Sparticles);

  SurfaceDerivative_x<2> Sdx{Sparticles,2,rCut,grid_spacing_surf};
  SurfaceDerivative_y<2> Sdy{Sparticles,2,rCut,grid_spacing_surf};
  SurfaceDerivative_z<2> Sdz{Sparticles,2,rCut,grid_spacing_surf};
  //Laplace_Beltrami<2> SLap{Sparticles,2,rCut,grid_spacing_surf};
  //Sdyy.DrawKernel<5>(Sparticles,0);
  //Sdzz.DrawKernel<5>(Sparticles,0);
/*  std::cout<<"SDXX:"<<std::endl;
  Sdxx.checkMomenta(Sparticles);
  std::cout<<"SDYY:"<<std::endl;
  Sdyy.checkMomenta(Sparticles);
  std::cout<<"SDZZ:"<<std::endl;
  Sdzz.checkMomenta(Sparticles);*/

  Sparticles.ghost_get<3>();
  Df=(Sdx(f[0])+Sdy(f[1])+Sdz(f[2]));
  //Df=SLap(f);
  auto it3 = Sparticles.getDomainIterator();
  double worst = 0.0;
  while (it3.isNext()) {
      auto p = it3.get();
      //Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
      if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
          worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
      }
      ++it3;
  }
        Sparticles.deleteGhost();
        //Sparticles.write("Sparticles");
        std::cout<<worst;
	BOOST_REQUIRE(worst < 0.03);
}
BOOST_AUTO_TEST_CASE(dcpse_surface_sphere) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n=40960;
  size_t n_sp=n;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n,n,n};
  double grid_spacing{boxSize/(sz[0]-1)};
  double grid_spacing_surf=std::sqrt(4.0*M_PI/n);//grid_spacing*30;
  double rCut{2.5 * grid_spacing_surf};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};

  constexpr int K = 1;
  // particles
  vector_dist_ws<3, double, aggregate<double,double,double[3],double,double[3],double>> Sparticles(0, domain,bc,ghost);
  // 1. particles on the Spherical surface
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    //std::vector<Vector3f> data;
    //GenerateSphere(1,data);
    for(int i=1;i<n_sp;i++)
        {
            double y = 1.0 - (i /double(n_sp - 1.0)) * 2.0;
            double radius = sqrt(1 - y * y);
            double Golden_theta = Golden_angle * i;
            double x = cos(Golden_theta) * radius;
            double z = sin(Golden_theta) * radius;
            Sparticles.add();
            Sparticles.getLastPos()[0] = x;
            Sparticles.getLastPos()[1] = y;
            Sparticles.getLastPos()[2] = z;
            double rm=sqrt(x*x+y*y+z*z);
            Sparticles.getLastProp<2>()[0] = x/rm;
            Sparticles.getLastProp<2>()[1] = y/rm;
            Sparticles.getLastProp<2>()[2] = z/rm;
            Sparticles.getLastProp<4>()[0] = 1.0 ;
            Sparticles.getLastProp<4>()[1] = std::atan2(sqrt(x*x+y*y),z);
            Sparticles.getLastProp<4>()[2] = std::atan2(y,x);
            if(i<=2*(K)+1)
            {Sparticles.getLastSubset(1);}
            else
            {Sparticles.getLastSubset(0);}
        }
    //std::cout << "n: " << n << " - grid spacing: " << grid_spacing << " - rCut: " << rCut << "Surf Normal spacing" << grid_spacing<<std::endl;
  }

  Sparticles.map();
  Sparticles.ghost_get<3>();

  vector_dist_subset<3,double,aggregate<double,double,double[3],double,double[3],double>> Sparticles_bulk(Sparticles,0);
  vector_dist_subset<3,double,aggregate<double,double,double[3],double,double[3],double>> Sparticles_boundary(Sparticles,1);
  auto &bulkIds=Sparticles_bulk.getIds();
  auto &bdrIds=Sparticles_boundary.getIds();
  std::unordered_map<const lm,double,key_hash,key_equal> Alm;
  //Setting max mode l_max
  //Setting amplitudes to 1
  for(int l=0;l<=K;l++){
      for(int m=-l;m<=l;m++){
          Alm[std::make_tuple(l,m)]=0;
      }
  }
  Alm[std::make_tuple(1,0)]=1;
  auto it2 = Sparticles.getDomainIterator();
  while (it2.isNext()) {
      auto p = it2.get();
      Point<3, double> xP = Sparticles.getProp<4>(p);
      /*double Sum=0;
      for(int m=-spL;m<=spL;++m)
      {
        Sum+=openfpm::math::Y(spL,m,xP[1],xP[2]);
      }*/
      //Sparticles.getProp<ANADF>(p) = Sum;//openfpm::math::Y(K,K,xP[1],xP[2]);openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);;
      Sparticles.getProp<3>(p)=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<1>(p)=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      ++it2;
  }
  auto f=getV<3>(Sparticles);
  auto Df=getV<0>(Sparticles);

  SurfaceDerivative_xx<2> Sdxx{Sparticles,2,rCut,grid_spacing_surf};
  SurfaceDerivative_yy<2> Sdyy{Sparticles,2,rCut,grid_spacing_surf};
  SurfaceDerivative_zz<2> Sdzz{Sparticles,2,rCut,grid_spacing_surf};
  //Laplace_Beltrami<2> SLap{Sparticles,2,rCut,grid_spacing_surf};
  //Sdyy.DrawKernel<5>(Sparticles,0);
  //Sdzz.DrawKernel<5>(Sparticles,0);
/*  std::cout<<"SDXX:"<<std::endl;
  Sdxx.checkMomenta(Sparticles);
  std::cout<<"SDYY:"<<std::endl;
  Sdyy.checkMomenta(Sparticles);
  std::cout<<"SDZZ:"<<std::endl;
  Sdzz.checkMomenta(Sparticles);*/

  Sparticles.ghost_get<3>();
  Df=(Sdxx(f)+Sdyy(f)+Sdzz(f));
  //Df=SLap(f);
  auto it3 = Sparticles.getDomainIterator();
  double worst = 0.0;
  while (it3.isNext()) {
      auto p = it3.get();
      //Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
      if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
          worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
      }
      ++it3;
  }
        Sparticles.deleteGhost();
        //Sparticles.write("Sparticles");
        std::cout<<worst;
	BOOST_REQUIRE(worst < 0.03);
}


BOOST_AUTO_TEST_CASE(dcpse_surface_sphere_old) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n=512;
  size_t n_sp=n;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n,n,n};
  double grid_spacing{boxSize/(sz[0]-1)};
  double grid_spacing_surf=grid_spacing*30;
  double rCut{2.5 * grid_spacing_surf};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};

  constexpr int K = 1;
  // particles
  vector_dist_ws<3, double, aggregate<double,double,double[3],double,double[3],double>> Sparticles(0, domain,bc,ghost);
  // 1. particles on the Spherical surface
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    //std::vector<Vector3f> data;
    //GenerateSphere(1,data);
    std::unordered_map<const lm,double,key_hash,key_equal> Alm;
          //Setting max mode l_max
          //Setting amplitudes to 1
          for(int l=0;l<=K;l++){
              for(int m=-l;m<=l;m++){
                  Alm[std::make_tuple(l,m)]=0;
              }
          }
    Alm[std::make_tuple(1,0)]=1;
    for(int i=1;i<n_sp;i++)
        {
            double y = 1.0 - (i /double(n_sp - 1.0)) * 2.0;
            double radius = sqrt(1 - y * y);
            double Golden_theta = Golden_angle * i;
            double x = cos(Golden_theta) * radius;
            double z = sin(Golden_theta) * radius;
            Sparticles.add();
            Sparticles.getLastPos()[0] = x;
            Sparticles.getLastPos()[1] = y;
            Sparticles.getLastPos()[2] = z;
            double rm=sqrt(x*x+y*y+z*z);
            Sparticles.getLastProp<2>()[0] = x/rm;
            Sparticles.getLastProp<2>()[1] = y/rm;
            Sparticles.getLastProp<2>()[2] = z/rm;
            Sparticles.getLastProp<4>()[0] = 1.0 ;
            Sparticles.getLastProp<4>()[1] = std::atan2(sqrt(x*x+y*y),z);
            Sparticles.getLastProp<4>()[2] = std::atan2(y,x);
            double m1=openfpm::math::sumY_Scalar<K>(1.0,std::atan2(sqrt(x*x+y*y),z),std::atan2(y,x),Alm);
            double m2=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(1.0,std::atan2(sqrt(x*x+y*y),z),std::atan2(y,x),Alm);
            Sparticles.getLastProp<3>()=m1;
            Sparticles.getLastProp<1>()=m2;
            Sparticles.getLastSubset(0);
            for(int j=1;j<=2;++j){
                Sparticles.add();
            Sparticles.getLastPos()[0] = x+j*grid_spacing_surf*x/rm;
            Sparticles.getLastPos()[1] = y+j*grid_spacing_surf*y/rm;
            Sparticles.getLastPos()[2] = z+j*grid_spacing_surf*z/rm;
            Sparticles.getLastProp<3>()=m1;
            Sparticles.getLastSubset(1);
            //Sparticles.getLastProp<1>(p)=m2;
            Sparticles.add();
            Sparticles.getLastPos()[0] = x-j*grid_spacing_surf*x/rm;
            Sparticles.getLastPos()[1] = y-j*grid_spacing_surf*y/rm;
            Sparticles.getLastPos()[2] = z-j*grid_spacing_surf*z/rm;
            Sparticles.getLastProp<3>()=m1;
            Sparticles.getLastSubset(1);
            //Sparticles.getLastProp<1>(p)=m2;
            }
        }
    //std::cout << "n: " << n << " - grid spacing: " << grid_spacing << " - rCut: " << rCut << "Surf Normal spacing" << grid_spacing<<std::endl;
  }

  Sparticles.map();
  Sparticles.ghost_get<3>();
  //Sparticles.write("SparticlesInit");

  vector_dist_subset<3,double,aggregate<double,double,double[3],double,double[3],double>> Sparticles_bulk(Sparticles,0);
  vector_dist_subset<3,double,aggregate<double,double,double[3],double,double[3],double>> Sparticles_boundary(Sparticles,1);
  auto &bulkIds=Sparticles_bulk.getIds();
  auto &bdrIds=Sparticles_boundary.getIds();
  /*auto it2 = Sparticles.getDomainIterator();
  while (it2.isNext()) {
      auto p = it2.get();
      Point<3, double> xP = Sparticles.getProp<4>(p);
      *//*double Sum=0;
      for(int m=-spL;m<=spL;++m)
      {
        Sum+=openfpm::math::Y(spL,m,xP[1],xP[2]);
      }*//*
      //Sparticles.getProp<ANADF>(p) = Sum;//openfpm::math::Y(K,K,xP[1],xP[2]);openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);;
      Sparticles.getProp<3>(p)=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<1>(p)=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      ++it2;
  }*/
  auto f=getV<3>(Sparticles);
  auto Df=getV<0>(Sparticles);

  //SurfaceDerivative_xx<2> Sdxx{Sparticles,2,rCut,grid_spacing_surf};
  //SurfaceDerivative_yy<2> Sdyy{Sparticles,2,rCut,grid_spacing_surf};
  //SurfaceDerivative_zz<2> Sdzz{Sparticles,2,rCut,grid_spacing_surf};
  Derivative_xx Sdxx{Sparticles,2,rCut};
  //std::cout<<"Dxx Done"<<std::endl;
  Derivative_yy Sdyy{Sparticles,2,rCut};
  //std::cout<<"Dyy Done"<<std::endl;
  Derivative_zz Sdzz{Sparticles,2,rCut};
  //std::cout<<"Dzz Done"<<std::endl;

  //Laplace_Beltrami<2> SLap{Sparticles,2,rCut,grid_spacing_surf};
  //SLap.DrawKernel<5>(Sparticles,73);
  //Sdxx.DrawKernel<5>(Sparticles,0);
  //Sdyy.DrawKernel<5>(Sparticles,0);
  //Sdzz.DrawKernel<5>(Sparticles,0);
/*  std::cout<<"SDXX:"<<std::endl;
  Sdxx.checkMomenta(Sparticles);
  std::cout<<"SDYY:"<<std::endl;
  Sdyy.checkMomenta(Sparticles);
  std::cout<<"SDZZ:"<<std::endl;
  Sdzz.checkMomenta(Sparticles);*/

  Sparticles.ghost_get<3>();
  Df=(Sdxx(f)+Sdyy(f)+Sdzz(f));
  //Df=SLap(f);
  //auto it3 = Sparticles_bulk.getDomainIterator();
  double worst = 0.0;
  for (int j = 0; j < bulkIds.size(); j++) {
      auto p = bulkIds.get<0>(j);
      //Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
        if (fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p)) > worst) {
                  worst = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
              }
  }
        Sparticles.deleteGhost();
        //Sparticles.write("SparticlesNoo");
        //std::cout<<"Worst: "<<worst<<std::endl;
        BOOST_REQUIRE(worst < 0.03);
}

/*BOOST_AUTO_TEST_CASE(dcpse_surface_sphere_proj) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n=512;
  size_t n_sp=n;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n,n,n};
  double grid_spacing{boxSize/(sz[0]-1)};
  double grid_spacing_surf=grid_spacing*30;
  double rCut{2.5 * grid_spacing_surf};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};

  constexpr int K = 1;
  // particles
  vector_dist_ws<3, double, aggregate<VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,double>> Sparticles(0, domain,bc,ghost);
  // 1. particles on the Spherical surface
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    //std::vector<Vector3f> data;
    //GenerateSphere(1,data);
    for(int i=1;i<n_sp;i++)
        {
            double y = 1.0 - (i /double(n_sp - 1.0)) * 2.0;
            double radius = sqrt(1 - y * y);
            double Golden_theta = Golden_angle * i;
            double x = cos(Golden_theta) * radius;
            double z = sin(Golden_theta) * radius;
            Sparticles.add();
            Sparticles.getLastPos()[0] = x;
            Sparticles.getLastPos()[1] = y;
            Sparticles.getLastPos()[2] = z;
            double rm=sqrt(x*x+y*y+z*z);
            Sparticles.getLastProp<2>()[0] = x/rm;
            Sparticles.getLastProp<2>()[1] = y/rm;
            Sparticles.getLastProp<2>()[2] = z/rm;
            Sparticles.getLastProp<4>()[0] = 1.0 ;
            Sparticles.getLastProp<4>()[1] = std::atan2(sqrt(x*x+y*y),z);
            Sparticles.getLastProp<4>()[2] = std::atan2(y,x);
            if(i<=2*(K)+1)
            {Sparticles.getLastSubset(1);}
            else
            {Sparticles.getLastSubset(0);}
        }
    //std::cout << "n: " << n << " - grid spacing: " << grid_spacing << " - rCut: " << rCut << "Surf Normal spacing" << grid_spacing<<std::endl;
  }

  Sparticles.map();
  Sparticles.ghost_get<3>();

  vector_dist_subset<3,double,aggregate<VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,double>> Sparticles_bulk(Sparticles,0);
  vector_dist_subset<3,double,aggregate<VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,VectorS<3,double>,double>> Sparticles_boundary(Sparticles,1);
  auto &bulkIds=Sparticles_bulk.getIds();
  auto &bdrIds=Sparticles_boundary.getIds();
  std::unordered_map<const lm,double,key_hash,key_equal> Alm;
  //Setting max mode l_max
  //Setting amplitudes to 1
  for(int l=0;l<=K;l++){
      for(int m=-l;m<=l;m++){
          Alm[std::make_tuple(l,m)]=0;
      }
  }
  Alm[std::make_tuple(1,0)]=1;
  auto it2 = Sparticles.getDomainIterator();
  while (it2.isNext()) {
      auto p = it2.get();
      Point<3, double> xP = Sparticles.getProp<4>(p);
      *//*double Sum=0;
      for(int m=-spL;m<=spL;++m)
      {
        Sum+=openfpm::math::Y(spL,m,xP[1],xP[2]);
      }*//*
      //Sparticles.getProp<ANADF>(p) = Sum;//openfpm::math::Y(K,K,xP[1],xP[2]);openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);;
      Sparticles.getProp<3>(p)[0]=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<3>(p)[1]=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<3>(p)[2]=openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<1>(p)[0]=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<1>(p)[1]=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      Sparticles.getProp<1>(p)[2]=-(K)*(K+1)*openfpm::math::sumY_Scalar<K>(xP[0],xP[1],xP[2],Alm);
      ++it2;
  }
  auto f=getV<3>(Sparticles);
  auto Df=getV<0>(Sparticles);

  SurfaceProjectedGradient<2> SGP{Sparticles,2,rCut,grid_spacing_surf};

  Sparticles.ghost_get<3>();
  Df=SGP(f);
  //Df=SLap(f);
  auto it3 = Sparticles.getDomainIterator();
  double worst = 0.0;
  while (it3.isNext()) {
      auto p = it3.get();
      //Sparticles.getProp<5>(p) = fabs(Sparticles.getProp<1>(p) - Sparticles.getProp<0>(p));
      if (fabs(Sparticles.getProp<1>(p)[0] - Sparticles.getProp<0>(p)[0]) > worst) {
          worst = fabs(Sparticles.getProp<1>(p)[0] - Sparticles.getProp<0>(p)[0]);
      }
      ++it3;
  }
        Sparticles.deleteGhost();
        //Sparticles.write("Sparticles");
        //std::cout<<worst;
        BOOST_REQUIRE(worst < 0.03);
}*/


 BOOST_AUTO_TEST_CASE(dcpse_surface_adaptive) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[3] = {81,81,1};
        Box<3, double> box({0, 0,-5}, {0.5, 0.5,5});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<3, double> ghost(spacing * 3.1);
        double rCut = 3.1  * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<3, double, aggregate<double,double,double,double,double,double,double[3]>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            domain.getLastPos()[2] = 0;

            ++it;
        }

        // Add multi res patch 1

        {
        const size_t sz2[3] = {40,40,1};
        Box<3,double> bx({0.25 + it.getSpacing(0)/4.0,0.25 + it.getSpacing(0)/4.0,-0.5},{sz2[0]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0, sz2[1]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0,0.5});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<3,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[1] = y;
            domain.getLastPos()[2] = 0;

            ++it2;
        }
        }

        // Add multi res patch 2

        {
        const size_t sz2[3] = {40,40,1};
        Box<3,double> bx({0.25 + 21.0*spacing/8.0,0.25 + 21.0*spacing/8.0,-5},{sz2[0]*spacing/4.0 + 0.25 + 21.0*spacing/8.0, sz2[1]*spacing/4.0 + 0.25 + 21*spacing/8.0,5});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<3,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();
            auto key = it2.get();
            double x = key.get(0) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[1] = y;
            domain.getLastPos()[2] = 0;

            ++it2;
        }
        }

        ///////////////////////

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        // Here fill me

        Box<3, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0,-5},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,5});

        Box<3, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,-5},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0,5});

        Box<3, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,-5},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0,5});

        Box<3, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,-5},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0,5});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        //vtk_box.write("vtk_box.vtk");

        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = domain.getPos(p);
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            }
            domain.getProp<6>(p)[0] = 0;
            domain.getProp<6>(p)[1] = 0;
            domain.getProp<6>(p)[2] = 1;

            ++it2;
        }

        domain.ghost_get<1,2,3>();
        SurfaceDerivative_xx<6> Dxx(domain, 2, rCut,3.9,support_options::ADAPTIVE);

/*        v=0;
        auto itNNN=domain.getDomainIterator();
        while(itNNN.isNext()){
            auto p=itNNN.get().getKey();
            Dxx.DrawKernel<0,decltype(domain)>(domain,p);
            domain.write_frame("Kernel",p);
            v=0;
            ++itNNN;
        }

*/
        //Dxx.DrawKernel<5,decltype(domain)>(domain,6161);
        //domain.write_frame("Kernel",6161);

        SurfaceDerivative_yy<6> Dyy(domain, 2, rCut,3.9,support_options::ADAPTIVE);
        SurfaceDerivative_zz<6> Dzz(domain, 2, rCut,3.9,support_options::ADAPTIVE);

        Dxx.save(domain,"Sdxx_test");
        Dyy.save(domain,"Sdyy_test");
        Dzz.save(domain,"Sdzz_test");

        domain.ghost_get<2>();
        sol=Dxx(anasol)+Dyy(anasol)+Dzz(anasol);
        domain.ghost_get<5>();

        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));

        }
        //std::cout << "Maximum Analytic Error: " << worst1 << std::endl;
        //domain.ghost_get<4>();
        //domain.write("Robin_anasol");
        BOOST_REQUIRE(worst1 < 0.03);

    }


    BOOST_AUTO_TEST_CASE(dcpse_surface_adaptive_load) {
//  int rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const size_t sz[3] = {81,81,1};
        Box<3, double> box({0, 0,-5}, {0.5, 0.5,5});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC,NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        Ghost<3, double> ghost(spacing * 3.1);
        double rCut = 3.1  * spacing;
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<3, double, aggregate<double,double,double,double,double,double,double[3]>> domain(0, box, bc, ghost);


        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        auto it = domain.getGridIterator(sz);
        while (it.isNext()) {
            domain.add();

            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            domain.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            domain.getLastPos()[1] = y;

            domain.getLastPos()[2] = 0;

            ++it;
        }

        // Add multi res patch 1

        {
        const size_t sz2[3] = {40,40,1};
        Box<3,double> bx({0.25 + it.getSpacing(0)/4.0,0.25 + it.getSpacing(0)/4.0,-0.5},{sz2[0]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0, sz2[1]*it.getSpacing(0)/2.0 + 0.25 + it.getSpacing(0)/4.0,0.5});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<3,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/2.0 + 0.25 + spacing/4.0;
            domain.getLastPos()[1] = y;
            domain.getLastPos()[2] = 0;

            ++it2;
        }
        }

        // Add multi res patch 2

        {
        const size_t sz2[3] = {40,40,1};
        Box<3,double> bx({0.25 + 21.0*spacing/8.0,0.25 + 21.0*spacing/8.0,-5},{sz2[0]*spacing/4.0 + 0.25 + 21.0*spacing/8.0, sz2[1]*spacing/4.0 + 0.25 + 21*spacing/8.0,5});
        openfpm::vector<size_t> rem;

        auto it = domain.getDomainIterator();

        while (it.isNext())
        {
        	auto k = it.get();

        	Point<3,double> xp = domain.getPos(k);

        	if (bx.isInside(xp) == true)
        	{
        		rem.add(k.getKey());
        	}

        	++it;
        }

        domain.remove(rem);

        auto it2 = domain.getGridIterator(sz2);
        while (it2.isNext()) {
            domain.add();

            auto key = it2.get();
            double x = key.get(0) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[0] = x;
            double y = key.get(1) * spacing/4.0 + 0.25 + 21*spacing/8.0;
            domain.getLastPos()[1] = y;
            domain.getLastPos()[2] = 0;

            ++it2;
        }
        }

        ///////////////////////

        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;
        openfpm::vector<aggregate<int>> ref_p;

        auto v = getV<0>(domain);
        auto RHS=getV<1>(domain);
        auto sol = getV<2>(domain);
        auto anasol = getV<3>(domain);
        auto err = getV<4>(domain);
        auto DCPSE_sol=getV<5>(domain);

        // Here fill me

        Box<3, double> up({box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0,-5},
                          {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0,5});

        Box<3, double> down({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0,-5},
                            {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0,5});

        Box<3, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,-5},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0,5});

        Box<3, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) + spacing / 2.0,-5},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0,5});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);

        // Create a writer and write
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        //vtk_box.write("vtk_box.vtk");


        auto it2 = domain.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = domain.getPos(p);
            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) = -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            } else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));

            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
                domain.getProp<1>(p) =  -2*M_PI*M_PI*sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
                domain.getProp<3>(p) = sin(M_PI*xp.get(0))*sin(M_PI*xp.get(1));
            }
            domain.getProp<6>(p)[0] = 0;
            domain.getProp<6>(p)[1] = 0;
            domain.getProp<6>(p)[2] = 1;

            ++it2;
        }

        domain.ghost_get<1,2,3>();
        SurfaceDerivative_xx<6> Dxx(domain, 2, rCut,3.9,support_options::LOAD);
        SurfaceDerivative_yy<6> Dyy(domain, 2, rCut,3.9,support_options::LOAD);
        SurfaceDerivative_zz<6> Dzz(domain, 2, rCut,3.9,support_options::LOAD);
        Dxx.load(domain,"Sdxx_test");
        Dyy.load(domain,"Sdyy_test");
        Dzz.load(domain,"Sdzz_test");

        domain.ghost_get<2>();
        sol=Dxx(anasol)+Dyy(anasol)+Dzz(anasol);
        domain.ghost_get<5>();

        double worst1 = 0.0;

        for(int j=0;j<bulk.size();j++)
        {   auto p=bulk.get<0>(j);
            if (fabs(domain.getProp<1>(p) - domain.getProp<2>(p)) >= worst1) {
                worst1 = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));
            }
            domain.getProp<4>(p) = fabs(domain.getProp<1>(p) - domain.getProp<2>(p));

        }
        //std::cout << "Maximum Analytic Error: " << worst1 << std::endl;

        //domain.ghost_get<4>();
        //domain.write("Robin_anasol");
        BOOST_REQUIRE(worst1 < 0.03);

    }

// when mass conservation holds on a deforming surface, it should change proportionally to the divergence of the velocity field. This test case checks the accuracy of the divergence operator on a growing sphere
BOOST_AUTO_TEST_CASE(dcpse_moving_surface_mass_conservation) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n_particles = 2048;
  size_t n_particles2 = 4096;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n_particles,n_particles,n_particles};
  double grid_spacing = std::sqrt(4.0*M_PI/n_particles);//{0.8/(std::pow(sz1[0],1.0/3.0)-1.0)};
  double grid_spacing2 = std::sqrt(4.0*M_PI/n_particles2);//{0.8/(std::pow(sz1[0],1.0/3.0)-1.0)};
  const size_t oporder = 5;
  double cutoff_factor = 3.5;
  double rCut{cutoff_factor * grid_spacing};
  double rCut2{cutoff_factor * grid_spacing2};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};
  Ghost<3,double> ghost2{rCut2 + grid_spacing2/8.0};
  // particles
  vector_dist<3,double, aggregate<double,double[3],double[3], double>> Sparticles(0,domain,bc,ghost);
  vector_dist<3,double, aggregate<double,double[3],double[3], double>> Sparticles2(0,domain,bc,ghost2);
  // properties: scalar_qty, normal, velocity, error
  // particles on the Spherical surface distributed with the Fibonacci sequence
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    // fill vector with resolution 1
    for(int i=0;i<n_particles;i++)
      {
	double y = 1.0 - (i /double(n_particles - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	Sparticles.add();
	Sparticles.getLastPos()[0] = x;
	Sparticles.getLastPos()[1] = y;
	Sparticles.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// fill unit surface normals
	Sparticles.getLastProp<1>()[0] = x/rm;
	Sparticles.getLastProp<1>()[1] = y/rm;
	Sparticles.getLastProp<1>()[2] = z/rm;
	// fill scalar field (constant concentration field)
	Sparticles.getLastProp<0>() = 1.0;
	// fill velocity field (unit normal)
	Sparticles.getLastProp<2>()[0] = x/rm;
	Sparticles.getLastProp<2>()[1] = y/rm;
	Sparticles.getLastProp<2>()[2] = z/rm;
      }
    // fill vector with resolution 2
    for(int i=0;i<n_particles2;i++)
      {
	double y = 1.0 - (i /double(n_particles2 - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	Sparticles2.add();
	Sparticles2.getLastPos()[0] = x;
	Sparticles2.getLastPos()[1] = y;
	Sparticles2.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// fill unit surface normals
	Sparticles2.getLastProp<1>()[0] = x/rm;
	Sparticles2.getLastProp<1>()[1] = y/rm;
	Sparticles2.getLastProp<1>()[2] = z/rm;
	// fill scalar field (constant concentration field)
	Sparticles2.getLastProp<0>() = 1.0;
	// fill velocity field (unit normal)
	Sparticles2.getLastProp<2>()[0] = x/rm;
	Sparticles2.getLastProp<2>()[1] = y/rm;
	Sparticles2.getLastProp<2>()[2] = z/rm;
      }
  }

  Sparticles.write("init1");
  Sparticles2.write("init2");

  Sparticles.map();
  Sparticles2.map();
  Sparticles.ghost_get<0,1,2>();
  Sparticles2.ghost_get<0,1,2>();

  // expressions for properties
  auto conc{getV<0>(Sparticles)};
  auto velo{getV<2>(Sparticles)};

  auto conc2{getV<0>(Sparticles2)};
  auto velo2{getV<2>(Sparticles2)};
  // time integration parameters
  double time = 0.0;
  double t_final = 0.4;//0.5*grid_spacing;
  double dt = 0.01*sqrt(4.0*M_PI/16384);
  //RK4 time integrator
  while(time < t_final)
  {
  	double updated_R1 = sqrt(Sparticles.getPos(0)[0]*Sparticles.getPos(0)[0] + Sparticles.getPos(0)[1]*Sparticles.getPos(0)[1] + Sparticles.getPos(0)[2]*Sparticles.getPos(0)[2]);
  	double updated_grid_spacing = sqrt(4.0*M_PI*updated_R1*updated_R1/n_particles);
  	double updated_R2 = sqrt(Sparticles2.getPos(0)[0]*Sparticles2.getPos(0)[0] + Sparticles2.getPos(0)[1]*Sparticles2.getPos(0)[1] + Sparticles2.getPos(0)[2]*Sparticles2.getPos(0)[2]);
  	double updated_grid_spacing2 = sqrt(4.0*M_PI*updated_R1*updated_R1/n_particles2);
	double updated_rCut = cutoff_factor*updated_grid_spacing;
	double updated_rCut2 = cutoff_factor*updated_grid_spacing2;
	SurfaceDerivative_x<1> Sdx{Sparticles,oporder,updated_rCut,updated_grid_spacing};
  	SurfaceDerivative_y<1> Sdy{Sparticles,oporder,updated_rCut,updated_grid_spacing};
 	SurfaceDerivative_z<1> Sdz{Sparticles,oporder,updated_rCut,updated_grid_spacing};
  	SurfaceDerivative_x<1> Sdx2{Sparticles2,oporder,updated_rCut2,updated_grid_spacing2};
  	SurfaceDerivative_y<1> Sdy2{Sparticles2,oporder,updated_rCut2,updated_grid_spacing2};
  	SurfaceDerivative_z<1> Sdz2{Sparticles2,oporder,updated_rCut2,updated_grid_spacing2};
  	auto k11 = -conc*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
  	auto k12 = -conc2*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
  	//auto k11 = -conc*2.0/updated_R1;
  	//auto k12 = -conc2*2.0/updated_R1;
	// integrate particle locations in time
	auto part = Sparticles.getDomainIterator();
	while(part.isNext())
	{
		auto a = part.get();
		for(int k = 0; k < 3; k++) Sparticles.getPos(a)[k] = Sparticles.getPos(a)[k] + 0.5*dt*Sparticles.template getProp<2>(a)[k];
		++part;
	}
	auto part2 = Sparticles2.getDomainIterator();
	while(part2.isNext())
	{
		auto a2 = part2.get();
		for (int k = 0; k < 3; k++) Sparticles2.getPos(a2)[k] = Sparticles2.getPos(a2)[k] + 0.5*dt*Sparticles2.template getProp<2>(a2)[k];
		++part2;
	}
  	Sparticles.map();
  	Sparticles.ghost_get<0,1,2>();
  	Sparticles2.map();
  	Sparticles2.ghost_get<0,1,2>();
	Sdx.update(Sparticles);
	Sdy.update(Sparticles);
	Sdz.update(Sparticles);
	Sdx2.update(Sparticles2);
	Sdy2.update(Sparticles2);
	Sdz2.update(Sparticles2);
	auto k21 = -(conc + 0.5*dt*k11)*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	auto k22 = -(conc2 + 0.5*dt*k12)*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
	auto k31 = -(conc + 0.5*dt*k21)*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	auto k32 = -(conc2 + 0.5*dt*k22)*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
  	//updated_R1 = sqrt(Sparticles.getPos(0)[0]*Sparticles.getPos(0)[0] + Sparticles.getPos(0)[1]*Sparticles.getPos(0)[1] + Sparticles.getPos(0)[2]*Sparticles.getPos(0)[2]);
	//auto k21 = -(conc + 0.5*dt*k11)*2.0/updated_R1;
	//auto k22 = -(conc2 + 0.5*dt*k12)*2.0/updated_R1;
	//auto k31 = -(conc + 0.5*dt*k21)*2.0/updated_R1;
	//auto k32 = -(conc2 + 0.5*dt*k22)*2.0/updated_R1;

	auto part12 = Sparticles.getDomainIterator();
	while(part12.isNext())
	{
		auto a = part12.get();
		for(int k = 0; k < 3; k++) Sparticles.getPos(a)[k] = Sparticles.getPos(a)[k] + 0.5*dt*Sparticles.template getProp<2>(a)[k];
		++part12;
	}
	auto part22 = Sparticles2.getDomainIterator();
	while(part22.isNext())
	{
		auto a2 = part22.get();
		for (int k = 0; k < 3; k++) Sparticles2.getPos(a2)[k] = Sparticles2.getPos(a2)[k] + 0.5*dt*Sparticles2.template getProp<2>(a2)[k];
		++part22;
	}
  	Sparticles.map();
  	Sparticles.ghost_get<0,1,2>();
  	Sparticles2.map();
  	Sparticles2.ghost_get<0,1,2>();
	Sdx.update(Sparticles);
	Sdy.update(Sparticles);
	Sdz.update(Sparticles);
	Sdx2.update(Sparticles2);
	Sdy2.update(Sparticles2);
	Sdz2.update(Sparticles2);
	auto k41 = -(conc + dt*k31)*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	auto k42 = -(conc2 + dt*k32)*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
  	//updated_R1 = sqrt(Sparticles.getPos(0)[0]*Sparticles.getPos(0)[0] + Sparticles.getPos(0)[1]*Sparticles.getPos(0)[1] + Sparticles.getPos(0)[2]*Sparticles.getPos(0)[2]);
	//auto k41 = -(conc + dt*k31)*2.0/updated_R1;
	//auto k42 = -(conc2 + dt*k32)*2.0/updated_R1;
	// integrate the concentration in time
	conc = conc + dt*(k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;
	conc2 = conc2 + dt*(k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
	//conc = conc - dt*conc*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	//conc2 = conc2 - dt*conc2*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
	Sdx.deallocate(Sparticles);
	Sdy.deallocate(Sparticles);
	Sdz.deallocate(Sparticles);
	Sdx2.deallocate(Sparticles2);
	Sdy2.deallocate(Sparticles2);
	Sdz2.deallocate(Sparticles2);
	time = time + dt;
	std::cout<<"time = "<<time<<std::endl;
  }
  // evaluate error at t_end
  auto it = Sparticles.getDomainIterator();
  double worst = 0.0;
  while (it.isNext()) 
  {
    auto p = it.get();

    double x = Sparticles.getPos(p)[0];
    double y = Sparticles.getPos(p)[1];
    double z = Sparticles.getPos(p)[2];

    //double theory =  2.0/sqrt(x*x + y*y + z*z);//1.0/(x*x + y*y + z*z);
    double theory = 1.0/(x*x + y*y + z*z);
    Sparticles.getProp<3>(p) = fabs(Sparticles.getProp<0>(p) - theory); // error

    if (fabs(Sparticles.getProp<3>(p)) > worst) {
      worst = fabs(Sparticles.getProp<3>(p));
    }
    ++it;
  }
  std::cout<<"having computed worst"<<std::endl;
  auto it2 = Sparticles2.getDomainIterator();
  double worst2 = 0.0;
  while (it2.isNext()) 
  {
    auto p = it2.get();

    double x = Sparticles2.getPos(p)[0];
    double y = Sparticles2.getPos(p)[1];
    double z = Sparticles2.getPos(p)[2];

    //double theory =  2.0/sqrt(x*x + y*y + z*z);//1.0/(x*x + y*y + z*z);
    double theory =  1.0/(x*x + y*y + z*z);
    Sparticles2.getProp<3>(p) = fabs(Sparticles2.getProp<0>(p) - theory); // error

    if (fabs(Sparticles2.getProp<3>(p)) > worst2) {
      worst2 = fabs(Sparticles2.getProp<3>(p));
    }
    ++it2;
  }
  Sparticles.write("after");
  Sparticles2.write("after2");
  std::cout<<"Linf mass conservation error with "<<n_particles<<" is: "<<worst<<std::endl;
  std::cout<<"Linf mass conservation error with "<<n_particles2<<" is: "<<worst2<<std::endl;
  std::cout<<"Convergence order is "<<std::log10(worst2/worst)/std::log10(grid_spacing2/grid_spacing)<<std::endl;
  std::cout<<"Having used an operator order = "<<oporder<<std::endl;
  Sparticles.deleteGhost();
  Sparticles.write("Sparticles_after");
  BOOST_REQUIRE(worst < 0.03);
}

// when mass conservation holds on a deforming surface, it should change proportionally to the divergence of the velocity field. This test case checks the accuracy of the divergence operator on a growing sphere
BOOST_AUTO_TEST_CASE(dcpse_moving_surface_mass_conservation) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n_particles = 2048;
  size_t n_particles2 = 4096;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz[3] = {n_particles,n_particles,n_particles};
  double grid_spacing = std::sqrt(4.0*M_PI/n_particles);//{0.8/(std::pow(sz1[0],1.0/3.0)-1.0)};
  double grid_spacing2 = std::sqrt(4.0*M_PI/n_particles2);//{0.8/(std::pow(sz1[0],1.0/3.0)-1.0)};
  const size_t oporder = 5;
  double cutoff_factor = 3.5;
  double rCut{cutoff_factor * grid_spacing};
  double rCut2{cutoff_factor * grid_spacing2};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut + grid_spacing/8.0};
  Ghost<3,double> ghost2{rCut2 + grid_spacing2/8.0};
  // particles
  vector_dist<3,double, aggregate<double,double[3],double[3], double>> Sparticles(0,domain,bc,ghost);
  vector_dist<3,double, aggregate<double,double[3],double[3], double>> Sparticles2(0,domain,bc,ghost2);
  // properties: scalar_qty, normal, velocity, error
  // particles on the Spherical surface distributed with the Fibonacci sequence
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    // fill vector with resolution 1
    for(int i=0;i<n_particles;i++)
      {
	double y = 1.0 - (i /double(n_particles - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	Sparticles.add();
	Sparticles.getLastPos()[0] = x;
	Sparticles.getLastPos()[1] = y;
	Sparticles.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// fill unit surface normals
	Sparticles.getLastProp<1>()[0] = x/rm;
	Sparticles.getLastProp<1>()[1] = y/rm;
	Sparticles.getLastProp<1>()[2] = z/rm;
	// fill scalar field (constant concentration field)
	Sparticles.getLastProp<0>() = 1.0;
	// fill velocity field (unit normal)
	Sparticles.getLastProp<2>()[0] = x/rm;
	Sparticles.getLastProp<2>()[1] = y/rm;
	Sparticles.getLastProp<2>()[2] = z/rm;
      }
    // fill vector with resolution 2
    for(int i=0;i<n_particles2;i++)
      {
	double y = 1.0 - (i /double(n_particles2 - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	Sparticles2.add();
	Sparticles2.getLastPos()[0] = x;
	Sparticles2.getLastPos()[1] = y;
	Sparticles2.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// fill unit surface normals
	Sparticles2.getLastProp<1>()[0] = x/rm;
	Sparticles2.getLastProp<1>()[1] = y/rm;
	Sparticles2.getLastProp<1>()[2] = z/rm;
	// fill scalar field (constant concentration field)
	Sparticles2.getLastProp<0>() = 1.0;
	// fill velocity field (unit normal)
	Sparticles2.getLastProp<2>()[0] = x/rm;
	Sparticles2.getLastProp<2>()[1] = y/rm;
	Sparticles2.getLastProp<2>()[2] = z/rm;
      }
  }

  Sparticles.write("init1");
  Sparticles2.write("init2");

  Sparticles.map();
  std::cout<<"after mapping particle set 1"<<std::endl;
  Sparticles2.map();
  std::cout<<"after mapping particle set 2"<<std::endl;
  Sparticles.ghost_get<0,1,2>();
  Sparticles2.ghost_get<0,1,2>();

  // surface dcpse parameters
  std::cout<<"Before initializing divergence operator on first particle set"<<std::endl;
  auto conc{getV<0>(Sparticles)};
  auto velo{getV<2>(Sparticles)};

  std::cout<<"Before initializing divergence operator on second particle set"<<std::endl;
  auto conc2{getV<0>(Sparticles2)};
  auto velo2{getV<2>(Sparticles2)};
  // time integration parameters
  double time = 0.0;
  double t_final = 0.4;//0.5*grid_spacing;
  double dt = 0.01*sqrt(4.0*M_PI/16384);
  while(time < t_final)
  {
  	double updated_R1 = sqrt(Sparticles.getPos(0)[0]*Sparticles.getPos(0)[0] + Sparticles.getPos(0)[1]*Sparticles.getPos(0)[1] + Sparticles.getPos(0)[2]*Sparticles.getPos(0)[2]);
  	double updated_grid_spacing = sqrt(4.0*M_PI*updated_R1*updated_R1/n_particles);
  	double updated_R2 = sqrt(Sparticles2.getPos(0)[0]*Sparticles2.getPos(0)[0] + Sparticles2.getPos(0)[1]*Sparticles2.getPos(0)[1] + Sparticles2.getPos(0)[2]*Sparticles2.getPos(0)[2]);
  	double updated_grid_spacing2 = sqrt(4.0*M_PI*updated_R1*updated_R1/n_particles2);
	double updated_rCut = cutoff_factor*updated_grid_spacing;
	double updated_rCut2 = cutoff_factor*updated_grid_spacing2;
	SurfaceDerivative_x<1> Sdx{Sparticles,oporder,updated_rCut,updated_grid_spacing};
  	SurfaceDerivative_y<1> Sdy{Sparticles,oporder,updated_rCut,updated_grid_spacing};
 	SurfaceDerivative_z<1> Sdz{Sparticles,oporder,updated_rCut,updated_grid_spacing};
  	SurfaceDerivative_x<1> Sdx2{Sparticles2,oporder,updated_rCut2,updated_grid_spacing2};
  	SurfaceDerivative_y<1> Sdy2{Sparticles2,oporder,updated_rCut2,updated_grid_spacing2};
  	SurfaceDerivative_z<1> Sdz2{Sparticles2,oporder,updated_rCut2,updated_grid_spacing2};
  	auto k11 = -conc*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
  	auto k12 = -conc2*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
  	//auto k11 = -conc*2.0/updated_R1;
  	//auto k12 = -conc2*2.0/updated_R1;
	//auto temppos1 = getV<PROP_POS>(Sparticles) + 0.5*dt*velo;
	//auto temppos2 = getV<PROP_POS>(Sparticles2) + 0.5*dt*velo2;
	// integrate particle locations in time
	auto part = Sparticles.getDomainIterator();
	while(part.isNext())
	{
		auto a = part.get();
		for(int k = 0; k < 3; k++) Sparticles.getPos(a)[k] = Sparticles.getPos(a)[k] + 0.5*dt*Sparticles.template getProp<2>(a)[k];
		++part;
	}
	auto part2 = Sparticles2.getDomainIterator();
	while(part2.isNext())
	{
		auto a2 = part2.get();
		for (int k = 0; k < 3; k++) Sparticles2.getPos(a2)[k] = Sparticles2.getPos(a2)[k] + 0.5*dt*Sparticles2.template getProp<2>(a2)[k];
		++part2;
	}
  	Sparticles.map();
  	Sparticles.ghost_get<0,1,2>();
  	Sparticles2.map();
  	Sparticles2.ghost_get<0,1,2>();
	Sdx.update(Sparticles);
	Sdy.update(Sparticles);
	Sdz.update(Sparticles);
	Sdx2.update(Sparticles2);
	Sdy2.update(Sparticles2);
	Sdz2.update(Sparticles2);
	auto k21 = -(conc + 0.5*dt*k11)*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	auto k22 = -(conc2 + 0.5*dt*k12)*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
	auto k31 = -(conc + 0.5*dt*k21)*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	auto k32 = -(conc2 + 0.5*dt*k22)*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
  	//updated_R1 = sqrt(Sparticles.getPos(0)[0]*Sparticles.getPos(0)[0] + Sparticles.getPos(0)[1]*Sparticles.getPos(0)[1] + Sparticles.getPos(0)[2]*Sparticles.getPos(0)[2]);
	//auto k21 = -(conc + 0.5*dt*k11)*2.0/updated_R1;
	//auto k22 = -(conc2 + 0.5*dt*k12)*2.0/updated_R1;
	//auto k31 = -(conc + 0.5*dt*k21)*2.0/updated_R1;
	//auto k32 = -(conc2 + 0.5*dt*k22)*2.0/updated_R1;

	auto part12 = Sparticles.getDomainIterator();
	while(part12.isNext())
	{
		auto a = part12.get();
		for(int k = 0; k < 3; k++) Sparticles.getPos(a)[k] = Sparticles.getPos(a)[k] + 0.5*dt*Sparticles.template getProp<2>(a)[k];
		++part12;
	}
	auto part22 = Sparticles2.getDomainIterator();
	while(part22.isNext())
	{
		auto a2 = part22.get();
		for (int k = 0; k < 3; k++) Sparticles2.getPos(a2)[k] = Sparticles2.getPos(a2)[k] + 0.5*dt*Sparticles2.template getProp<2>(a2)[k];
		++part22;
	}
  	Sparticles.map();
  	Sparticles.ghost_get<0,1,2>();
  	Sparticles2.map();
  	Sparticles2.ghost_get<0,1,2>();
	Sdx.update(Sparticles);
	Sdy.update(Sparticles);
	Sdz.update(Sparticles);
	Sdx2.update(Sparticles2);
	Sdy2.update(Sparticles2);
	Sdz2.update(Sparticles2);
	auto k41 = -(conc + dt*k31)*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	auto k42 = -(conc2 + dt*k32)*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
  	//updated_R1 = sqrt(Sparticles.getPos(0)[0]*Sparticles.getPos(0)[0] + Sparticles.getPos(0)[1]*Sparticles.getPos(0)[1] + Sparticles.getPos(0)[2]*Sparticles.getPos(0)[2]);
	//auto k41 = -(conc + dt*k31)*2.0/updated_R1;
	//auto k42 = -(conc2 + dt*k32)*2.0/updated_R1;
	// integrate the concentration in time
	conc = conc + dt*(k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;
	conc2 = conc2 + dt*(k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
	//conc = conc - dt*conc*(Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]));
	//conc2 = conc2 - dt*conc2*(Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]));
	Sdx.deallocate(Sparticles);
	Sdy.deallocate(Sparticles);
	Sdz.deallocate(Sparticles);
	Sdx2.deallocate(Sparticles2);
	Sdy2.deallocate(Sparticles2);
	Sdz2.deallocate(Sparticles2);
	time = time + dt;
	std::cout<<"time = "<<time<<std::endl;
  }
  //std::cout<<"Before computing divergence with first particle set"<<std::endl;
  //conc = Sdx(velo[0]) + Sdy(velo[1]) + Sdz(velo[2]);
  //std::cout<<"Before computing divergence with second particle set"<<std::endl;
  //conc2 = Sdx2(velo2[0]) + Sdy2(velo2[1]) + Sdz2(velo2[2]);
  auto it = Sparticles.getDomainIterator();
  double worst = 0.0;
  while (it.isNext()) 
  {
    auto p = it.get();

    double x = Sparticles.getPos(p)[0];
    double y = Sparticles.getPos(p)[1];
    double z = Sparticles.getPos(p)[2];

    //double theory =  2.0/sqrt(x*x + y*y + z*z);//1.0/(x*x + y*y + z*z);
    double theory = 1.0/(x*x + y*y + z*z);
    Sparticles.getProp<3>(p) = fabs(Sparticles.getProp<0>(p) - theory); // error

    if (fabs(Sparticles.getProp<3>(p)) > worst) {
      worst = fabs(Sparticles.getProp<3>(p));
    }
    ++it;
  }
  std::cout<<"having computed worst"<<std::endl;
  auto it2 = Sparticles2.getDomainIterator();
  double worst2 = 0.0;
  while (it2.isNext()) 
  {
    auto p = it2.get();

    double x = Sparticles2.getPos(p)[0];
    double y = Sparticles2.getPos(p)[1];
    double z = Sparticles2.getPos(p)[2];

    //double theory =  2.0/sqrt(x*x + y*y + z*z);//1.0/(x*x + y*y + z*z);
    double theory =  1.0/(x*x + y*y + z*z);
    Sparticles2.getProp<3>(p) = fabs(Sparticles2.getProp<0>(p) - theory); // error

    if (fabs(Sparticles2.getProp<3>(p)) > worst2) {
      worst2 = fabs(Sparticles2.getProp<3>(p));
    }
    ++it2;
  }
  Sparticles.write("after");
  Sparticles2.write("after2");
  std::cout<<"Linf Laplace error with roughly "<<n_particles<<" is: "<<worst<<std::endl;
  std::cout<<"Linf Laplace error with roughly "<<n_particles2<<" is: "<<worst2<<std::endl;
  std::cout<<"Convergence order is "<<std::log10(worst2/worst)/std::log10(grid_spacing2/grid_spacing)<<std::endl;
  std::cout<<"Operator order = "<<oporder<<std::endl;
  Sparticles.deleteGhost();
  Sparticles.write("Sparticles_after");
  BOOST_REQUIRE(worst < 0.03);
}

BOOST_AUTO_TEST_CASE(dcpse_surface_p2p_interpolation_sphere_scalar) {
  auto & v_cl = create_vcluster();
  timer tt;
  tt.start();
  size_t n_from1=4096;
  size_t n_from2=8192;
  size_t n_to=256;
  // Domain
  double boxP1{-1.5}, boxP2{1.5};
  double boxSize{boxP2 - boxP1};
  size_t sz1[3] = {n_from1,n_from1,n_from1};
  size_t sz2[3] = {n_from2,n_from2,n_from2};
  size_t szTo[3] = {n_to,n_to,n_to};
  double grid_spacing1 = std::sqrt(4.0*M_PI/n_from1);//{0.8/(std::pow(sz1[0],1.0/3.0)-1.0)};
  double grid_spacing2 = std::sqrt(4.0*M_PI/n_from2);//{0.8/(std::pow(sz2[0],1.0/3.0)-1.0)};
  double grid_spacingTo = std::sqrt(4.0*M_PI/n_to);//{0.8/(std::pow(szTo[0],1.0/3.0)-1.0)};
  double grid_spacing_surf2=grid_spacing2;
  double grid_spacing_surf1=grid_spacing1;
  double grid_spacing_surfTo=grid_spacingTo;
  double cutoff_factor = 3.5;
  double rCut2{cutoff_factor * grid_spacing_surf2};
  double rCut1{cutoff_factor * grid_spacing_surf1};
  double rCutTo{cutoff_factor * grid_spacing_surfTo};

  Box<3,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost1{rCut1 + grid_spacing1/8.0};
  Ghost<3,double> ghost2{rCut2 + grid_spacing2/8.0};
  Ghost<3,double> ghostTo{rCutTo + grid_spacingTo/8.0};
  // particles
  vector_dist<3,double, aggregate<double,double[3]>> SparticlesFrom1(0,domain,bc,ghost1);
  vector_dist<3,double, aggregate<double,double[3]>> SparticlesFrom2(0,domain,bc,ghost2);
  // properties: scalar_qty, normal, error
  vector_dist<3,double, aggregate<double,double,double,double>> SparticlesTo(0,domain,bc,ghostTo);
  // properties: scalar obtained from interpolation from data with resolution 1,
  //		scalar obtained from interpolation from data with resolution 2,
  //		error of scalar 1, error of scalar 2
  // particles on the Spherical surface distributed with the Fibonacci sequence
  double Golden_angle=M_PI * (3.0 - sqrt(5.0));
  if (v_cl.rank() == 0) {
    // fill vector with resolution 1
    for(int i=0;i<n_from1;i++)
      {
	double y = 1.0 - (i /double(n_from1 - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	SparticlesFrom1.add();
	SparticlesFrom1.getLastPos()[0] = x;
	SparticlesFrom1.getLastPos()[1] = y;
	SparticlesFrom1.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// fill unit surface normals
	SparticlesFrom1.getLastProp<1>()[0] = x/rm;
	SparticlesFrom1.getLastProp<1>()[1] = y/rm;
	SparticlesFrom1.getLastProp<1>()[2] = z/rm;
	// fill scalar field (spherical harmonic 2,0)
	//SparticlesFrom1.getLastProp<0>() = std::sqrt(5.0/(16.0*M_PI)) * (3*z*z - 1.0);
        // spherical harmonic 2,1
	//SparticlesFrom1.getLastProp<0>() = 0.5*std::sqrt(15.0/M_PI)*x*z;
	// spherical harmonic 2,2
	//SparticlesFrom1.getLastProp<0>() = 0.25*std::sqrt(15.0/M_PI)*(x*x - y*y);
	// spherical harmonic 3,2
	SparticlesFrom1.getLastProp<0>() = 0.25*std::sqrt(105.0/M_PI)*(x*x - y*y)*z;
	// spherical harmonic 0,0
	//SparticlesFrom1.getLastProp<0>() = 0.5/std::sqrt(M_PI);
      }
    // fill vector with resolution 2
    for(int i=0;i<n_from2;i++)
      {
	double y = 1.0 - (i /double(n_from2 - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	SparticlesFrom2.add();
	SparticlesFrom2.getLastPos()[0] = x;
	SparticlesFrom2.getLastPos()[1] = y;
	SparticlesFrom2.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// fill unit surface normals
	SparticlesFrom2.getLastProp<1>()[0] = x/rm;
	SparticlesFrom2.getLastProp<1>()[1] = y/rm;
	SparticlesFrom2.getLastProp<1>()[2] = z/rm;
	// fill scalar field (spherical harmonic)
	//SparticlesFrom2.getLastProp<0>() = std::sqrt(5.0/(16.0*M_PI)) * (3*z*z - 1.0);
        // spherical harmonic 2,1
	//SparticlesFrom2.getLastProp<0>() = 0.5*std::sqrt(15.0/M_PI)*x*z;
      	// spherical harmonic 2,2
	//SparticlesFrom2.getLastProp<0>() = 0.25*std::sqrt(15.0/M_PI)*(x*x - y*y);
	// spherical harmonic 3,2
	SparticlesFrom2.getLastProp<0>() = 0.25*std::sqrt(105.0/M_PI)*(x*x - y*y)*z;
	// spherical harmonic 0,0
      	//SparticlesFrom2.getLastProp<0>() = 0.5/std::sqrt(M_PI);
      }
      // fill vector with positions at which surface interpolation is supposed to be performed
    	for(int i=0;i<((int)(n_to-1));i++)
    {
	double y = 1.0 - (i /double(n_to - 1.0)) * 2.0;
	double radius = sqrt(1 - y * y);
	double Golden_theta = Golden_angle * i;
	double x = cos(Golden_theta) * radius;
	double z = sin(Golden_theta) * radius;
	SparticlesTo.add();
	SparticlesTo.getLastPos()[0] = x;
	SparticlesTo.getLastPos()[1] = y;
	SparticlesTo.getLastPos()[2] = z;
	double rm=sqrt(x*x+y*y+z*z);
	// initialize scalar fields as 0.0
	SparticlesTo.getLastProp<0>() = 0.0;//std::sqrt(3.0/(4.0*M_PI)) * z;
	SparticlesTo.getLastProp<1>() = 0.0;//std::sqrt(3.0/(4.0*M_PI)) * z;
      }
  }

  SparticlesFrom1.write("from1_before");
  SparticlesFrom2.write("from2_before");
  SparticlesTo.write("to_before");

  SparticlesFrom1.map();
  SparticlesFrom1.ghost_get<0,1>();
  SparticlesFrom2.map();
  SparticlesFrom2.ghost_get<0,1>();
  SparticlesTo.map();
  SparticlesTo.ghost_get<0>();

  const size_t oporder = 5;

  PPInterpolation<decltype(SparticlesFrom1),decltype(SparticlesTo), 1> ppSurface(SparticlesFrom1,SparticlesTo,oporder,rCut1,grid_spacing_surf1, support_options::RADIUS);
  ppSurface.p2p<0,0>();
  PPInterpolation<decltype(SparticlesFrom2),decltype(SparticlesTo), 1> ppSurface2(SparticlesFrom2,SparticlesTo,oporder,rCut2,grid_spacing_surf2, support_options::RADIUS);
  ppSurface2.p2p<0,1>();

  auto it = SparticlesTo.getDomainIterator();
  double worst = 0.0;
  double worst2 = 0.0;
  while (it.isNext()) {
    auto p = it.get();

    double x = SparticlesTo.getPos(p)[0];
    double y = SparticlesTo.getPos(p)[1];
    double z = SparticlesTo.getPos(p)[2];

    //double sphericalHarmonic = std::sqrt(5.0/(16.0*M_PI)) * (3*z*z - 1.0);
    // spherical harmonic 2,1
    // double sphericalHarmonic = 0.5*std::sqrt(15.0/M_PI)*x*z;
    // spherical harmonic 2,2
    //double sphericalHarmonic = 0.25*std::sqrt(15.0/M_PI)*(x*x - y*y);
    // spherical harmonic 3,2
    double sphericalHarmonic = 0.25*std::sqrt(105.0/M_PI)*(x*x - y*y)*z;
    //double sphericalHarmonic = 0.5/std::sqrt(M_PI);

    SparticlesTo.getProp<2>(p) = fabs(SparticlesTo.getProp<0>(p) - sphericalHarmonic); // error
    SparticlesTo.getProp<3>(p) = fabs(SparticlesTo.getProp<1>(p) - sphericalHarmonic); // error

    if (fabs(SparticlesTo.getProp<0>(p) - sphericalHarmonic) > worst) {
      worst = fabs(SparticlesTo.getProp<0>(p) - sphericalHarmonic);
    }
    if (fabs(SparticlesTo.getProp<1>(p) - sphericalHarmonic) > worst2) {
      worst2 = fabs(SparticlesTo.getProp<1>(p) - sphericalHarmonic);
    }
    ++it;
  }
  std::cout<<"Linf interpolation error with h_from1=1/"<<n_from1<<" to h_to=1/"<<n_to<<" is: "<<worst<<std::endl;
  std::cout<<"Linf interpolation error with h_from2=1/"<<n_from2<<" to h_to=1/"<<n_to<<" is: "<<worst2<<std::endl;
  std::cout<<"Convergence order is "<<std::log10(worst2/worst)/std::log10(std::sqrt((float)n_from1/(float)n_from2))<<std::endl;
  std::cout<<"Operator order = "<<oporder<<std::endl;
  SparticlesTo.deleteGhost();
  SparticlesFrom1.deleteGhost();
  SparticlesTo.write("SparticlesTo_after");
  SparticlesFrom1.write("SparticlesFrom1_after");
  BOOST_REQUIRE(worst < 0.03);
  BOOST_REQUIRE(worst2 < 0.03);
}

BOOST_AUTO_TEST_CASE(dcpse_surface_p2p_interpolation_plane_scalar) {

  auto & v_cl = create_vcluster();

  size_t n{10};
  double rCut{3.1};
  
  // Domain and simulation parameters
  Box<3,double> domain{{-1,-1,-1},{1,1,1}};
  size_t sz[3] = {n,n,n};
  double grid_spacing{2.0/(n-1)};
  size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};
  Ghost<3,double> ghost{rCut};

  vector_dist<3,double,aggregate<double,double[3],double>> part_from{0,domain,bc,ghost};
  vector_dist<3,double,aggregate<double,double[3],double>> part_to{0,domain,bc,ghost};
  // props: scalar_qty, normal, error

  // Create particles_from in a grid-like manner
  if (v_cl.rank() == 0) {

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {

	part_from.add();
	
	part_from.getLastPos()[0] = -1 + i*grid_spacing;
	part_from.getLastPos()[1] = -1 + j*grid_spacing;
	part_from.getLastPos()[2] = 0;
	
	part_from.getLastProp<0>() = std::fabs(part_from.getLastPos()[0]); // scalar_qty
	part_from.getLastProp<1>()[0] = 0; // normal_x
	part_from.getLastProp<1>()[1] = 0; // normal_y
	part_from.getLastProp<1>()[2] = 1; // normal_z
	part_from.getLastProp<2>() = 0; // error
      }
  }
  part_from.map();
  part_from.ghost_get<0,1>();

  // Create particles_to in a grid-like manner + random noise
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist_threshold{-1.0,1.0};
  
  if (v_cl.rank() == 0) {

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {

	part_to.add();
	
	part_to.getLastPos()[0] = -1 + i*grid_spacing + dist_threshold(gen) * 0.07; // 7% noise
	part_to.getLastPos()[1] = -1 + j*grid_spacing + dist_threshold(gen) * 0.07;
	part_to.getLastPos()[2] = 0;
	
	part_to.getLastProp<0>() = 0; // scalar_qty
	part_to.getLastProp<1>()[0] = 0; // normal_x
	part_to.getLastProp<1>()[1] = 0; // normal_y
	part_to.getLastProp<1>()[2] = 1; // normal_z
	part_to.getLastProp<2>() = 0; // error
      }
  }
  part_to.map();
  part_to.ghost_get<0,1>();

  // Interpolate
  PPInterpolation<decltype(part_from),decltype(part_to),1> ppSurface(part_from,part_to,2,rCut,grid_spacing);
  ppSurface.p2p<0,0>();

  // Compute maximum error
  auto it = part_to.getDomainIterator();
  double worst = 0.0;
  while (it.isNext()) {
    auto key{it.get()};

    part_to.getProp<2>(key) = std::fabs(part_to.getProp<0>(key) - part_to.getPos(key)[0]); // error

    if (part_to.getProp<2>(key) > worst)
      worst = part_to.getProp<2>(key);
    ++it;
  }

  // Write particles
  part_from.deleteGhost();
  part_from.write("surface_p2p_interp_plane_part_from");
  part_to.deleteGhost();
  part_to.write("surface_p2p_interp_plane_part_to");

  BOOST_REQUIRE(worst < 0.03);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
