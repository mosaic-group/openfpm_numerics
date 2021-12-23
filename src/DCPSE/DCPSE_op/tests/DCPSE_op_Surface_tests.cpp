//
// Created by Abhinav Singh on 15.11.21.
//

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
BOOST_AUTO_TEST_CASE(dcpse_surface_sphere) {
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
  Sdxx.DrawKernel<5>(Sparticles,0);
  //Sdyy.DrawKernel<5>(Sparticles,0);
  //Sdzz.DrawKernel<5>(Sparticles,0);
/*  std::cout<<"SDXX:"<<std::endl;
  Sdxx.checkMomenta(Sparticles);
  std::cout<<"SDYY:"<<std::endl;
  Sdyy.checkMomenta(Sparticles);
  std::cout<<"SDZZ:"<<std::endl;
  Sdzz.checkMomenta(Sparticles);*/


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
        //std::cout<<worst;
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
  Sdxx.DrawKernel<5>(Sparticles,0);
  //Sdyy.DrawKernel<5>(Sparticles,0);
  //Sdzz.DrawKernel<5>(Sparticles,0);
/*  std::cout<<"SDXX:"<<std::endl;
  Sdxx.checkMomenta(Sparticles);
  std::cout<<"SDYY:"<<std::endl;
  Sdyy.checkMomenta(Sparticles);
  std::cout<<"SDZZ:"<<std::endl;
  Sdzz.checkMomenta(Sparticles);*/


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
BOOST_AUTO_TEST_SUITE_END()

#endif
#endif