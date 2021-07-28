//
// Created by Abhinav Singh on 01.12.20.
//
#include "config.h"
#include <type_traits>
#include <cstring>
#include "util/common.hpp"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "OdeIntegrators/boost_vector_algebra_ofp.hpp"

typedef texp_v<double> state_type;
const double a = 2.8e-4;
const double b = 5e-3;
const double tau = .1;
const double k = .005;

void *vectorGlobal;
typedef vector_dist<2, double, aggregate<double,double,double,double,double,double>> vector_type;

template<typename op_type>
struct Fitz
{
    double a = 2.8e-4;
    double b = 5e-3;
    double tau = .1;
    double k = 0;
    op_type &Lap;


    //Constructor
    Fitz(op_type &Lap):Lap(Lap)
    {}

    void operator()( const state_type_2d_ofp &x , state_type_2d_ofp &dxdt , const double t ) const
    {
        vector_type &Temp= *(vector_type *) vectorGlobal;
        auto u=getV<4>(Temp);
        auto v=getV<5>(Temp);
        u=x.data.get<0>();
        v=x.data.get<1>();

        Temp.ghost_get<4,5>();
        dxdt.data.get<0>() = a*Lap(u) + u - u*v - v + k;
        dxdt.data.get<1>() = (b*Lap(v)+ u-v) / tau;
    }
};

void Exponential_struct_ofp( const state_type_2d_ofp &x , state_type_2d_ofp &dxdt , const double t )
{
    dxdt.data.get<0>() = x.data.get<0>();
    dxdt.data.get<1>() = 2.0*x.data.get<1>();
}
void Exponential_struct_ofp2( const state_type_3d_ofp &x , state_type_3d_ofp &dxdt , const double t )
{
    dxdt.data.get<0>() = x.data.get<0>();
    dxdt.data.get<1>() = 2.0*x.data.get<1>();
}


void Exponential( const state_type &x , state_type &dxdt , const double t )
{
    dxdt = x;
}
void sigmoid( const state_type &x , state_type &dxdt , const double t )
{
    dxdt = x*(1-x);
}

BOOST_AUTO_TEST_SUITE(odeInt_BASE_tests)
BOOST_AUTO_TEST_CASE(odeint_base_test1) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
        Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
        size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            Particles.getLastProp<0>() = xp0*yp0*exp(0);
            Particles.getLastProp<1>() = xp0*yp0*exp(0.4);
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();
        auto Init = getV<0>(Particles);
        auto Sol = getV<1>(Particles);
        auto OdeSol = getV<2>(Particles);

        state_type x0;
        x0=Init;
        // The rhs of x' = f(x)

        double t=0,tf=0.4;
        const double dt=0.1;

        //This doesnt work Why?
        //size_t steps=boost::numeric::odeint::integrate(Exponential,x0,0.0,tf,dt);

        size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< state_type >(),Exponential,x0,0.0,tf,dt);

        OdeSol=x0;
        auto it2 = Particles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst) {
                worst = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it2;
        }

        //std::cout<<worst<<std::endl;
        BOOST_REQUIRE(worst < 1e-6);

        x0=Init;
        boost::numeric::odeint::runge_kutta4< state_type > rk4;
        while (t<tf)
        {
            rk4.do_step(Exponential,x0,t,dt);
            OdeSol=x0;
            Particles.write_frame("OdeSol",int(t/dt));
            t+=dt;
        }

        OdeSol=x0;
        auto it3 = Particles.getDomainIterator();
        double worst2 = 0.0;
        while (it3.isNext()) {
            auto p = it3.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst2) {
                worst2 = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it3;
        }
        //std::cout<<worst2<<std::endl;
        BOOST_REQUIRE(worst < 1e-6);
        BOOST_REQUIRE_EQUAL(worst,worst2);
}

    BOOST_AUTO_TEST_CASE(odeint_base_test_STRUCT_ofp) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
        Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
        size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double,double,double,double,double,double>> Particles(0, box, bc, ghost);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            Particles.getLastProp<0>() = xp0*yp0*exp(0);
            Particles.getLastProp<1>() = xp0*yp0*exp(0.4);
            Particles.getLastProp<2>() = xp0*yp0*exp(0);
            Particles.getLastProp<3>() = xp0*yp0*exp(0.8);
            ++it;
        }
        Particles.map();
        Particles.ghost_get<0>();
        auto Init1 = getV<0>(Particles);
        auto Sol1 = getV<1>(Particles);
        auto Init2 = getV<2>(Particles);
        auto Sol2 = getV<3>(Particles);
        auto OdeSol1 = getV<4>(Particles);
        auto OdeSol2 = getV<5>(Particles);

        state_type_3d_ofp x0;//(Init1,Init2);
        x0.data.get<0>()=Init1;
        x0.data.get<1>()=Init2;
        x0.data.get<2>()=Init1;

        // The rhs of x' = f(x)
        double t=0,tf=0.4;
        const double dt=0.1;

        //size_t steps=boost::numeric::odeint::integrate(Exponential_struct,x0,0.0,tf,dt);

        //size_t steps=boost::numeric::odeint::integrate_const(
        //        boost::numeric::odeint::runge_kutta4< state_type_struct_ofp<vector_type>,double,state_type_struct_ofp<vector_type>,double,boost::numeric::odeint::vector_space_algebra_ofp > ()
        //                ,Exponential_struct_ofp,x0,0.0,tf,dt);

        //typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type_struct_ofp<vector_type>,double,state_type_struct_ofp<vector_type>,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;
        typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type_3d_ofp,double,state_type_3d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;
        integrate_adaptive( stepper_type() , Exponential_struct_ofp2 , x0 , t , tf , dt);

        OdeSol1=x0.data.get<0>();
        OdeSol2=x0.data.get<1>();
        auto it2 = Particles.getDomainIterator();
        double worst = 0.0;
        double worst2 = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<4>(p)) > worst) {
                worst = fabs(Particles.getProp<1>(p) - Particles.getProp<4>(p));
            }
            if (fabs(Particles.getProp<3>(p) - Particles.getProp<5>(p)) > worst2) {
                worst2 = fabs(Particles.getProp<3>(p) - Particles.getProp<5>(p));
            }
            ++it2;
        }




        std::cout<<worst<<std::endl;
        std::cout<<worst2<<std::endl;
        //BOOST_REQUIRE(worst < 1e-6);

/*        x0.u=Init1;
        x0.v=Init2;

        //x0=Init;
        boost::numeric::odeint::runge_kutta4< state_type_struct_ofp > rk4;
        while (t<tf)
        {
            rk4.do_step(Exponential_struct,x0,t,dt);
            OdeSol1=x0.u;
            OdeSol2=x0.v;
            Particles.write_frame("OdeSol",int(t/dt));
            t+=dt;
        }

        OdeSol1=x0.u;
        OdeSol2=x0.v;
        auto it3 = Particles.getDomainIterator();
        double worst3 = 0.0;
        double worst4 = 0.0;
        while (it3.isNext()) {
            auto p = it3.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<4>(p)) > worst3) {
                worst3 = fabs(Particles.getProp<1>(p) - Particles.getProp<4>(p));
            }
            if (fabs(Particles.getProp<3>(p) - Particles.getProp<5>(p)) > worst4) {
                worst4 = fabs(Particles.getProp<3>(p) - Particles.getProp<5>(p));
            }
            ++it3;
        }
        std::cout<<worst3<<std::endl;
        std::cout<<worst4<<std::endl;*/
        //BOOST_REQUIRE(worst < 1e-6);
        //BOOST_REQUIRE_EQUAL(worst,worst2);
    }





    BOOST_AUTO_TEST_CASE(odeint_base_test2) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
        Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
        size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

        double t=0.0,tf=0.5;
        const double dt=0.1;

        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            Particles.getLastProp<0>() = xp0*yp0*exp(t);
            Particles.getLastProp<1>() = xp0*yp0*exp(tf);
            ++it;
        }
        Particles.map();
        Particles.ghost_get<0>();
        auto Init = getV<0>(Particles);
        auto Sol = getV<1>(Particles);
        auto OdeSol = getV<2>(Particles);

        state_type x0;
        x0=Init;
        // The rhs of x' = f(x)




        //This doesnt work Why?
        //size_t steps=boost::numeric::odeint::integrate(Exponential,x0,0.0,tf,dt);
        typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type > > stepper_type;
        integrate_adaptive( stepper_type() , Exponential , x0 , t , tf , dt);
        //size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< state_type >(),Exponential,x0,0.0,tf,dt);

        OdeSol=x0;
        auto it2 = Particles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst) {
                worst = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it2;
        }

        //std::cout<<worst<<std::endl;
        BOOST_REQUIRE(worst < 1e-6);
/*
        x0=Init;
        boost::numeric::odeint::runge_kutta4< state_type > rk4;
        for( size_t i=0 ; i<int(tf/dt) ; ++i,t+=dt )
        {
            rk4.do_step(Exponential,x0,t,dt);
            OdeSol=x0;
            Particles.write_frame("OdeSol",int(t/dt));
            t+=dt;
        }

        OdeSol=x0;
        auto it3 = Particles.getDomainIterator();
        double worst2 = 0.0;
        while (it3.isNext()) {
            auto p = it3.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst2) {
                worst2 = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it3;
        }
        //std::cout<<worst2<<std::endl;
        BOOST_REQUIRE(worst < 1e-6);
        BOOST_REQUIRE_EQUAL(worst,worst2);*/
    }
    BOOST_AUTO_TEST_CASE(odeint_base_test3) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
        Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
        size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

        double t=0.0,tf=0.5;
        const double dt=0.1;

        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            Particles.getLastProp<0>() = xp0*yp0*1/(1+exp(-t));
            Particles.getLastProp<1>() = xp0*yp0*1/(1+exp(-tf));
            ++it;
        }
        Particles.map();
        Particles.ghost_get<0>();
        auto Init = getV<0>(Particles);
        auto Sol = getV<1>(Particles);
        auto OdeSol = getV<2>(Particles);

        state_type x0;
        x0=Init;
        // The rhs of x' = f(x)
        //size_t steps=boost::numeric::odeint::integrate(Exponential,x0,0.0,tf,dt);
        //typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type_struct_ofp2<vector_type>,double,state_type_struct_ofp2<vector_type>,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;
        //integrate_adaptive( stepper_type() , sigmoid , x0 , t , tf , dt);
        //size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< state_type >(),Exponential,x0,0.0,tf,dt);

        OdeSol=x0;
        auto it2 = Particles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst) {
                worst = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it2;
        }

        std::cout<<worst<<std::endl;
        Particles.write("SigmoidODE");

        //BOOST_REQUIRE(worst < 1e-6);
/*
        x0=Init;
        boost::numeric::odeint::runge_kutta4< state_type > rk4;
        for( size_t i=0 ; i<int(tf/dt) ; ++i,t+=dt )
        {
            rk4.do_step(Exponential,x0,t,dt);
            OdeSol=x0;
            Particles.write_frame("OdeSol",int(t/dt));
            t+=dt;
        }

        OdeSol=x0;
        auto it3 = Particles.getDomainIterator();
        double worst2 = 0.0;
        while (it3.isNext()) {
            auto p = it3.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst2) {
                worst2 = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it3;
        }
        //std::cout<<worst2<<std::endl;
        BOOST_REQUIRE(worst < 1e-6);
        BOOST_REQUIRE_EQUAL(worst,worst2);*/
    }

    BOOST_AUTO_TEST_CASE(odeint_diffusion) {
        size_t edgeSemiSize = 40;
        const size_t sz[2] = {edgeSemiSize+1,edgeSemiSize+1};
        Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
        size_t bc[2] = { PERIODIC, PERIODIC };
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

        double t=0.0,tf=0.5;
        const double dt=0.1;

        auto it = Particles.getGridIterator(sz);
        while (it.isNext())
        {
            Particles.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double xp0 = k0 * spacing[0];
            Particles.getLastPos()[0] = xp0;
            mem_id k1 = key.get(1);
            double yp0 = k1 * spacing[1];
            Particles.getLastPos()[1] = yp0;
            Particles.getLastProp<0>() = 0.0;
            Particles.getLastProp<1>() = 0.0;
            if (xp0==0.5 && yp0==0.5)
            {Particles.getLastProp<0>() = 1.0;}
            ++it;
        }
        Particles.map();
        Particles.ghost_get<0>();
        auto Init = getV<0>(Particles);
        auto Sol = getV<1>(Particles);
        auto OdeSol = getV<2>(Particles);

        state_type x0;
        x0=Init;
        // The rhs of x' = f(x)
        //size_t steps=boost::numeric::odeint::integrate(Exponential,x0,0.0,tf,dt);
        typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type > > stepper_type;
        integrate_adaptive( stepper_type() , sigmoid , x0 , t , tf , dt);
        //size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< state_type >(),Exponential,x0,0.0,tf,dt);

        OdeSol=x0;
        auto it2 = Particles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst) {
                worst = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it2;
        }

        std::cout<<worst<<std::endl;
        Particles.write("OdeInt-Diffusion");
    }

BOOST_AUTO_TEST_CASE(dcpse_op_react_diff_test) {
        size_t edgeSemiSize = 20;
        const size_t sz[2] = {2 * edgeSemiSize+1, 2 * edgeSemiSize+1};
        Box<2, double> box({0, 0}, {1.0, 1.0});
        size_t bc[2] = {PERIODIC, PERIODIC};
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        Ghost<2, double> ghost(spacing[0] * 3);
        double rCut = 3.0 * spacing[0];
        BOOST_TEST_MESSAGE("Init vector_dist...");
        double sigma2 = spacing[0] * spacing[1] / (2 * 4);

        //properties: u, v, du, dv
        vector_dist<2, double, aggregate<double, double, double, double,double,double>> domain(0, box,
        bc,
        ghost);

        //Init_DCPSE(domain)
        BOOST_TEST_MESSAGE("Init domain...");

        std::srand(std::time(nullptr));
        const double max_rand = 2147483647;

        auto it = domain.getGridIterator(sz);
        size_t pointId = 0;
        size_t counter = 0;
        double minNormOne = 999;
        while (it.isNext()) {
            domain.add();
            auto key = it.get();
            mem_id k0 = key.get(0);
            double x = k0 * spacing[0];
            domain.getLastPos()[0] = x;//+ gaussian(rng);
            mem_id k1 = key.get(1);
            double y = k1 * spacing[1];
            domain.getLastPos()[1] = y;//+gaussian(rng);
            // Here fill the function value
            domain.template getLastProp<0>() = 0.0;
            domain.template getLastProp<1>() = 0.0;
            domain.template getLastProp<2>() = 0;
            domain.template getLastProp<3>() = 0;
            if (x==0.5 && y==0.5){
                domain.template getLastProp<0>() = 1.0;
                domain.template getLastProp<1>() = 1.0;
            }
            ++counter;
            ++it;
        }
        BOOST_TEST_MESSAGE("Sync domain across processors...");

        domain.map();
        domain.ghost_get<0>();

        domain.write_frame("ReactDiffTest",0);

        //Derivative_x Dx(domain, 2, rCut);
        //Derivative_y Dy(domain, 2, rCut);
        //Gradient Grad(domain, 2, rCut);
        vectorGlobal=(void *) &domain;

        Laplacian Lap(domain, 2, rCut);

        auto u = getV<0>(domain);
        auto v = getV<1>(domain);
        auto fu = getV<2>(domain);
        auto fv = getV<3>(domain);

        Fitz<Laplacian> System(Lap);
        state_type_2d_ofp x0;
        x0.data.get<0>()=u;
        x0.data.get<1>()=v;

        double dt = 0.001;
        double t = 0.0;
        double tf = 0.5;
        domain.write("ReactDiffTestInit");
        typedef boost::numeric::odeint::controlled_runge_kutta< boost::numeric::odeint::runge_kutta_cash_karp54< state_type_2d_ofp,double,state_type_2d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp>> stepper_type;
/*Debug code
 * int i=0;
        while(i==1){
            sleep(500);
}*/


        integrate_adaptive( stepper_type() , System , x0 , t , tf , dt);
        fu=x0.data.get<0>();
        fv=x0.data.get<1>();
        domain.write("ReactDiffTest");
/*        //size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4<state_type_struct_ofp>(),Fitz,x0,0.0,tf,dt);
        boost::numeric::odeint::runge_kutta4< state_type_struct_ofp<vector_type>,double,state_type_struct_ofp<vector_type>,double,boost::numeric::odeint::vector_space_algebra_ofp> rk4;
        x0.data.get<0>() = u;
        x0.data.get<1>() = v;
        int time_steps = 5;
        for (int step = 1; step <= time_steps; step++)
        {
            std::cout << "step " << step << std::endl;
            rk4.do_step(System,x0,t,dt);
            fu=x0.data.get<0>();
            fv=x0.data.get<1>();
            domain.write_frame("ReactDiff",step);
        }
        Lap.deallocate(domain);*/
}


BOOST_AUTO_TEST_SUITE_END()
