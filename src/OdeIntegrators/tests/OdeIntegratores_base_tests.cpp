//
// Created by Abhinav Singh on 01.12.20.
//
#include "config.h"


#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"

typedef texp_v<double> state_type;

void Exponential( const state_type &x , state_type &dxdt , const double t )
{
    dxdt = x;
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


BOOST_AUTO_TEST_SUITE_END()
