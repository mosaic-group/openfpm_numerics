//
// Created by abhinav on 2/28/23.
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
#include "OdeIntegrators/OdeIntegrators.hpp"
//#include "DCPSE/DCPSE_op/DCPSE_op.hpp"

typedef state_type_1d_ofp_gpu state_type;
//const double a = 2.8e-4;
//const double b = 5e-3;
//const double tau = .1;
//const double k = .005;

void ExponentialGPU( const state_type &x , state_type &dxdt , const double t )
{
    dxdt.data.get<0>() = x.data.get<0>();
    //x.data.get<0>().getVector().deviceToHost<0>();
    //dxdt.data.get<0>().getVector().deviceToHost<0>();
}

BOOST_AUTO_TEST_SUITE(odeInt_BASE_tests)

BOOST_AUTO_TEST_CASE(odeint_base_test_gpu)
        {
        size_t edgeSemiSize = 512;
        const size_t sz[2] = {edgeSemiSize,edgeSemiSize };
        Box<2, double> box({ 0, 0 }, { 1.0, 1.0 });
        size_t bc[2] = { NON_PERIODIC, NON_PERIODIC };
        double spacing[2];
        spacing[0] = 1.0 / (sz[0] - 1);
        spacing[1] = 1.0 / (sz[1] - 1);
        double rCut = 3.9 * spacing[0];
        Ghost<2, double> ghost(rCut);
        BOOST_TEST_MESSAGE("Init vector_dist...");

        vector_dist_gpu<2, double, aggregate<double, double,double>> Particles(0, box, bc, ghost);

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
            Particles.getLastProp<0>() = xp0*yp0*exp(-5);
            Particles.getLastProp<1>() = xp0*yp0*exp(5);
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();
        Particles.hostToDeviceProp<0,1,2>();
        auto Init = getV<0,comp_dev>(Particles);
        auto Sol = getV<1,comp_dev>(Particles);
        auto OdeSol = getV<2,comp_dev>(Particles);

        state_type x0;
        x0.data.get<0>()=Init;
        x0.data.get<0>().getVector().deviceToHost<0>();
        // The rhs of x' = f(x)

        double t0=-5,tf=5;
        const double dt=0.01;


        //This doesnt work Why?
        //size_t steps=boost::numeric::odeint::integrate(Exponential,x0,0.0,tf,dt);
        timer tt;
        tt.start();
        size_t steps=boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra_ofp_gpu,boost::numeric::odeint::ofp_operations>(),ExponentialGPU,x0,t0,tf,dt);
        tt.stop();
        OdeSol=x0.data.get<0>();
        Particles.deviceToHostProp<0,1,2>();
        auto it2 = Particles.getDomainIterator();
        double worst = 0.0;
        while (it2.isNext()) {
            auto p = it2.get();
            if (fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p)) > worst) {
                worst = fabs(Particles.getProp<1>(p) - Particles.getProp<2>(p));
            }
            ++it2;
        }
        std::cout<<"WCT:"<<tt.getwct()<<std::endl;
        std::cout<<"CPU:"<<tt.getcputime()<<std::endl;
        std::cout<<worst<<std::endl;
        BOOST_REQUIRE(worst < 1e-6);
        }
BOOST_AUTO_TEST_SUITE_END()
