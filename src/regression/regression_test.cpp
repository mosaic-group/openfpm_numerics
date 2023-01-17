/*
 * Unit tests for the regression module
 * author : Sachin (sthekke@mpi-cbg.de)
 * date : 28.04.2022
 *
 */

#ifndef OPENFPM_NUMERICS_REGRESSION_TEST_HPP_
#define OPENFPM_NUMERICS_REGRESSION_TEST_HPP_


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "regression.hpp"
#include "Vector/vector_dist.hpp"
#include "DMatrix/EMatrix.hpp"

BOOST_AUTO_TEST_SUITE( Regression_test )



BOOST_AUTO_TEST_CASE ( Regression_domain_initialization )
{
    Box<2,float> domain({0.0,0.0},{1.0,1.0});
    // Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<2,float> g(0.01);

    using vectorType = vector_dist<2,float, aggregate<double> >;
    vectorType vd(4096,domain,bc,g);
    // the scalar is the element at position 0 in the aggregate
    const int scalar = 0;

    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto key = it.get();
        double posx = (double)rand() / RAND_MAX;
        double posy = (double)rand() / RAND_MAX;

        vd.getPos(key)[0] = posx;
        vd.getPos(key)[1] = posy;

        // Use synthetic data x*y
        vd.template getProp<scalar>(key) = sin(posx*posy);
        ++it;
    }    
    vd.map();

    double size = 0.075;
    auto NN = vd.getCellList(size);
    auto dom = new RegressionDomain<vectorType, decltype(NN)>(vd, {0.8,0.8}, size, NN);
    
    std::cout<<"Initialized domain with "<<dom->getNumParticles()<<" particles.\n";

    auto model = new RegressionModel<2, 0, vectorType, EMatrixXd, EVectorXd>(vd, dom, 6, 2.0);

    double max_err = -1.0;
    for(double x = 0.75; x < 0.85;x+=0.01)
    {
        for(double y = 0.75; y < 0.85; y+=0.01)
        {
            Point<2, double> pos{x,y};
            double val = model->eval(pos);
            double actual = sin(x*y);
            double err = std::abs(actual - val);
            if (err > max_err) max_err = err;
        }
    }

    double tolerance = 1e-7;
    bool check;
    if (std::abs(max_err) < tolerance)
        check = true;
    else
        check = false;
    std::cout<<"Max err = "<<max_err<<"\n";
    BOOST_TEST( check );

    if(model)
        delete model;
    
    if(dom)
        delete dom;

}



BOOST_AUTO_TEST_CASE ( Regression_without_domain_initialization)
{
    Box<2,float> domain({0.0,0.0},{1.0,1.0});
    // Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<2,float> g(0.01);

    using vectorType = vector_dist<2,float, aggregate<double> >;
    vectorType vd(4096,domain,bc,g);
    // the scalar is the element at position 0 in the aggregate
    const int scalar = 0;

    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto key = it.get();
        double posx = (double)rand() / RAND_MAX;
        double posy = (double)rand() / RAND_MAX;

        vd.getPos(key)[0] = posx;
        vd.getPos(key)[1] = posy;

        // Use synthetic data x*y
        vd.template getProp<scalar>(key) = sin(posx*posy);
        ++it;
    }    
    vd.map();
    
    auto model = new RegressionModel<2, 0, vectorType, EMatrixXd, EVectorXd>(vd, static_cast<unsigned int>(10));

    double max_err = -1.0;
    for(double x = 0.75; x < 0.85;x+=0.01)
    {
        for(double y = 0.75; y < 0.85; y+=0.01)
        {
            Point<2, double> pos{x,y};
            double val = model->eval(pos);
            double actual = sin(x*y);
            double err = std::abs(actual - val);
            if (err > max_err) max_err = err;
        }
    }

    double tolerance = 1e-7;
    bool check;
    if (std::abs(max_err) < tolerance)
        check = true;
    else
        check = false;
    std::cout<<"Max err = "<<max_err<<"\n";
    BOOST_TEST( check );

    if(model)
        delete model;

}


BOOST_AUTO_TEST_SUITE_END()

#endif /* OPENFPM_NUMERICS_REGRESSION_TEST_HPP_ */
