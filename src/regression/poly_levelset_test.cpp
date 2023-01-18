/*
 * Unit tests for the Regression module PolyLevelset submodule
 * author : Sachin (sthekke@mpi-cbg.de)
 * date : 18.01.2023
 *
 */


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "poly_levelset.hpp"
#include "Vector/vector_dist.hpp"
#include "DMatrix/EMatrix.hpp"

BOOST_AUTO_TEST_SUITE( Regression_test )



BOOST_AUTO_TEST_CASE ( PolyLevelset_Sphere )
{
    Box<3,float> domain({-2.0,-2.0,-2.0},{2.0,2.0,2.0});

    // Here we define the boundary conditions of our problem
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

    // extended boundary around the domain, and the processor domain
    Ghost<3,float> g(0.01);

    using vectorType = vector_dist<3,float, aggregate<double, double> >;

    vectorType vd(1024,domain,bc,g);

    constexpr int mean_curvature = 0;
    constexpr int gauss_curvature = 1;

    // Initialize points on sphere
    auto it = vd.getDomainIterator();
    while (it.isNext())
    {
        auto key = it.get();
        double theta = ((double)rand() / RAND_MAX) * M_PI;
        double phi = ((double)rand() / RAND_MAX) * 2.0 * M_PI;

        vd.getPos(key)[0] = cos(theta) * sin(phi);
        vd.getPos(key)[1] = cos(theta) * cos(phi);
        vd.getPos(key)[2] = sin(theta);

        vd.template getProp<mean_curvature>(key) = 0.0;
        vd.template getProp<gauss_curvature>(key) = 0.0;
        
        ++it;
    }    
    vd.map();
    
    auto model = new PolyLevelset<3, vectorType>(vd, 1e-4);

    double max_err = -1.0;
    auto it2 = vd.getDomainIterator();
    while (it2.isNext())
    {
        auto key = it2.get();
        vd.template getProp<mean_curvature>(key) = model->estimate_mean_curvature_at(vd.getPos(key));
        // vd.template getProp<gauss_curvature>(key) = model->estimate_gauss_curvature_at(vd.getPos(key));

        double val = vd.getProp<mean_curvature>(key);
        double actual = 1.0;
        double err = std::abs(actual - val);
        if (err > max_err) max_err = err;

        ++it2;
    }


    double tolerance = 1e-4;
    bool check;
    if (std::abs(max_err) < tolerance)
        check = true;
    else
        check = false;
    std::cout<<"Max err (poly level) = "<<max_err<<"\n";
    BOOST_TEST( check );

    if(model)
        delete model;

}


BOOST_AUTO_TEST_SUITE_END()

