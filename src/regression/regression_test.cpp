/*
 * Unit tests for the regression module
 * author : Sachin (sthekke@mpi-cbg.de)
 * date : 28.04.2022
 *
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "regression.hpp"
#include "Vector/vector_dist.hpp"
#include "DMatrix/EMatrix.hpp"

BOOST_AUTO_TEST_SUITE( Regression_test )



BOOST_AUTO_TEST_CASE ( Regression_local )
{
    Box<2,float> domain({0.0,0.0},{1.0,1.0});
    // Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<2,float> g(0.01);

    using vectorType = vector_dist<2,float, aggregate<double> >;
    vectorType vd(2048,domain,bc,g);
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

    auto it2 = vd.getDomainIterator();
    auto NN = vd.getCellList(0.1);
    /*
    auto support = RegressionSupport<vectorType, decltype(NN)>(vd, it2, 10, N_PARTICLES, NN);
    
    std::cout<<"Initialized domain with "<<support.getNumParticles()<<" particles.\n";
    
    auto model = RegressionModel<2, 0, vectorType>(vd, dom, 6, 2.0);

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

    double tolerance = 1e-5;
    bool check;
    if (std::abs(max_err) < tolerance)
        check = true;
    else
        check = false;
    std::cout<<"Max err = "<<max_err<<"\n";
    BOOST_TEST( check );
    */
}



BOOST_AUTO_TEST_CASE ( Regression_without_domain_initialization)
{
    Box<2,float> domain({0.0,0.0},{1.0,1.0});
    // Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};
    // extended boundary around the domain, and the processor domain
    Ghost<2,float> g(0.01);

    using vectorType = vector_dist<2,float, aggregate<double> >;
    
    vectorType vd_orig(1024,domain,bc,g);

    vectorType vd_test(vd_orig.getDecomposition(), 200);
    
    Vcluster<> & v_cl = create_vcluster();
    long int N_prc = v_cl.getProcessingUnits();

    if (v_cl.getProcessUnitID() == 0)
        std::cout<<"Running regression test with "<<N_prc<<" procs.\n";

    // the scalar is the element at position 0 in the aggregate
    const int scalar = 0;

    // Creating a grid with synthetic data
    {
        auto it = vd_orig.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            double posx = (double)rand() / RAND_MAX;
            double posy = (double)rand() / RAND_MAX;

            vd_orig.getPos(key)[0] = posx;
            vd_orig.getPos(key)[1] = posy;

            // Use synthetic data x*y
            vd_orig.template getProp<scalar>(key) = sin(posx*posy);
            ++it;
        }    
        vd_orig.map();
    }

    // Creating test points
    {
        auto it = vd_test.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            double posx = (double)rand() / RAND_MAX;
            double posy = (double)rand() / RAND_MAX;

            vd_test.getPos(key)[0] = posx;
            vd_test.getPos(key)[1] = posy;

            vd_test.template getProp<scalar>(key) = 0.0;
            
            ++it;
        }    
        vd_test.map();
    }
    
    auto model = RegressionModel<2, 0>(vd_orig, 1e-6);

    double max_err = -1.0;
    // Checking the error
    {
        auto it = vd_test.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            Point<2, double>  pos = {vd_test.getPos(key)[0], vd_test.getPos(key)[1]};
            double val = model.eval(pos);
            double actual = sin(pos[0]*pos[1]);
            double err = std::abs(actual - val);
            if (err > max_err) max_err = err;
           
            vd_test.template getProp<scalar>(key) = val;
        
            ++it;
        }    
        vd_test.ghost_get<scalar>();
    }

    v_cl.max(max_err);
	v_cl.execute();
	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Maximum error: " << max_err << "\n";

    double tolerance = 1e-5;
    bool check;
    if (std::abs(max_err) < tolerance)
        check = true;
    else
        check = false;

    BOOST_TEST( check );


}


BOOST_AUTO_TEST_SUITE_END()
