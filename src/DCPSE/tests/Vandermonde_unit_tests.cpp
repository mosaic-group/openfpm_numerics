//
// Created by tommaso on 22/03/19.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <Vector/vector_dist.hpp>
#include <DCPSE/MonomialBasis.hpp>
#include <DCPSE/VandermondeRowBuilder.hpp>
#include <DCPSE/Vandermonde.hpp>
#include "DMatrix/EMatrix.hpp"

BOOST_AUTO_TEST_SUITE(Vandermonde_tests)

// If EIGEN is not present, EMatrix is not available and we don't need to build this test
#ifdef HAVE_EIGEN

    BOOST_AUTO_TEST_CASE(VandermondeRowBuilder_AllOnes_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> row(1, mb.size());
        Point<2, double> x({2, 2});
        double eps = 2;
        VandermondeRowBuilder<2, double> vrb(mb);
        vrb.buildRow(row, 0, x / eps);
        // For the way the row has been constructed, it should be composed of only 1s
        bool isRowAllOnes = true;
        for (int i = 0; i < mb.size(); ++i)
        {
            isRowAllOnes = isRowAllOnes && (row(0, i) == 1);
        }
        BOOST_REQUIRE(isRowAllOnes);
    }

    BOOST_AUTO_TEST_CASE(VandermondeRowBuilder_OneZero_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> row(1, mb.size());
        Point<2, double> x({1, 0});
        double eps = 1;
        VandermondeRowBuilder<2, double> vrb(mb);
        vrb.buildRow(row, 0, x / eps);
        // For the way the row has been constructed, it should be composed of only 1s
        bool areValuesOk = true;
        for (int i = 0; i < mb.size(); ++i)
        {
            bool isThereY = mb.getElement(i).getExponent(1) > 0;
            bool curCheck = (row(0, i) == !isThereY);
            areValuesOk = areValuesOk && curCheck;
        }
        BOOST_REQUIRE(areValuesOk);
    }

    BOOST_AUTO_TEST_CASE(VandermondeRowBuilder_ZeroOne_test)
    {
        MonomialBasis<2> mb({1, 0}, 2);
        EMatrix<double, Eigen::Dynamic, Eigen::Dynamic> row(1, mb.size());
        Point<2, double> x({0, 1});
        double eps = 1;
        VandermondeRowBuilder<2, double> vrb(mb);
        vrb.buildRow(row, 0, x / eps);
        // For the way the row has been constructed, it should be composed of only 1s
        bool areValuesOk = true;
        for (int i = 0; i < mb.size(); ++i)
        {
            bool isThereX = mb.getElement(i).getExponent(0) > 0;
            bool curCheck = (row(0, i) == !isThereX);
            areValuesOk = areValuesOk && curCheck;
        }
        BOOST_REQUIRE(areValuesOk);
    }

#endif // HAVE_EIGEN

BOOST_AUTO_TEST_SUITE_END()
