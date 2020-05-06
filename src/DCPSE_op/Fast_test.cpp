//
// Created by Abhinav Singh on 25.04.20.
//
#include "config.h"

#define BOOST_TEST_DYN_LINK

#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "DCPSE_op.hpp"
#include "DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "EqnsStruct.hpp"
BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests)
    BOOST_AUTO_TEST_CASE(dcpse_Kernels) {
        const size_t sz[2] = {31,31};
        Box<2, double> box({0, 0}, {10,10});
        size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut =spacing*(3.1);
        double ord = 2;
        double sampling = 1.9;
        double rCut2 = 3.1*spacing;
        double ord2 = 2;
        double sampling2 = 1.9;

        double sigma2 = spacing * spacing / (2 * 4);
        std::mt19937 rng{7};
        std::normal_distribution<> gaussian{0, sigma2};

        Ghost<2, double> ghost(rCut);
        std::cout<<spacing<<std::endl;
        //                                  sf    W   DW        RHS    Wnew  V
        vector_dist<2, double, aggregate<double,double,double,double,double,VectorS<2, double>>> Particles(0, box, bc, ghost);
        auto f1 = getV<0>(Particles);
        auto f2 = getV<1>(Particles);
        auto f3 = getV<2>(Particles);
        auto f4 = getV<3>(Particles);
        auto f5 = getV<4>(Particles);
        auto f6 = getV<5>(Particles);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            Particles.getLastProp<0>() =0;
            Particles.getLastProp<1>() =0;
            Particles.getLastProp<2>() =0.0;
            Particles.getLastProp<3>() =0.0;
            Particles.getLastProp<4>() =0.0;
            Particles.getLastProp<5>()[0] =x;
            Particles.getLastProp<5>()[1] =y;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> up_p;
        openfpm::vector<aggregate<int>> dw_p;
        openfpm::vector<aggregate<int>> l_p;
        openfpm::vector<aggregate<int>> r_p;

        openfpm::vector<aggregate<int>> up_p1;
        openfpm::vector<aggregate<int>> dw_p1;
        openfpm::vector<aggregate<int>> l_p1;
        openfpm::vector<aggregate<int>> r_p1;

        openfpm::vector<aggregate<int>> BULK;



        // Here fill up the boxes for particle detection.

        Box<2, double> up({box.getLow(0) + spacing / 2.0, box.getHigh(1) - spacing / 2.0},
                          {box.getHigh(0) -  spacing / 2.0, box.getHigh(1) + spacing / 2.0});

        Box<2, double> up_d({box.getLow(0) +  spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0},
                            {box.getHigh(0) -  spacing / 2.0, box.getHigh(1) - spacing / 2.0});

        Box<2, double> down({box.getLow(0) +  spacing / 2.0, box.getLow(1) - spacing / 2.0},
                            {box.getHigh(0) -  spacing / 2.0, box.getLow(1) + spacing / 2.0});

        Box<2, double> down_u({box.getLow(0) +   spacing / 2.0, box.getLow(1) + spacing / 2.0},
                              {box.getHigh(0) -  spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0});

        Box<2, double> left({box.getLow(0) - spacing / 2.0, box.getLow(1) -  spacing / 2.0},
                            {box.getLow(0) + spacing / 2.0, box.getHigh(1) +  spacing / 2.0});

        Box<2, double> left_r({box.getLow(0) + spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                              {box.getLow(0) + 3 * spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});

        Box<2, double> right({box.getHigh(0) - spacing / 2.0, box.getLow(1) -  spacing / 2.0},
                             {box.getHigh(0) + spacing / 2.0, box.getHigh(1) +  spacing / 2.0});

        Box<2, double> right_l({box.getHigh(0) - 3 * spacing / 2.0, box.getLow(1) + 3 * spacing / 2.0},
                               {box.getHigh(0) - spacing / 2.0, box.getHigh(1) - 3 * spacing / 2.0});
        Box<2, double> Bulk_box({box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0},
                                {box.getHigh(0) + spacing / 2.0, box.getHigh(1)  +spacing / 2.0});


        openfpm::vector<Box<2, double>> boxes;
        boxes.add(up);
        boxes.add(up_d);
        boxes.add(down);
        boxes.add(down_u);
        boxes.add(left);
        boxes.add(left_r);
        boxes.add(right);
        boxes.add(right_l);
        boxes.add(Bulk_box);
        VTKWriter<openfpm::vector<Box<2, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("box_sf.vtk");
        // Particles.write_frame("Re1000-1e-3-Lid_sf",0);

        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<2, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p) =0;//-M_PI*M_PI*sin(M_PI*xp[0])*sin(M_PI*xp[1]);
            Particles.getProp<1>(p) =0;//sin(M_PI*xp[0])*sin(M_PI*xp[1]);
            Particles.getProp<2>(p) =0.0;
            Particles.getProp<3>(p) =0.0;
            Particles.getProp<4>(p) =0.0;
            BULK.add();
            BULK.last().get<0>() = p.getKey();

            if (up.isInside(xp) == true) {
                up_p.add();
                up_p.last().get<0>() = p.getKey();
            }
            else if (down.isInside(xp) == true) {
                dw_p.add();
                dw_p.last().get<0>() = p.getKey();
            }
            else if (left.isInside(xp) == true) {
                l_p.add();
                l_p.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                r_p.add();
                r_p.last().get<0>() = p.getKey();
            }
            else if (Bulk_box.isInside(xp) == true)
            {
                if (up_d.isInside(xp) == true) {
                    up_p1.add();
                    up_p1.last().get<0>() = p.getKey();
                }
                else if (down_u.isInside(xp) == true) {
                    dw_p1.add();
                    dw_p1.last().get<0>() = p.getKey();
                }else if (left_r.isInside(xp) == true) {
                    l_p1.add();
                    l_p1.last().get<0>() = p.getKey();
                } else if (right_l.isInside(xp) == true) {
                    r_p1.add();
                    r_p1.last().get<0>() = p.getKey();
                }
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }

        for(int j=0;j<bulk.size();j++) {
            auto p = bulk.get<0>(j);
            Particles.getPos(p)[0]=Particles.getPos(p)[0]+gaussian(rng);
            Particles.getPos(p)[1]=Particles.getPos(p)[1]+gaussian(rng);
        }



        Laplacian Lap(Particles, ord2, rCut2,sampling,support_options::RADIUS);
        Derivative_xy Dxy(Particles,ord, rCut, sampling,support_options::RADIUS);
        Derivative_xy Dxy2(Particles,ord2, rCut2, sampling2,support_options::RADIUS);
        Derivative_x Dx(Particles,ord, rCut, sampling,support_options::RADIUS);
        Derivative_y Dy(Particles,ord, rCut, sampling,support_options::RADIUS);

        std::cout<<"Dx"<<std::endl;
        Dx.checkMomenta(Particles);
        std::cout<<"Dy"<<std::endl;
        Dy.checkMomenta(Particles);
        std::cout<<"Dxy"<<std::endl;
        Dxy.checkMomenta(Particles);
        std::cout<<"Lap"<<std::endl;
        Lap.checkMomenta(Particles);
        auto its2 = Particles.getDomainIterator();
        int ctr=0;
        while (its2.isNext()) {
            auto p = its2.get();
            Dx.DrawKernel<0>(Particles, p.getKey());
            Dy.DrawKernel<1>(Particles, p.getKey());
            Dxy.DrawKernel<2>(Particles, p.getKey());
            Dxy2.DrawKernel<3>(Particles, p.getKey());
            Lap.DrawKernel<4>(Particles, p.getKey());
            Particles.write_frame("Kernel_moved",ctr);
            f1=0;
            f2=0;
            f3=0;
            f4=0;
            f5=0;
            ++its2;
            ctr++;
        }


        for(int j=0;j<bulk.size();j++) {
            auto p = bulk.get<0>(j);
            Particles.getPos(p)[0]=Particles.getProp<5>(p)[0];
            Particles.getPos(p)[1]=Particles.getProp<5>(p)[1];
        }

        Particles.map();
        Dx.update(Particles);
        Dy.update(Particles);
        Dxy.update(Particles);
        auto Dyx = Dxy;
        Lap.update(Particles);

        std::cout<<"Dx"<<std::endl;
        Dx.checkMomenta(Particles);
        std::cout<<"Dy"<<std::endl;
        Dy.checkMomenta(Particles);
        std::cout<<"Dxy"<<std::endl;
        Dxy.checkMomenta(Particles);
        std::cout<<"Lap"<<std::endl;
        Lap.checkMomenta(Particles);


        auto its = Particles.getDomainIterator();
        ctr=0;
        while (its.isNext()) {
            auto p = its.get();
            Dx.DrawKernel<0>(Particles, p.getKey());
            Dy.DrawKernel<1>(Particles, p.getKey());
            Dxy.DrawKernel<2>(Particles, p.getKey());
            Dxy2.DrawKernel<3>(Particles, p.getKey());
            Lap.DrawKernel<4>(Particles, p.getKey());
            Particles.write_frame("Kernel_unmoved",ctr);
            f1=0;
            f2=0;
            f3=0;
            f4=0;
            ++its;
            ctr++;
        }

    }

BOOST_AUTO_TEST_SUITE_END()
