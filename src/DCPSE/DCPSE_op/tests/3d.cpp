//
// Created by Abhinav Singh on 15.06.20.
//

#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40

#include "config.h"

#define BOOST_TEST_DYN_LINK


#include "util/util_debug.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "../DCPSE_op.hpp"
#include "../DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "../EqnsStruct.hpp"
#include "Decomposition/Distribution/SpaceDistribution.hpp"
template<typename particle_type,typename particle_type2>
void indexUpdate(
        particle_type &Particles,particle_type2 &Particles_subset,openfpm::vector<aggregate<int>> &Boundary,openfpm::vector<aggregate<int>> &bulk,
        Box<3, double> &up, Box<3, double> &down, Box<3, double> &left, Box<3, double> &right,Box<3, double> &front, Box<3, double> &back) {
    Boundary.clear();
    bulk.clear();
    Particles_subset.clear();

    auto it2 = Particles.getDomainIterator();
    while (it2.isNext()) {
        auto p = it2.get();
        Point<3, double> xp = Particles.getPos(p);
        if (front.isInside(xp) == true) {
            Boundary.add();
            Boundary.last().get<0>() = p.getKey();
        } else if (back.isInside(xp) == true) {
            Boundary.add();
            Boundary.last().get<0>() = p.getKey();
        } else if (left.isInside(xp) == true) {
            Boundary.add();
            Boundary.last().get<0>() = p.getKey();
        } else if (right.isInside(xp) == true) {
            Boundary.add();
            Boundary.last().get<0>() = p.getKey();
        } else if (up.isInside(xp) == true) {
            Boundary.add();
            Boundary.last().get<0>() = p.getKey();
        } else if (down.isInside(xp) == true) {
            Boundary.add();
            Boundary.last().get<0>() = p.getKey();
        } else {
            bulk.add();
            bulk.last().get<0>() = p.getKey();
        }
        ++it2;
    }
    for (int i = 0; i < bulk.size(); i++) {
        Particles_subset.add();
        Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
        Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
    }

}

BOOST_AUTO_TEST_SUITE(dcpse_op_suite_tests3)

    BOOST_AUTO_TEST_CASE(Active3dSimple) {
        timer tt2;
        tt2.start();
        size_t grd_sz = 21;
        double dt = 1e-3;
        double V_err_eps = 5e-2;
        double boxsize = 10;
        const size_t sz[3] = {grd_sz, grd_sz, grd_sz};
        Box<3, double> box({0, 0, 0}, {boxsize, boxsize, boxsize});
        size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
        double Lx = box.getHigh(0);
        double Ly = box.getHigh(1);
        double Lz = box.getHigh(2);
        double spacing = box.getHigh(0) / (sz[0] - 1);
        double rCut = 3.9 * spacing;
        double rCut2 = 3.9 * spacing;
        int ord = 2;
        int ord2 = 2;
        double sampling_factor = 4.0;
        double sampling_factor2 = 2.4;
        Ghost<3, double> ghost(rCut);
        auto &v_cl = create_vcluster();
        /*                                 pol                          V           vort              Ext                 Press     strain       stress                      Mfield,   dPol                      dV         RHS                  f1     f2     f3    f4     f5     f6       H               V_t      div    H_t                                                                                      delmu */
        vector_dist<3, double, aggregate<VectorS<3, double>, VectorS<3, double>, double[3][3], VectorS<3, double>, double, double[3][3], double[3][3], VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, double, double, double, double, double, double, double, VectorS<3, double>, double, double, double[3], double[3], double[3], double[3], double[3], double[3], double, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, double, double, double, double[3][3],double[3][3]>> Particles(
                0, box, bc, ghost);


        double x0, y0, z0, x1, y1, z1;
        x0 = box.getLow(0);
        y0 = box.getLow(1);
        z0 = box.getLow(2);
        x1 = box.getHigh(0);
        y1 = box.getHigh(1);
        z1 = box.getHigh(2);

        auto it = Particles.getGridIterator(sz);
        while (it.isNext()) {
            Particles.add();
            auto key = it.get();
            double x = key.get(0) * it.getSpacing(0);
            Particles.getLastPos()[0] = x;
            double y = key.get(1) * it.getSpacing(1);
            Particles.getLastPos()[1] = y;
            double z = key.get(2) * it.getSpacing(1);
            Particles.getLastPos()[2] = z;
            ++it;
        }

        Particles.map();
        Particles.ghost_get<0>();

        openfpm::vector<aggregate<int>> bulk;
        openfpm::vector<aggregate<int>> Boundary;

        constexpr int x = 0;
        constexpr int y = 1;
        constexpr int z = 2;


        constexpr int Polarization = 0;
        constexpr int Velocity = 1;
        constexpr int Vorticity = 2;
        constexpr int ExtForce = 3;
        constexpr int Pressure = 4;
        constexpr int Strain_rate = 5;
        constexpr int Stress = 6;
        constexpr int MolField = 7;
        auto Pos = getV<PROP_POS>(Particles);

        auto Pol = getV<Polarization>(Particles);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto dPol = getV<8>(Particles);
        auto dV = getV<9>(Particles);
        auto RHS = getV<10>(Particles);
        auto f1 = getV<11>(Particles);
        auto f2 = getV<12>(Particles);
        auto f3 = getV<13>(Particles);
        auto f4 = getV<14>(Particles);
        auto f5 = getV<15>(Particles);
        auto f6 = getV<16>(Particles);
        auto H = getV<17>(Particles);
        auto V_t = getV<18>(Particles);
        auto div = getV<19>(Particles);
        auto H_t = getV<20>(Particles);
        auto Df1 = getV<21>(Particles);
        auto Df2 = getV<22>(Particles);
        auto Df3 = getV<23>(Particles);
        auto Df4 = getV<24>(Particles);
        auto Df5 = getV<25>(Particles);
        auto Df6 = getV<26>(Particles);
        auto delmu = getV<27>(Particles);
        auto k1 = getV<28>(Particles);
        auto k2 = getV<29>(Particles);
        auto k3 = getV<30>(Particles);
        auto k4 = getV<31>(Particles);
        auto H_p_b = getV<32>(Particles);
        auto FranckEnergyDensity = getV<33>(Particles);
        auto r = getV<34>(Particles);
        auto q = getV<35>(Particles);


        double eta = 1.0;
        double nu = -0.5;
        double gama = 0.1;
        double zeta = 0.07;
        double Ks = 1.0;
        double Kt = 1.1;
        double Kb = 1.5;
        double lambda = 0.1;
        //double delmu = -1.0;
        g = 0;
        delmu = -1.0;
        P = 0;
        V = 0;
        // Here fill up the boxes for particle Boundary detection.
        Particles.ghost_get<ExtForce, 27>(SKIP_LABELLING);

        // Here fill up the boxes for particle detection.

        Box<3, double> up(
                {box.getLow(0) - spacing / 2.0, box.getHigh(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> down(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getLow(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> left(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getLow(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> right(
                {box.getHigh(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        Box<3, double> front(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getLow(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getLow(2) + spacing / 2.0});

        Box<3, double> back(
                {box.getLow(0) - spacing / 2.0, box.getLow(1) - spacing / 2.0, box.getHigh(2) - spacing / 2.0},
                {box.getHigh(0) + spacing / 2.0, box.getHigh(1) + spacing / 2.0, box.getHigh(2) + spacing / 2.0});

        openfpm::vector<Box<3, double>> boxes;
        boxes.add(up);
        boxes.add(down);
        boxes.add(left);
        boxes.add(right);
        boxes.add(front);
        boxes.add(back);
        VTKWriter<openfpm::vector<Box<3, double>>, VECTOR_BOX> vtk_box;
        vtk_box.add(boxes);
        vtk_box.write("boxes_3d.vtk");
        auto it2 = Particles.getDomainIterator();

        while (it2.isNext()) {
            auto p = it2.get();
            Point<3, double> xp = Particles.getPos(p);
            Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * xp[x] - Lx) / Lx) - sin((2 * xp[y] - Ly) / Ly)));
            Particles.getProp<0>(p)[z] = 0;
            if (front.isInside(xp) == true) {
                Boundary.add();
                Boundary.last().get<0>() = p.getKey();
            } else if (back.isInside(xp) == true) {
                Boundary.add();
                Boundary.last().get<0>() = p.getKey();
            } else if (left.isInside(xp) == true) {
                Boundary.add();
                Boundary.last().get<0>() = p.getKey();
            } else if (right.isInside(xp) == true) {
                Boundary.add();
                Boundary.last().get<0>() = p.getKey();
            } else if (up.isInside(xp) == true) {
                Boundary.add();
                Boundary.last().get<0>() = p.getKey();
            } else if (down.isInside(xp) == true) {
                Boundary.add();
                Boundary.last().get<0>() = p.getKey();
            } else {
                bulk.add();
                bulk.last().get<0>() = p.getKey();
            }
            ++it2;
        }

        //vector_dist_subset<3, double, aggregate<VectorS<3, double>, VectorS<3, double>, double[3][3], VectorS<3, double>, double, double[3][3], double[3][3], VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, double, double, double, double, double, double, double, VectorS<3, double>, double, double, double[3], double[3], double[3], double[3], double[3], double[3], double, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, VectorS<3, double>, double, double, double, double[3][3],double[3][3]>> Particles_subset(Particles, bulk);
        vector_dist<3, double, aggregate<double, VectorS<3, double>>> Particles_subset(Particles.getDecomposition(), 0);
        for (int i = 0; i < bulk.size(); i++) {
            Particles_subset.add();
            Particles_subset.getLastPos()[0] = Particles.getPos(bulk.template get<0>(i))[0];
            Particles_subset.getLastPos()[1] = Particles.getPos(bulk.template get<0>(i))[1];
        }

        Particles_subset.map();
        Particles_subset.ghost_get<0>();

        auto Grad_bulk = getV<1>(Particles_subset);
        auto P_bulk = getV<0>(Particles_subset);


        //Particles.write_frame("Active3d_Parts", 0);

        //Particles_subset.write("Pars");
        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dx(Particles_subset, ord,
                                                                                                 rCut, sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dy(Particles_subset, ord,
                                                                                                 rCut, sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_z Dz(Particles, ord, rCut, sampling_factor, support_options::RADIUS), Bulk_Dz(Particles_subset, ord,
                                                                                                 rCut, sampling_factor,
                                                                                                 support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yz Dyz(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_xz Dxz(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        auto Dyx = Dxy;
        auto Dzy = Dyz;
        auto Dzx = Dxz;

        Derivative_xx Dxx(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);//, Dxx2(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);
        Derivative_zz Dzz(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);

        V_t = V;
        auto px=Pol[x];
        auto py=Pol[y];
        auto pz=Pol[z];

        auto dxpx=Dx(Pol[x]);
        auto dxpy=Dx(Pol[y]);
        auto dxpz=Dx(Pol[z]);
        auto dypx=Dy(Pol[x]);
        auto dypy=Dy(Pol[y]);
        auto dypz=Dy(Pol[z]);
        auto dzpx=Dz(Pol[x]);
        auto dzpy=Dz(Pol[y]);
        auto dzpz=Dz(Pol[z]);
        auto dxxpx=Dxx(Pol[x]);
        auto dxxpy=Dxx(Pol[y]);
        auto dxxpz=Dxx(Pol[z]);
        auto dyypx=Dyy(Pol[x]);
        auto dyypy=Dyy(Pol[y]);
        auto dyypz=Dyy(Pol[z]);
        auto dzzpx=Dzz(Pol[x]);
        auto dzzpy=Dzz(Pol[y]);
        auto dzzpz=Dzz(Pol[z]);
        auto dxypx=Dxy(Pol[x]);
        auto dxypy=Dxy(Pol[y]);
        auto dxypz=Dxy(Pol[z]);
        auto dxzpx=Dxz(Pol[x]);
        auto dxzpy=Dxz(Pol[y]);
        auto dxzpz=Dxz(Pol[z]);
        auto dyzpx=Dyz(Pol[x]);
        auto dyzpy=Dyz(Pol[y]);
        auto dyzpz=Dyz(Pol[z]);
        auto dxhx=Dx(h[x]);
        auto dxhy=Dx(h[y]);
        auto dxhz=Dx(h[z]);
        auto dyhx=Dy(h[x]);
        auto dyhy=Dy(h[y]);
        auto dyhz=Dy(h[z]);
        auto dzhx=Dz(h[x]);
        auto dzhy=Dz(h[y]);
        auto dzhz=Dz(h[z]);

        auto dxqxx=Dx(Pol[x]*Pol[x]-1/3*(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]));
        auto dyqxy=Dy(Pol[x]*Pol[x]);
        auto dzqxz=Dz(Pol[x]*Pol[z]);
        auto dxqyx=Dx(Pol[y]*Pol[x]);
        auto dyqyy=Dy(Pol[y]*Pol[y]-1/3*(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]));
        auto dzqyz=Dz(Pol[y]*Pol[z]);
        auto dxqzx=Dx(Pol[z]*Pol[x]);
        auto dyqzy=Dy(Pol[z]*Pol[y]);
        auto dzqzz=Dz(Pol[z]*Pol[z]-1/3*(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]));


        eq_id vx, vy, vz;

        vx.setId(0);
        vy.setId(1);
        vz.setId(2);
        timer tt;
        double V_err = 1, V_err_old;
        int n = 0;
        int nmax = 300;
        int ctr = 0, errctr, Vreset = 0;
        double tim = 0;
        double tf = 1.024;
        div = 0;
        double sum, sum1, sum_k;
        while (tim <= tf) {
            tt.start();
            petsc_solver<double> solverPetsc;
            solverPetsc.setSolver(KSPGMRES);
            //solverPetsc.setRestart(250);
            solverPetsc.setPreconditioner(PCJACOBI);
            Particles.ghost_get<Polarization>(SKIP_LABELLING);



            FranckEnergyDensity = 0.5*Ks*(dxpx + dypy + dzpz)*(dxpx + dypy + dzpz) +
                                  0.5*Kt*((dypz - dzpy)*px + (-dxpz + dzpx)*py + (dxpy - dypx)*pz)*((dypz - dzpy)*px + (-dxpz + dzpx)*py + (dxpy - dypx)*pz) +
                                  0.5*Kb*((-dxpz*px + dzpx*px - dypz*py + dzpy*py)*(-dxpz*px + dzpx*px - dypz*py + dzpy*py) +
                                        (dxpy*py - dypx*py + dxpz*pz - dzpx*pz)*(dxpy*py - dypx*py + dxpz*pz - dzpx*pz) +
                                        (-dxpy*px + dypx*px + dypz*pz - dzpy*pz)*(-dxpy*px + dypx*px + dypz*pz - dzpy*pz));

            h[x]=Ks*(dxxpx + dxypy + dxzpz) +
                    Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                    py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                    px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
                    Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                    px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

            h[y]= Ks*(dxypx + dyypy + dyzpz) +
                    Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                    py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                    px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
                    Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                    px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));


            h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
                    Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                    py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                    px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
                    Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                    px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));

            //Particles.ghost_get<33>(SKIP_LABELLING);
            sigma[x][x] =
                    -dxpx*(dxpx + dypy + dzpz)*Ks -
                    dxpy*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dxpz*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dxpz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dxpy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
            sigma[x][y] =
                    -dxpy*(dxpx + dypy + dzpz)*Ks - dxpz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) +
                    dxpx*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dxpz*Kt*px*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) +
                    dxpx*Kb*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) -
                    dxpz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) +
                    dxpx*Kb*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz);

            sigma[x][z] =
                    -dxpz*(dxpx + dypy + dzpz)*Ks +
                    dxpy*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dxpx*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dxpx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dxpy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));
            sigma[y][x] =
                    -dypx*(dxpx + dypy + dzpz)*Ks +
                    dypz*Kt*py*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dypy*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    0.5*dypz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dypy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
            sigma[y][y] =
                    -dypy*(dxpx + dypy + dzpz)*Ks -
                    dypz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) -
                    dypz*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dypx*Kt*pz*(-dypz*px + dzpy*px + dxpz*py - dzpx*py -dxpy*pz + dypx*pz) -
                    dypx*Kb*py*(-dxpy*py + dypx*py + (-dxpz + dzpx)*pz) -
                    dypx*Kb*px*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) -
                    dypz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz);
            sigma[y][z] =
                    -dypz*(dxpx + dypy + dzpz)*Ks +
                    dypy*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dypx*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py -dxpy*pz + dypx*pz) -
                    0.5*dypx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dypy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));
            sigma[z][x] =
                    -dzpx*(dxpx + dypy + dzpz)*Ks -
                    dzpz*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) +
                    dzpy*Kt*pz*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dzpz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dzpy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
            sigma[z][y] =
                    -dzpy*(dxpx + dypy + dzpz)*Ks -
                    dzpz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) -
                    dzpz*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dzpx*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dzpx*Kb*py*(-dxpy*py + dypx*py + (-dxpz + dzpx)*pz) -
                    dzpz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) +
                    dzpx*Kb*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz);

            sigma[z][z] =
                    -dzpz*(dxpx + dypy + dzpz)*Ks -
                    dzpx*Kt*py*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dzpy*Kt*px*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dzpx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dzpy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));


            Particles.ghost_get<Stress>(SKIP_LABELLING);
            tt.start();
            Particles.ghost_get<Polarization>(SKIP_LABELLING);



            FranckEnergyDensity = 0.5*Ks*(dxpx + dypy + dzpz)*(dxpx + dypy + dzpz) +
                                  0.5*Kt*((dypz - dzpy)*px + (-dxpz + dzpx)*py + (dxpy - dypx)*pz)*((dypz - dzpy)*px + (-dxpz + dzpx)*py + (dxpy - dypx)*pz) +
                                  0.5*Kb*((-dxpz*px + dzpx*px - dypz*py + dzpy*py)*(-dxpz*px + dzpx*px - dypz*py + dzpy*py) +
                                          (dxpy*py - dypx*py + dxpz*pz - dzpx*pz)*(dxpy*py - dypx*py + dxpz*pz - dzpx*pz) +
                                          (-dxpy*px + dypx*px + dypz*pz - dzpy*pz)*(-dxpy*px + dypx*px + dypz*pz - dzpy*pz));

            h[x]=Ks*(dxxpx + dxypy + dxzpz) +
                 Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                     py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                     px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
                 Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                     px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

            h[y]= Ks*(dxypx + dyypy + dyzpz) +
                  Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                      py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                      px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
                  Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                      px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));


            h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
                 Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                     py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                     px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
                 Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                     px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));

            //Particles.ghost_get<33>(SKIP_LABELLING);
            sigma[x][x] =
                    -dxpx*(dxpx + dypy + dzpz)*Ks -
                    dxpy*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dxpz*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dxpz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dxpy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
            sigma[x][y] =
                    -dxpy*(dxpx + dypy + dzpz)*Ks - dxpz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) +
                    dxpx*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dxpz*Kt*px*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) +
                    dxpx*Kb*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) -
                    dxpz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) +
                    dxpx*Kb*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz);

            sigma[x][z] =
                    -dxpz*(dxpx + dypy + dzpz)*Ks +
                    dxpy*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dxpx*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dxpx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dxpy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));
            sigma[y][x] =
                    -dypx*(dxpx + dypy + dzpz)*Ks +
                    dypz*Kt*py*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dypy*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    0.5*dypz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dypy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
            sigma[y][y] =
                    -dypy*(dxpx + dypy + dzpz)*Ks -
                    dypz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) -
                    dypz*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dypx*Kt*pz*(-dypz*px + dzpy*px + dxpz*py - dzpx*py -dxpy*pz + dypx*pz) -
                    dypx*Kb*py*(-dxpy*py + dypx*py + (-dxpz + dzpx)*pz) -
                    dypx*Kb*px*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) -
                    dypz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz);
            sigma[y][z] =
                    -dypz*(dxpx + dypy + dzpz)*Ks +
                    dypy*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dypx*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py -dxpy*pz + dypx*pz) -
                    0.5*dypx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dypy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));
            sigma[z][x] =
                    -dzpx*(dxpx + dypy + dzpz)*Ks -
                    dzpz*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) +
                    dzpy*Kt*pz*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dzpz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dzpy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
            sigma[z][y] =
                    -dzpy*(dxpx + dypy + dzpz)*Ks -
                    dzpz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) -
                    dzpz*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                    dzpx*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dzpx*Kb*py*(-dxpy*py + dypx*py + (-dxpz + dzpx)*pz) -
                    dzpz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) +
                    dzpx*Kb*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz);

            sigma[z][z] =
                    -dzpz*(dxpx + dypy + dzpz)*Ks -
                    dzpx*Kt*py*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                    dzpy*Kt*px*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                    0.5*dzpx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                    0.5*dzpy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));


            Particles.ghost_get<Stress>(SKIP_LABELLING);
            Particles.ghost_get<MolField>(SKIP_LABELLING);

            dV[x] = 0.5*(-dypy*h[x] + dypx*h[y] + dyhy*px - dyhx*py) + 0.5*(-dzpz*h[x] + dzpx*h[z] + dzhz*px - dzhx*pz) +
                    zeta*delmu*(dxqxx + dyqxy + dzqxz)  +
                    nu*(1/2*(-dypy*h[x] + dypx*h[y] + dyhy*px - dyhx*py) + 1/3*(-dxpx*h[x] - dxpy*h[y] - dxpz*h[z] - dxhx*px - dxhy*py -dxhz*pz) + 1/2*(-dzpz*h[x] + dzpx*h[z] + dzhz*px - dzhx*pz))+
                    Dx(sigma[x][x])+Dy(sigma[x][y])+Dz(sigma[x][z]);


            dV[y] = 0.5*(dxpy*h[x] - dxpx*h[y] - dxhy*px + dxhx*py) + 0.5*(-dzpz*h[y] + dzpy*h[z] + dzhz*py - dzhy*pz) +
                    zeta*delmu*(dxqyx + dyqyy + dzqyz) +
                    nu*(1/2*(dxpy*h[x] - dxpx*h[y] - dxhy*px + dxhx*py) + 1/3*(-dypx*h[x] - dypy*h[y] - dypz*h[z] - dyhx*px - dyhy*py - dyhz*pz) + 1/2*(-dzpz*h[y] + dzpy*h[z] + dzhz*py - dzhy*pz))+
                    Dx(sigma[y][x])+Dy(sigma[y][y])+Dz(sigma[y][z]);

            dV[z]=0.5*(dxpz*h[x] - dxpx*h[z] - dxhz*px + dxhx*pz) + 0.5*(dypz*h[y] - dypy*h[z] - dyhz*py + dyhy*pz) +
                  zeta*delmu*(dxqzx + dyqzy +dzqzz)  +
                  nu*(1/2*(dxpz*h[x] - dxpx*h[z] - dxhz*px + dxhx*pz)+ 1/2*(dypz*h[y] - dypy*h[z] - dyhz*py + dyhy*pz) + 1/3*(-dzpx*h[x] - dzpy*h[y] - dzpz*h[z] - dzhx*px - dzhy*py -dzhz*pz))+
                  Dx(sigma[z][x])+Dy(sigma[z][y])+Dz(sigma[z][z]);

            Particles.ghost_get<9>(SKIP_LABELLING);


            //Particles.write("PolarI");
            //Velocity Solution n iterations



            auto Stokes1 = eta * (2*Dxx(V[x]) + Dxy(V[y]) + Dyy(V[x]) + Dxz(V[z]) + Dzz(V[x]));
            auto Stokes2 = eta * (2*Dyy(V[y]) + Dxx(V[y]) + Dxy(V[x]) + Dyz(V[z]) + Dzz(V[y]));
            auto Stokes3 = eta * (2*Dzz(V[z]) + Dxx(V[z]) + Dxz(V[x]) + Dyy(V[z]) + Dyz(V[y]));
            std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
            tt.start();
            V_err = 1;
            n = 0;
            errctr = 0;
            if (Vreset == 1) {
                P_bulk = 0;
                P = 0;
                Vreset = 0;
            }
            P = 0;
            P_bulk = 0;

            while (V_err >= V_err_eps && n <= nmax) {
                RHS[x] = dV[x];
                RHS[y] = dV[y];
                RHS[z] = dV[z];
                Particles_subset.ghost_get<0>(SKIP_LABELLING);
                Grad_bulk[x] = Bulk_Dx(P_bulk);
                Grad_bulk[y] = Bulk_Dy(P_bulk);
                Grad_bulk[z] = Bulk_Dz(P_bulk);
                for (int i = 0; i < bulk.size(); i++) {
                    Particles.template getProp<10>(bulk.template get<0>(i))[x] += Particles_subset.getProp<1>(i)[x];
                    Particles.template getProp<10>(bulk.template get<0>(i))[y] += Particles_subset.getProp<1>(i)[y];
                    Particles.template getProp<10>(bulk.template get<0>(i))[z] += Particles_subset.getProp<1>(i)[z];
                }
                Particles.ghost_get<10>(SKIP_LABELLING);
                DCPSE_scheme<equations3d3, decltype(Particles)> Solver(Particles);
                Solver.impose(Stokes1, bulk, RHS[0], vx);
                Solver.impose(Stokes2, bulk, RHS[1], vy);
                Solver.impose(Stokes3, bulk, RHS[2], vz);
                Solver.impose(V[x], Boundary, 0, vx);
                Solver.impose(V[y], Boundary, 0, vy);
                Solver.impose(V[z], Boundary, 0, vx);
                Solver.solve_with_solver(solverPetsc, V[x], V[y], V[z]);
                Particles.ghost_get<Velocity>(SKIP_LABELLING);
                div = -(Dx(V[x]) + Dy(V[y])+Dz(V[z]));
                P = P + div;
                for (int i = 0; i < bulk.size(); i++) {
                    Particles_subset.getProp<0>(i) = Particles.template getProp<4>(bulk.template get<0>(i));
                }
                for (int j = 0; j < Boundary.size(); j++) {
                    auto p = Boundary.get<0>(j);
                    Particles.getProp<4>(p) = 0;
                }
                sum = 0;
                sum1 = 0;
                for (int j = 0; j < bulk.size(); j++) {
                    auto p = bulk.get<0>(j);
                    sum += (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) *
                           (Particles.getProp<18>(p)[0] - Particles.getProp<1>(p)[0]) +
                           (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]) *
                           (Particles.getProp<18>(p)[1] - Particles.getProp<1>(p)[1]) +
                           (Particles.getProp<18>(p)[2] - Particles.getProp<1>(p)[2]) *
                           (Particles.getProp<18>(p)[2] - Particles.getProp<1>(p)[2])

                            ;
                    sum1 += Particles.getProp<1>(p)[0] * Particles.getProp<1>(p)[0] +
                            Particles.getProp<1>(p)[1] * Particles.getProp<1>(p)[1]+
                            Particles.getProp<1>(p)[2] * Particles.getProp<1>(p)[2];
                }
                sum = sqrt(sum);
                sum1 = sqrt(sum1);

                v_cl.sum(sum);
                v_cl.sum(sum1);
                v_cl.execute();
                V_t = V;
                Particles.ghost_get<1, 4, 18>(SKIP_LABELLING);
                V_err_old = V_err;
                V_err = sum / sum1;
                if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
                    errctr++;
                    //alpha_P -= 0.1;
                } else {
                    errctr = 0;
                }
                if (n > 3) {
                    if (errctr > 3) {
                        std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN ERROR" << std::endl;
                        Vreset = 1;
                        break;
                    } else {
                        Vreset = 0;
                    }
                }
                n++;
                //Particles.write_frame("V_debug", n);
                if (v_cl.rank() == 0) {
                    std::cout << "Rel l2 cgs err in V = " << V_err << " at " << n << std::endl;
                }
            }
            tt.stop();
            u[x][x] = Dx(V[x]);
            u[x][y] = 0.5 * (Dx(V[y]) + Dy(V[x]));
            u[x][z] = 0.5 * (Dx(V[z]) + Dz(V[x]));
            u[y][x] = 0.5 * (Dy(V[x]) + Dx(V[y]));
            u[y][y] = Dy(V[y]);
            u[y][z] = 0.5 * (Dy(V[z]) + Dz(V[y]));
            u[z][x] = 0.5 * (Dz(V[x]) + Dx(V[z]));
            u[z][y] = 0.5 * (Dz(V[y]) + Dy(V[z]));
            u[z][z] = Dz(V[z]);

            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds with " << n
                          << " iterations. dt is set to " << dt
                          << std::endl;
            }

            W[x][x] = 0;
            W[x][y] = 0.5 * (Dy(V[x]) - Dx(V[y]));
            W[x][z] = 0.5 * (Dz(V[x]) - Dx(V[z]));
            W[y][x] = 0.5 * (Dx(V[y]) - Dy(V[x]));
            W[y][y] = 0;
            W[y][z] = 0.5 * (Dz(V[y]) - Dy(V[z]));
            W[z][x] = 0.5 * (Dx(V[z]) - Dz(V[x]));
            W[z][y] = 0.5 * (Dy(V[z]) - Dz(V[y]));
            W[z][z] = 0;

            Particles.deleteGhost();
            Particles.write_frame("Polar3d", ctr);
            Particles.ghost_get<0>();
            ctr++;
            auto lambda = -1/(9*gama)*(-3*h[x]*Pol[x]-3*h[y]*Pol[y]-3*h[z]*Pol[z]+
                    gama*nu*(Pol[x]*Pol[x]*u[x][x]+Pol[y]*Pol[y]*u[y][y]+Pol[z]*Pol[z]*u[z][z]+
                    Pol[x]*(Pol[y]*(u[x][y]+u[y][x])+Pol[z]*(u[x][z]+u[z][x])) +Pol[y]*Pol[z]*(u[y][z]+u[z][y])))/(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]);
            dPol=Pol;
            Particles.ghost_get<8>(SKIP_LABELLING);
            k1[x] = h[x]/gama-nu*(Pol[x]*u[x][x]+Pol[y]*u[x][y]+Pol[z]*u[x][z]) + lambda*Pol[x]/delmu + (W[x][x]*Pol[x]+W[x][y]*Pol[y]+W[x][z]*Pol[z]);
            k1[y] = h[y]/gama-nu*(Pol[x]*u[y][x]+Pol[y]*u[y][y]+Pol[z]*u[y][z]) + lambda*Pol[y]/delmu + (W[y][x]*Pol[x]+W[y][y]*Pol[y]+W[y][z]*Pol[z]);
            k1[z] = h[z]/gama-nu*(Pol[x]*u[z][x]+Pol[y]*u[z][y]+Pol[z]*u[z][z]) + lambda*Pol[z]/delmu + (W[z][x]*Pol[x]+W[z][y]*Pol[y]+W[z][z]*Pol[z]);
            Pol = dPol + (0.5 * dt) * k1;
    /*        for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
            }*/
            for (int j = 0; j < Boundary.size(); j++) {
                auto p = Boundary.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[z] = 0.0;
            }
            Particles.ghost_get<0>(SKIP_LABELLING);

            h[x]=Ks*(dxxpx + dxypy + dxzpz) +
                 Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                     py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                     px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
                 Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                     px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

            h[y]= Ks*(dxypx + dyypy + dyzpz) +
                  Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                      py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                      px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
                  Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                      px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));


            h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
                 Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                     py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                     px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
                 Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                     px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));


            k2[x] = h[x]/gama-nu*(Pol[x]*u[x][x]+Pol[y]*u[x][y]+Pol[z]*u[x][z]) + lambda*Pol[x]/delmu + (W[x][x]*Pol[x]+W[x][y]*Pol[y]+W[x][z]*Pol[z]);
            k2[y] = h[y]/gama-nu*(Pol[x]*u[y][x]+Pol[y]*u[y][y]+Pol[z]*u[y][z]) + lambda*Pol[y]/delmu + (W[y][x]*Pol[x]+W[y][y]*Pol[y]+W[y][z]*Pol[z]);
            k2[z] = h[z]/gama-nu*(Pol[x]*u[z][x]+Pol[y]*u[z][y]+Pol[z]*u[z][z]) + lambda*Pol[z]/delmu + (W[z][x]*Pol[x]+W[z][y]*Pol[y]+W[z][z]*Pol[z]);
            Pol = dPol + (0.5 * dt) * k2;
       /*     for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
            }*/
            for (int j = 0; j < Boundary.size(); j++) {
                auto p = Boundary.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[z] = 0.0;

            }


            Particles.ghost_get<0>(SKIP_LABELLING);
            h[x]=Ks*(dxxpx + dxypy + dxzpz) +
                 Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                     py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                     px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
                 Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                     px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

            h[y]= Ks*(dxypx + dyypy + dyzpz) +
                  Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                      py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                      px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
                  Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                      px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));


            h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
                 Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                     py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                     px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
                 Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                     px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));


            k3[x] = h[x]/gama-nu*(Pol[x]*u[x][x]+Pol[y]*u[x][y]+Pol[z]*u[x][z]) + lambda*Pol[x]/delmu + (W[x][x]*Pol[x]+W[x][y]*Pol[y]+W[x][z]*Pol[z]);
            k3[y] = h[y]/gama-nu*(Pol[x]*u[y][x]+Pol[y]*u[y][y]+Pol[z]*u[y][z]) + lambda*Pol[y]/delmu + (W[y][x]*Pol[x]+W[y][y]*Pol[y]+W[y][z]*Pol[z]);
            k3[z] = h[z]/gama-nu*(Pol[x]*u[z][x]+Pol[y]*u[z][y]+Pol[z]*u[z][z]) + lambda*Pol[z]/delmu + (W[z][x]*Pol[x]+W[z][y]*Pol[y]+W[z][z]*Pol[z]);
            Pol = dPol + (dt)*k3;
            /*for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
            }*/
            for (int j = 0; j < Boundary.size(); j++) {
                auto p = Boundary.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[z] = 0.0;

            }

            Particles.ghost_get<0>(SKIP_LABELLING);
            h[x]=Ks*(dxxpx + dxypy + dxzpz) +
                 Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                     py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                     px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
                 Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                     px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

            h[y]= Ks*(dxypx + dyypy + dyzpz) +
                  Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                      py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                      px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
                  Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                      px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));


            h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
                 Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                     py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                     px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
                 Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                     px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));

            k4[x] = h[x]/gama-nu*(Pol[x]*u[x][x]+Pol[y]*u[x][y]+Pol[z]*u[x][z]) + lambda*Pol[x]/delmu + (W[x][x]*Pol[x]+W[x][y]*Pol[y]+W[x][z]*Pol[z]);
            k4[y] = h[y]/gama-nu*(Pol[x]*u[y][x]+Pol[y]*u[y][y]+Pol[z]*u[y][z]) + lambda*Pol[y]/delmu + (W[y][x]*Pol[x]+W[y][y]*Pol[y]+W[y][z]*Pol[z]);
            k4[z] = h[z]/gama-nu*(Pol[x]*u[z][x]+Pol[y]*u[z][y]+Pol[z]*u[z][z]) + lambda*Pol[z]/delmu + (W[z][x]*Pol[x]+W[z][y]*Pol[y]+W[z][z]*Pol[z]);

            Pol = dPol + (dt / 6.0) * (k1 + (2.0 * k2) + (2.0 * k3) + k4);
           /* for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
            }*/
            for (int j = 0; j < Boundary.size(); j++) {
                auto p = Boundary.get<0>(j);
                Particles.getProp<0>(p) = (Particles.getProp<0>(p) == 0) ? 1 : Particles.getProp<0>(p);
                Particles.getProp<0>(p)[x] = sin(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[y] = cos(2 * M_PI * (cos((2 * Particles.getPos(p)[x] - Lx) / Lx) -
                                                             sin((2 * Particles.getPos(p)[y] - Ly) / Ly)));
                Particles.getProp<0>(p)[z] = 0.0;
            }

            k1 = V;
            k2 = 0.5 * dt * k1 + V;
            k3 = 0.5 * dt * k2 + V;
            k4 = dt * k3 + V;
            Pos = Pos + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);

            Particles.map();

            Particles.ghost_get<0, ExtForce, 27>();
            indexUpdate(Particles,Particles_subset,Boundary, bulk, up, down, left,right,front,back);
            Particles_subset.map();
            Particles_subset.ghost_get<0>();

            //Particles_subset.write("debug");

            tt.start();
            Dx.update(Particles);
            Dy.update(Particles);
            Dz.update(Particles);
            Dxy.update(Particles);
            Dxz.update(Particles);
            Dyz.update(Particles);
            auto Dyx = Dxy;
            auto Dzy = Dyz;
            auto Dzx = Dxz;
            Dxx.update(Particles);
            Dyy.update(Particles);
            Dzz.update(Particles);

            Bulk_Dx.update(Particles_subset);
            Bulk_Dy.update(Particles_subset);
            Bulk_Dz.update(Particles_subset);


            tt.stop();
            if (v_cl.rank() == 0) {
                std::cout << "Updation of operators took " << tt.getwct() << " seconds." << std::endl;
                std::cout << "Time step " << ctr - 1 << " : " << tim << " over." << std::endl;
                std::cout << "----------------------------------------------------------" << std::endl;
            }

            tim += dt;
        }

        Particles.deleteGhost();
        Particles.write("Polar_Last");

        Dx.deallocate(Particles);
        Dy.deallocate(Particles);
        Dz.deallocate(Particles);
        Dxy.deallocate(Particles);
        Dxz.deallocate(Particles);
        Dyz.deallocate(Particles);
        Dxx.deallocate(Particles);
        Dyy.deallocate(Particles);
        Dzz.deallocate(Particles);
        Bulk_Dx.deallocate(Particles_subset);
        Bulk_Dy.deallocate(Particles_subset);
        Bulk_Dz.deallocate(Particles_subset);
        Particles.deleteGhost();
        tt2.stop();
        if (v_cl.rank() == 0) {
            std::cout << "The simulation took " << tt2.getcputime() << "(CPU) ------ " << tt2.getwct()
                      << "(Wall) Seconds.";
        }
    }


BOOST_AUTO_TEST_SUITE_END()