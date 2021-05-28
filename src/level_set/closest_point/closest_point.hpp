//
// Created by sachin on 9/9/20.
// Contains the closest point estimation and a few other level set helper functions that 
// are based on Algoim framework (for stencil interpolation etc.) to be used in OpenFPM
//

#ifndef __CLOSEST_POINT_HPP__
#define __CLOSEST_POINT_HPP__

#include "algoim_hocp.hpp"

constexpr int algoim_padding = 3;

template<typename GridType, typename GridKeyType, const unsigned int DIM, const unsigned int wrapping_field>
struct AlgoimWrapper
{
    GridType &gd;
    int patch_id;
    AlgoimWrapper(GridType& ls_phi, const int pid) : gd(ls_phi), patch_id(pid) {}
    double operator() (const blitz::TinyVector<int,DIM> idx) const
    {
        long int local_key[DIM];
        
        auto ghost_offset = gd.getLocalGridsInfo().get(patch_id).Dbox.getKP1();
        for (int dim = 0; dim < DIM; ++dim)
            local_key[dim] = idx(dim) - algoim_padding;

        // Generate OpenFPM grid_key object from local grid indices
        GridKeyType grid_key(patch_id, grid_key_dx<DIM> (local_key) + ghost_offset);
        
        return gd.template get<wrapping_field>(grid_key);
    }
};

// Compute the closest point in the local patch
template<typename GridType, typename GridKeyType, typename DomainType, const unsigned int DIM, const unsigned int ORDER, const unsigned int phi_field, const unsigned int cp_field> //, const unsigned int cp_error_field>
void estimateClosestPoint3D(GridType &gd, DomainType &domain, const double nb_gamma)
{
    using Poly = typename Algoim::StencilPoly<DIM,ORDER>::T_Poly;

    blitz::TinyVector<double,DIM> dx = {gd.spacing(0), gd.spacing(1), gd.spacing(2)};
    auto &patches = gd.getLocalGridsInfo();
    Vcluster<> &v_cl = create_vcluster();

    for(int i = 0; i < patches.size();i++)
    {
        auto p_xlo = patches.get(i).Dbox.getLow(0) + patches.get(i).origin[0];
        auto p_xhi = patches.get(i).Dbox.getHigh(0) + patches.get(i).origin[0];
        auto p_ylo = patches.get(i).Dbox.getLow(1) + patches.get(i).origin[1];
        auto p_yhi = patches.get(i).Dbox.getHigh(1) + patches.get(i).origin[1];
        auto p_zlo = patches.get(i).Dbox.getLow(2) + patches.get(i).origin[2];
        auto p_zhi = patches.get(i).Dbox.getHigh(2) + patches.get(i).origin[2];

        AlgoimWrapper<GridType, GridKeyType, DIM, phi_field> phiwrap(gd, i);

        // Find all cells containing the interface and construct the high-order polynomials
        std::vector<Algoim::detail::CellPoly<DIM,Poly>> cells;

        blitz::TinyVector<int,DIM> ext = {static_cast<int>(p_xhi - p_xlo + 1 + 2*algoim_padding), static_cast<int>(p_yhi - p_ylo + 1 + 2*algoim_padding), static_cast<int>(p_zhi - p_zlo + 1 + 2*algoim_padding)};

        Algoim::detail::createCellPolynomials(ext, phiwrap, dx, false, cells);

        std::vector<blitz::TinyVector<double,DIM>> points;
        std::vector<int> pointcells;
        Algoim::detail::samplePolynomials<DIM,Poly>(cells, 2, dx, 0.0, points, pointcells);

        Algoim::KDTree<double,DIM> kdtree(points);

        // Pass everything to the closest point computation engine
        Algoim::ComputeHighOrderCP<DIM,Poly> hocp(nb_gamma < std::numeric_limits<double>::max() ? nb_gamma*nb_gamma : std::numeric_limits<double>::max(), // squared bandradius
                                        0.5*blitz::max(dx), // amount that each polynomial overlaps / size of the bounding ball in Newton's method
                                        Algoim::sqr(std::max(1.0e-14, std::pow(blitz::max(dx), Poly::order))), // tolerance to determine convergence
                                        cells, kdtree, points, pointcells, dx, 0.0);

        auto it = gd.getSubDomainIterator({p_xlo, p_ylo, p_zlo},{p_xhi, p_yhi, p_zhi});
        while(it.isNext())
        {
            auto key = it.get();
            if(std::abs(gd.template get<phi_field>(key)) <= nb_gamma)
            {
                auto key_g = gd.getGKey(key);
                // NOTE: This is not the real grid coordinates, but internal coordinates for algoim
                double patch_posx = (key_g.get(0) - p_xlo + algoim_padding) * gd.spacing(0);
                double patch_posy = (key_g.get(1) - p_ylo + algoim_padding) * gd.spacing(1);
                double patch_posz = (key_g.get(2) - p_zlo + algoim_padding) * gd.spacing(2);
                blitz::TinyVector<double,DIM> patch_pos = {patch_posx, patch_posy, patch_posz}, cp;

                if (hocp.compute(patch_pos, cp))
                {
                    gd.template get<cp_field>(key)[0] = cp(0);
                    gd.template get<cp_field>(key)[1] = cp(1);
                    gd.template get<cp_field>(key)[2] = cp(2);
                }
            }
            ++it;
        }
    }
    return;
}

template<typename GridType, typename GridKeyType, int DIM, int ORDER, const unsigned int cp_field, const unsigned int extend_field, const unsigned int extend_field_temp, const unsigned int phi_field>
void extendLSFields3D(GridType &gd, const double nb_gamma)
{
    using Poly = typename Algoim::StencilPoly<DIM,ORDER>::T_Poly;
    auto &patches = gd.getLocalGridsInfo();
    blitz::TinyVector<double,DIM> dx = {gd.spacing(0), gd.spacing(1), gd.spacing(2)};

    for(int i = 0; i < patches.size();i++)
    {
        auto p_xlo = patches.get(i).Dbox.getLow(0) + patches.get(i).origin[0];
        auto p_xhi = patches.get(i).Dbox.getHigh(0) + patches.get(i).origin[0];
        auto p_ylo = patches.get(i).Dbox.getLow(1) + patches.get(i).origin[1];
        auto p_yhi = patches.get(i).Dbox.getHigh(1) + patches.get(i).origin[1];
        auto p_zlo = patches.get(i).Dbox.getLow(2) + patches.get(i).origin[2];
        auto p_zhi = patches.get(i).Dbox.getHigh(2) + patches.get(i).origin[2];

        auto it = gd.getSubDomainIterator({p_xlo, p_ylo, p_zlo},{p_xhi, p_yhi, p_zhi});

        while(it.isNext())
        {
            auto key = it.get();
            if(std::abs(gd.template get<phi_field>(key)) < nb_gamma)
            {
                double cpx = gd.template get<cp_field>(key)[0];
                double cpy = gd.template get<cp_field>(key)[1];
                double cpz = gd.template get<cp_field>(key)[2];

                blitz::TinyVector<int,DIM> coord = {static_cast<int>(floor(cpx / gd.spacing(0))), static_cast<int>(floor(cpy / gd.spacing(1))), static_cast<int>(floor(cpz / gd.spacing(2)))};
                blitz::TinyVector<double,DIM> pos = {cpx - coord(0)*gd.spacing(0), cpy - coord(1)*gd.spacing(1), cpz - coord(2)*gd.spacing(2)};   
                
                AlgoimWrapper<GridType, GridKeyType, DIM, extend_field> fieldwrap(gd,i);
                Poly field_poly = Poly(coord, fieldwrap, dx);
                gd.template get<extend_field_temp>(key) = field_poly(pos);
            }
            ++it;
        }
    }
    
    // Copy the results to the actual variable
    auto it = gd.getDomainIterator();
    while(it.isNext())
    {
        auto key = it.get();
        if(std::abs(gd.template get<phi_field>(key)) < nb_gamma)
            gd.template get<extend_field>(key) = gd.template get<extend_field_temp>(key);
        ++it;
    }
}

template<typename GridType, typename GridKeyType, int DIM, int ORDER, const unsigned int cp_field, const unsigned int phi_field>
void reinitializeLS3D(GridType &gd, const double nb_gamma)
{
    using Poly = typename Algoim::StencilPoly<DIM,ORDER>::T_Poly;
    auto &patches = gd.getLocalGridsInfo();
    blitz::TinyVector<double,DIM> dx = {gd.spacing(0), gd.spacing(1), gd.spacing(2)};

    for(int i = 0; i < patches.size();i++)
    {
        auto p_xlo = patches.get(i).Dbox.getLow(0) + patches.get(i).origin[0];
        auto p_xhi = patches.get(i).Dbox.getHigh(0) + patches.get(i).origin[0];
        auto p_ylo = patches.get(i).Dbox.getLow(1) + patches.get(i).origin[1];
        auto p_yhi = patches.get(i).Dbox.getHigh(1) + patches.get(i).origin[1];
        auto p_zlo = patches.get(i).Dbox.getLow(2) + patches.get(i).origin[2];
        auto p_zhi = patches.get(i).Dbox.getHigh(2) + patches.get(i).origin[2];

        auto it = gd.getSubDomainIterator({p_xlo, p_ylo, p_zlo},{p_xhi, p_yhi, p_zhi});

        while(it.isNext())
        {
            auto key = it.get();
            if(std::abs(gd.template get<phi_field>(key)) < nb_gamma)
            {
                double sign_fn = (gd.template get<phi_field>(key) >= 0.0)?1.0:-1.0;
                auto key_g = gd.getGKey(key);
                // NOTE: This is not the real grid coordinates, but internal coordinates for algoim
                double patch_posx = (key_g.get(0) - p_xlo + algoim_padding) * gd.spacing(0);
                double patch_posy = (key_g.get(1) - p_ylo + algoim_padding) * gd.spacing(1);
                double patch_posz = (key_g.get(2) - p_zlo + algoim_padding) * gd.spacing(2);
                
                double cpx = gd.template get<cp_field>(key)[0];
                double cpy = gd.template get<cp_field>(key)[1];
                double cpz = gd.template get<cp_field>(key)[2];

                gd.template get<phi_field>(key) = sign_fn*sqrt((patch_posx - cpx)*(patch_posx - cpx) + (patch_posy - cpy)*(patch_posy - cpy) + (patch_posz - cpz)*(patch_posz - cpz));
            }
            ++it;
        }
    }
}

#endif //__CLOSEST_POINT_HPP__
