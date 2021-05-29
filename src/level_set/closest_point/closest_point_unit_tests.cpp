/* Unit tests for the closest point methods
 * Date : 29 May 2021
 * Author : sachin
 */

#include<iostream>
#include <boost/test/unit_test_log.hpp>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "VCluster/VCluster.hpp"
#include "Vector/Vector.hpp"
#include "FiniteDifference/util/common.hpp"
#include "FiniteDifference/eq.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "closest_point.hpp"


constexpr int narrow_band_half_width = 8;

typedef struct EllipseParameters{
    double origin[3];
    double radiusA;
    double radiusB;
    double radiusC;
    double eccentricity;
} EllipseParams;

// Generate an ellipsoid initial levelset signed distance function
template<typename grid_type, typename domain_type, const unsigned int phi_field>
void initializeLSEllipsoid(grid_type &gd, const domain_type &domain,  const EllipseParams &params)
{
    auto it = gd.getDomainIterator();
    double dx = gd.getSpacing()[0];
    double dy = gd.getSpacing()[1];
    double dz = gd.getSpacing()[2];
    while(it.isNext())
    {
        auto key = it.get();
        auto key_g = gd.getGKey(key);
        
        double posx = key_g.get(0)*dx + domain.getLow(0);
        double posy = key_g.get(1)*dy + domain.getLow(1);
        double posz = key_g.get(2)*dz + domain.getLow(2);
        
        // NOTE: This is not the SDF. It is just an implicit function whose zero contour is an ellipsoid.
        double phi_val = 1.0 - sqrt(((posx - params.origin[0])/params.radiusA)*((posx - params.origin[0])/params.radiusA) + ((posy - params.origin[1])/params.radiusB)*((posy - params.origin[1])/params.radiusB) + ((posz - params.origin[2])/params.radiusC)*((posz - params.origin[2])/params.radiusC));
        gd.template get<phi_field>(key) = phi_val;
        ++it;
    }
}

BOOST_AUTO_TEST_SUITE( closest_point_test )


BOOST_AUTO_TEST_CASE( closest_point_unit_sphere )
{

    constexpr int SIM_DIM = 3;
    constexpr int POLY_ORDER = 5;
    constexpr int SIM_GRID_SIZE = 128;
    constexpr double PI = 3.141592653589793;

    // Fields - phi, cp
    using GridDist = grid_dist_id<SIM_DIM,double,aggregate<double,double[SIM_DIM]>>;
    using GridKey = grid_dist_key_dx<SIM_DIM>;

    // Grid size on each dimension
    const long int sz[SIM_DIM] = {SIM_GRID_SIZE, SIM_GRID_SIZE, SIM_GRID_SIZE};
    const size_t szu[SIM_DIM] = {(size_t) sz[0], (size_t) sz[1], (size_t) sz[2]};

    // 3D physical domain
    Box<SIM_DIM,double> domain({-1.5,-1.5,-1.5},{1.5,1.5,1.5});

    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int z = 2;

    // Alias for properties on the grid
    constexpr int phi = 0;
    constexpr int cp = 1;
    
    double nb_gamma = 0.0;

    // Create VCluster
    Vcluster<> &v_cl = create_vcluster();

	periodicity<SIM_DIM> grid_bc = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};
	// Ghost in grid units
	Ghost <SIM_DIM, long int> grid_ghost(2*narrow_band_half_width);
	GridDist gdist(szu, domain, grid_ghost, grid_bc);

    EllipseParams params;
    params.origin[x] = 0.0;
    params.origin[y] = 0.0;
    params.origin[z] = 0.0;
    params.radiusA = 1.0;
    params.radiusB = 1.0;
    params.radiusC = 1.0;

    initializeLSEllipsoid<GridDist, Box<SIM_DIM,double>, phi>(gdist, domain, params);
    nb_gamma = narrow_band_half_width * gdist.spacing(0);

    gdist.template ghost_get<phi>();

    estimateClosestPoint<GridDist, GridKey, Box<SIM_DIM,double>, SIM_DIM, POLY_ORDER, phi, cp>(gdist, domain, nb_gamma);
    gdist.template ghost_get<cp>();

    // Estimate error in closest point estimation
    auto &patches = gdist.getLocalGridsInfo();
    double max_error_cp = -1.0;
    for(int i = 0; i < patches.size();i++)
    {
        auto p_xlo = patches.get(i).Dbox.getLow(0) + patches.get(i).origin[0];
        auto p_xhi = patches.get(i).Dbox.getHigh(0) + patches.get(i).origin[0];
        auto p_ylo = patches.get(i).Dbox.getLow(1) + patches.get(i).origin[1];
        auto p_yhi = patches.get(i).Dbox.getHigh(1) + patches.get(i).origin[1];
        auto p_zlo = patches.get(i).Dbox.getLow(2) + patches.get(i).origin[2];
        auto p_zhi = patches.get(i).Dbox.getHigh(2) + patches.get(i).origin[2];

        auto it = gdist.getSubDomainIterator({p_xlo, p_ylo, p_zlo}, {p_xhi, p_yhi, p_zhi});
        while(it.isNext())
        {
            auto key = it.get();

            if(std::abs(gdist.template get<phi>(key)) < nb_gamma)
            {
                auto key_g = gdist.getGKey(key);

                double cpx = gdist.template get<cp>(key)[x];
                double cpy = gdist.template get<cp>(key)[y];
                double cpz = gdist.template get<cp>(key)[z];

                double estim_px = domain.getLow(x) + (p_xlo - algoim_padding)*gdist.spacing(x) + cpx;
                double estim_py = domain.getLow(y) + (p_ylo - algoim_padding)*gdist.spacing(y) + cpy;
                double estim_pz = domain.getLow(z) + (p_zlo - algoim_padding)*gdist.spacing(z) + cpz;
        
                double posx = key_g.get(0)*gdist.spacing(0) + domain.getLow(0);
                double posy = key_g.get(1)*gdist.spacing(1) + domain.getLow(1);
                double posz = key_g.get(2)*gdist.spacing(2) + domain.getLow(2);
                    
                double norm = sqrt(posx*posx + posy*posy + posz*posz);

                double exact_px = posx / norm;
                double exact_py = posy / norm;
                double exact_pz = posz / norm;

                max_error_cp = std::max({std::abs(estim_px - exact_px), std::abs(estim_py - exact_py), std::abs(estim_pz - exact_pz), max_error_cp});
            }
            ++it;
        }
    }
    std::cout<<"Closest_point error : " + std::to_string(max_error_cp)<<std::endl;
    double tolerance = 1e-5;
    bool check;
    if (std::abs(max_error_cp) < tolerance)
        check = true;
    else
        check = false;
    
    BOOST_TEST( check );

}


BOOST_AUTO_TEST_SUITE_END()