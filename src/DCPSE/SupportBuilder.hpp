//
// Created by tommaso on 25/03/19.
//
// Modified by Abhinav and Pietro

#ifndef OPENFPM_PDATA_SUPPORTBUILDER_HPP
#define OPENFPM_PDATA_SUPPORTBUILDER_HPP

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>
#include <utility>

enum support_options
{
    RADIUS,
    LOAD,
    ADAPTIVE,             //used in SurfaceDcpse
    AT_LEAST_N_PARTICLES  //used in regression module
};

#endif //OPENFPM_PDATA_SUPPORTBUILDER_HPP
