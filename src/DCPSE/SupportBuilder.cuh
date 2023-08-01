//
// Created by Serhii
//

#ifndef OPENFPM_PDATA_SUPPORTBUILDER_CUH
#define OPENFPM_PDATA_SUPPORTBUILDER_CUH

#include <Space/Shape/Point.hpp>
#include <Vector/vector_dist.hpp>
#include "Support.hpp"
#include <utility>
#include "SupportBuilder.hpp"


template <unsigned int dim>
__device__ __host__ bool nextCell(size_t (&offset)[dim], size_t maxOffset) {
    size_t i = 0;

    while (i < dim) {
        if ((++offset[i++])/maxOffset)
            for (size_t j = 0; j < i; ++j)
                offset[j] = 0;
        else
            return true;
    }
    return false;
}

template<unsigned int dim, typename T, typename particles_type, typename CellList_type, typename supportSize_type>
__global__ void gatherSupportSize_gpu(
    particles_type particles, CellList_type cl, supportSize_type supportSize, T rCut) {
    auto p_key = GET_PARTICLE(particles);
    Point<dim, T> pos = particles.getPos(p_key);
    auto cell = cl.getCellGrid(pos);

    size_t grSize[dim]; cl.getGrid().getSize(grSize);
    size_t offset[dim]; for (int i = 0; i < dim; ++i) offset[i] = 0;    
    grid_key_dx<dim> middle; for (int i = 0; i < dim; ++i) middle.set_d(i,1);     

    size_t N = 0;
    do {
        auto key=grid_key_dx<dim>(offset); key=cell+key-middle;

        for (size_t i = 0; i < dim; ++i)
            if (key.value(i) < 0 || key.value(i) >= grSize[i])
                continue;

        mem_id id = cl.getGrid().LinId(key);
        const size_t cellLinId = static_cast<size_t>(id);
        const size_t elemsInCell = cl.getNelements(cellLinId);

        for (size_t k = 0; k < elemsInCell; ++k) {
            size_t el = cl.get(cellLinId, k);

            if (p_key == el) continue;
            if (pos.distance(particles.getPosOrig(el)) < rCut) ++N;
        }

    } while (nextCell<dim>(offset, 2+1));

    supportSize.get(p_key) = N;
}

template<unsigned int dim, typename T, typename particles_type, typename CellList_type, typename supportKey_type>
__global__ void assembleSupport_gpu(particles_type particles, CellList_type cl, supportKey_type supportSize, supportKey_type supportKeys1D, T rCut) {
    auto p_key = GET_PARTICLE(particles);
    Point<dim, T> pos = particles.getPos(p_key);
    auto cell = cl.getCellGrid(pos);

    size_t  supportKeysSize = supportSize.get(p_key+1)-supportSize.get(p_key);
    size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[supportSize.get(p_key)];

    size_t grSize[dim]; cl.getGrid().getSize(grSize);
    size_t offset[dim]; for (int i = 0; i < dim; ++i) offset[i] = 0;    
    grid_key_dx<dim> middle; for (int i = 0; i < dim; ++i) middle.set_d(i,1);     

    size_t N = 0;
    do {
        auto key=grid_key_dx<dim>(offset); key=cell+key-middle;

        for (size_t i = 0; i < dim; ++i)
            if (key.value(i) < 0 || key.value(i) >= grSize[i])
                continue;

        mem_id id = cl.getGrid().LinId(key);
        const size_t cellLinId = static_cast<size_t>(id);
        const size_t elemsInCell = cl.getNelements(cellLinId);

        for (size_t k = 0; k < elemsInCell; ++k) {
            size_t el = cl.get(cellLinId, k);

            if (p_key == el) continue;
            if (pos.distance(particles.getPosOrig(el)) < rCut) supportKeys[N++] = el;
        }

    } while (nextCell<dim>(offset, 2+1));
}



template<typename vector_type>
class SupportBuilderGPU
{
private:
    vector_type &domain;
    typename vector_type::stype rCut;

public:
    SupportBuilderGPU(vector_type &domain, typename vector_type::stype rCut)
        : domain(domain), rCut(rCut) {}

    void getSupport(size_t N, openfpm::vector_custd<size_t>& kerOffsets, openfpm::vector_custd<size_t>& supportKeys1D, 
        size_t& maxSupport, size_t& supportKeysTotalN)
    {
        domain.hostToDevicePos();
        auto it = domain.getDomainIteratorGPU(512);
        typedef CellList_gen<vector_type::dims, typename vector_type::stype, Process_keys_lin, Mem_fast<CudaMemory>, shift<vector_type::dims, typename vector_type::stype>> params;
        // auto NN = domain.getCellListGPU(rCut);
        auto NN = domain.template getCellList<params>(rCut);
        NN.hostToDevice();


        // +1 to allow getting size from cumulative sum: "size[i+1] - size[i]"
        kerOffsets.resize(N+1);
        gatherSupportSize_gpu<vector_type::dims><<<it.wthr,it.thr>>>(domain.toKernel(), NN.toKernel(), kerOffsets.toKernel(), rCut);
        kerOffsets.template deviceToHost();

        supportKeysTotalN = 0; maxSupport = 0;

        for (size_t i = 0; i < N; ++i) {   
            size_t sz = kerOffsets.get(i);
            kerOffsets.get(i) = supportKeysTotalN;
            supportKeysTotalN += sz;
            if (maxSupport < sz) maxSupport = sz;
        }
        kerOffsets.get(N) = supportKeysTotalN;

        supportKeys1D.resize(supportKeysTotalN);
        kerOffsets.template hostToDevice();
        assembleSupport_gpu<vector_type::dims><<<it.wthr,it.thr>>>(domain.toKernel(), NN.toKernel(), kerOffsets.toKernel(), supportKeys1D.toKernel(), rCut);
        supportKeys1D.template deviceToHost();
    }
};


#endif //OPENFPM_PDATA_SUPPORTBUILDER_CUH
