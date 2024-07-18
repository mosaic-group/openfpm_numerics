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


template<unsigned int dim, typename T, typename particles_type, typename CellList_type, typename supportSize_type>
__global__ void gatherSupportSize_gpu(
	particles_type particles,
	CellList_type cellList,
	supportSize_type supportSize,
	T rCut)
{
	auto p = GET_PARTICLE(particles);
	Point<dim, T> pos = particles.getPos(p);

	size_t N = 0;
	auto Np = cellList.getNNIteratorBox(cellList.getCell(pos));
	while (Np.isNext())
	{
		auto q = Np.get(); ++Np;

		if (p == q) continue;
		if (pos.distance(particles.getPos(q)) < rCut) ++N;
	}

	supportSize.get(p) = N;
}


template<unsigned int dim, typename T, typename particles_type, typename CellList_type, typename supportKey_type>
__global__ void assembleSupport_gpu(
	particles_type particles,
	CellList_type cellList,
	supportKey_type supportSize,
	supportKey_type supportKeys1D,
	T rCut)
{
	auto p = GET_PARTICLE(particles);
	Point<dim, T> pos = particles.getPos(p);

	size_t  supportKeysSize = supportSize.get(p+1)-supportSize.get(p);
	size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[supportSize.get(p)];

	size_t N = 0;
	auto Np = cellList.getNNIteratorBox(cellList.getCell(pos));
	while (Np.isNext())
	{
		auto q = Np.get(); ++Np;

		if (p == q) continue;
		if (pos.distance(particles.getPos(q)) < rCut) supportKeys[N++] = q;
	}
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

	void getSupport(
		size_t N,
		openfpm::vector_custd<size_t>& kerOffsets,
		openfpm::vector_custd<size_t>& supportKeys1D,
		size_t& maxSupport,
		size_t& supportKeysTotalN)
	{
		domain.hostToDevicePos();
		auto it = domain.getDomainIteratorGPU(512);
		typedef CellList<vector_type::dims, typename vector_type::stype, Mem_fast<CudaMemory>, shift<vector_type::dims, typename vector_type::stype>> params;
		auto NN = domain.getCellListGPU(rCut);
		domain.updateCellListGPU(NN);

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
