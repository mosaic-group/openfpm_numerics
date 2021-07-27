//
// Created by Serhii
//
#ifndef OPENFPM_PDATA_DCPSE_CUH
#define OPENFPM_PDATA_DCPSE_CUH

#if defined(__NVCC__) && defined(HAVE_EIGEN)

#include "Vector/vector_dist.hpp"
#include "MonomialBasis.hpp"
#include "DMatrix/EMatrix.hpp"
#include "SupportBuilder.hpp"
#include "Support.hpp"
#include "Vandermonde.hpp"
#include "Vandermonde.cuh"
#include "DcpseDiagonalScalingMatrix.hpp"
#include "DcpseRhs.hpp"
#include "DcpseRhs.cuh"

#include <chrono>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>


template<unsigned int dim, typename particles_type, typename T, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename calcKernels_type>
__global__ void calcKernels_gpu(particles_type, monomialBasis_type, supportKey_type, T**, localEps_type, size_t, calcKernels_type);

template<unsigned int dim, typename T, typename particles_type, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename matrixOffset_type, typename matrix_type>
__global__ void assembleLocalMatrices_gpu( particles_type, Point<dim, unsigned int>, unsigned int, monomialBasis_type, supportKey_type, supportKey_type,
    T**, T**, localEps_type, localEps_type, matrixOffset_type, matrix_type, matrixOffset_type, matrix_type, matrix_type BMat);

template <unsigned int dim, typename T, typename monomialBasis_type>
__device__ T computeKernel_gpu(Point<dim, T>&, const T*, const monomialBasis_type&);


template<unsigned int dim, typename vector_type>
class Dcpse_gpu {
    // static_assert(std::is_floating_point<typename vector_type::stype>, "CUBLAS supports only float or double");
public:
    typedef typename vector_type::stype T;
    typedef typename vector_type::value_type part_type;
    typedef vector_type vtype;

    // This works in this way:
    // 1) User constructs this by giving a domain of points (where one of the properties is the value of our f),
    //    the signature of the differential operator and the error order bound.
    // 2) The machinery for assembling and solving the linear system for coefficients starts...
    // 3) The user can then call an evaluate(point) method to get the evaluation of the differential operator
    //    on the given point.
private:
    const Point<dim, unsigned int> differentialSignature;
    const unsigned int differentialOrder;
    MonomialBasis<dim> monomialBasis;

    openfpm::vector_custd<size_t> localSupportRefs; // Each MPI rank has just access to the local ones
    openfpm::vector<openfpm::vector<size_t>> localSupportKeys; // Each MPI rank has just access to the local ones

    openfpm::vector_custd<T> localEps; // Each MPI rank has just access to the local ones
    openfpm::vector_custd<T> localEpsInvPow; // Each MPI rank has just access to the local ones

    openfpm::vector<size_t> kerOffsets;
    openfpm::vector_custd<T> calcKernels;

    vector_type & particles;
    double rCut;
    unsigned int convergenceOrder;
    double supportSizeFactor;

    support_options opt;

public:

    // Here we require the first element of the aggregate to be:
    // 1) the value of the function f on the point
    Dcpse_gpu(vector_type &particles,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1,
          support_options opt = support_options::N_PARTICLES)
        :particles(particles),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            opt(opt)
    {
        // This 
        particles.ghost_get_subset();
        if (supportSizeFactor < 1) 
        {
            initializeAdaptive(particles, convergenceOrder, rCut);
        } 
        else 
        {
            initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
        }
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k)
    {
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        auto & keys = localSupportKeys.get(k);
        for (int i = 0; i < keys.size(); i++)
        {
            size_t xqK = keys.get(i);
            particles.template getProp<prp>(xqK) += calcKernels.get(kerOff+i);
        }
    }

    template<unsigned int prp>
    void DrawKernelNN(vector_type &particles, int k)
    {
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        auto & keys = localSupportKeys.get(k);
        for (int i = 0; i < keys.size(); i++)
        {
            size_t xqK = keys.get(i);
            particles.template getProp<prp>(xqK) = 1.0;
        }
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k, int i)
    {
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        auto & keys = localSupportKeys.get(k);
        for (int i = 0; i < keys.size(); i++)
        {
            size_t xqK = keys.get(i);
            particles.template getProp<prp>(xqK)[i] += calcKernels.get(kerOff+i);
        }
    }

    void checkMomenta(vector_type &particles)
    {
        openfpm::vector<aggregate<double,double>> momenta;
        openfpm::vector<aggregate<double,double>> momenta_accu;

        momenta.resize(monomialBasis.size());
        momenta_accu.resize(monomialBasis.size());

        for (int i = 0; i < momenta.size(); i++)
        {
            momenta.template get<0>(i) =  3000000000.0;
            momenta.template get<1>(i) = -3000000000.0;
        }

        size_t N = particles.size_local_orig();
        for (size_t j = 0; j < N; ++j)
        {
            double eps = localEps.get(j);

            for (int i = 0; i < momenta.size(); i++)
            {
                momenta_accu.template get<0>(i) =  0.0;
            }

            size_t xpK = localSupportRefs.get(j);
            Point<dim, T> xp = particles.getPos(xpK);
            size_t kerOff = kerOffsets.get(xpK);
            auto & keys = localSupportKeys.get(j);

            for (int i = 0; i < keys.size(); i++)
            {
                size_t xqK = keys.get(i);
                Point<dim, T> xq = particles.getPosOrig(xqK);
                Point<dim, T> normalizedArg = (xp - xq) / eps;

                auto ker = calcKernels.get(kerOff+i);

                int counter = 0;
                size_t N = monomialBasis.getElements().size();

                for (size_t i = 0; i < N; ++i)
                {
                    const Monomial<dim> &m = monomialBasis.getElement(i);

                    T mbValue = m.evaluate(normalizedArg);
                    momenta_accu.template get<0>(counter) += mbValue * ker;

                    ++counter;
                }
            }

            for (int i = 0; i < momenta.size(); i++)
            {
                if (momenta_accu.template get<0>(i) < momenta.template get<0>(i))
                {
                    momenta.template get<0>(i) = momenta_accu.template get<0>(i);
                }

                if (momenta_accu.template get<1>(i) > momenta.template get<1>(i))
                {
                    momenta.template get<1>(i) = momenta_accu.template get<0>(i);
                }
            }
        }

        for (int i = 0; i < momenta.size(); i++)
        {
            std::cout << "MOMENTA: " << monomialBasis.getElements()[i] << "Min: " << momenta.template get<0>(i) << "  " << "Max: " << momenta.template get<1>(i) << std::endl;
        }
    }

    /**
     * Computes the value of the differential operator on all the particles,
     * using the f values stored at the fValuePos position in the aggregate
     * and storing the resulting Df values at the DfValuePos position in the aggregate.
     * @tparam fValuePos Position in the aggregate of the f values to use.
     * @tparam DfValuePos Position in the aggregate of the Df values to store.
     * @param particles The set of particles to iterate over.
     */
    template<unsigned int fValuePos, unsigned int DfValuePos>
    void computeDifferentialOperator(vector_type &particles) {
        char sign = 1;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        size_t N = particles.size_local_orig();
        for (size_t j = 0; j < N; ++j)
        {
            double epsInvPow = localEpsInvPow.get(j);

            T Dfxp = 0;
            size_t xpK = localSupportRefs.get(j);
            Point<dim, typename vector_type::stype> xp = particles.getPos(xpK);
            T fxp = sign * particles.template getProp<fValuePos>(xpK);
            size_t kerOff = kerOffsets.get(xpK);
            auto & keys = localSupportKeys.get(j);

            for (int i = 0; i < keys.size(); i++)
            {
                size_t xqK = keys.get(i);
                T fxq = particles.template getProp<fValuePos>(xqK);

                Dfxp += (fxq + fxp) * calcKernels.get(kerOff+i);
            }
            Dfxp *= epsInvPow;
       
            particles.template getProp<DfValuePos>(xpK) = Dfxp;
        }
    }


    /*! \brief Get the number of neighbours
     *
     * \return the number of neighbours
     *
     */
    inline int getNumNN(const vect_dist_key_dx &key)
    {
        return localSupportKeys.get(key.getKey()).size();
    }

    /*! \brief Get the coefficent j (Neighbour) of the particle key
     *
     * \param key particle
     * \param j neighbour
     *
     * \return the coefficent
     *
     */
    inline T getCoeffNN(const vect_dist_key_dx &key, int j)
    {
        size_t base = kerOffsets.get(key.getKey());
        return calcKernels.get(base + j);
    }

    /*! \brief Get the number of neighbours
     *
     * \return the number of neighbours
     *
     */
    inline size_t getIndexNN(const vect_dist_key_dx &key, int j)
    {
        return localSupportKeys.get(key.getKey()).get(j);
    }


    inline T getSign()
    {
        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        return sign;
    }

    T getEpsilonInvPrefactor(const vect_dist_key_dx &key)
    {
        return localEpsInvPow.get(key.getKey());
    }

    /**
     * Computes the value of the differential operator for one particle for o1 representing a scalar
     *
     * \param key particle
     * \param o1 source property
     * \return the selected derivative
     *
     */
    template<typename op_type>
    auto computeDifferentialOperator(const vect_dist_key_dx &key,
                                     op_type &o1) -> decltype(is_scalar<std::is_fundamental<decltype(o1.value(
            key))>::value>::analyze(key, o1)) {

        typedef decltype(is_scalar<std::is_fundamental<decltype(o1.value(key))>::value>::analyze(key, o1)) expr_type;

        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        double eps = localEps.get(key.getKey());
        double epsInvPow = localEpsInvPow.get(key.getKey());

        auto &particles = o1.getVector();

        expr_type Dfxp = 0;
        size_t xpK = localSupportRefs.get(key.getKey());
        Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key);
        size_t kerOff = kerOffsets.get(xpK);
        auto & keys = localSupportKeys.get(key.getKey());
        for (int i = 0; i < keys.size(); i++)
        {
            size_t xqK = keys.get(i);
            expr_type fxq = o1.value(vect_dist_key_dx(xqK));
            Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+i);
        }
        Dfxp = Dfxp * epsInvPow;
        //
        //T trueDfxp = particles.template getProp<2>(xpK);
        // Store Dfxp in the right position
        return Dfxp;
    }

    /**
     * Computes the value of the differential operator for one particle for o1 representing a vector
     *
     * \param key particle
     * \param o1 source property
     * \param i component
     * \return the selected derivative
     *
     */
    template<typename op_type>
    auto computeDifferentialOperator(const vect_dist_key_dx &key,
                                     op_type &o1,
                                     int i) -> typename decltype(is_scalar<std::is_fundamental<decltype(o1.value(
            key))>::value>::analyze(key, o1))::coord_type {

        typedef typename decltype(is_scalar<std::is_fundamental<decltype(o1.value(key))>::value>::analyze(key, o1))::coord_type expr_type;

        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        double eps = localEps.get(key.getKey());
        double epsInvPow = localEpsInvPow(key.getKey());

        auto &particles = o1.getVector();

        expr_type Dfxp = 0;
        size_t xpK = localSupportRefs.get(key.getKey());

        Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key)[i];
        size_t kerOff = kerOffsets.get(xpK);
        auto & keys = localSupportKeys.get(key.getKey());
        for (int j = 0 ; j < keys.size() ; j++)
        {
            size_t xqK = keys.get(j);
            expr_type fxq = o1.value(vect_dist_key_dx(xqK))[i];
            Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+j);
        }
        Dfxp = Dfxp * epsInvPow;
        //
        //T trueDfxp = particles.template getProp<2>(xpK);
        // Store Dfxp in the right position
        return Dfxp;
    }

    void initializeUpdate(vector_type &particles)
    {
        localSupportKeys.clear();
        localSupportRefs.clear();
        localEps.clear();
        localEpsInvPow.clear();
        calcKernels.clear();
        kerOffsets.clear();

        initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
    }

private:

    void initializeAdaptive(vector_type &particles,
                            unsigned int convergenceOrder,
                            T rCut) {
        // Still need to be tested
        SupportBuilder<vector_type>
                supportBuilder(particles, differentialSignature, rCut);
        unsigned int requiredSupportSize = monomialBasis.size();

        localSupportKeys.resize(particles.size_local_orig());
        localSupportRefs.resize(particles.size_local_orig());
        localEps.resize(particles.size_local_orig());
        localEpsInvPow.resize(particles.size_local_orig());
        kerOffsets.resize(particles.size_local_orig());
        kerOffsets.fill(-1);

        size_t offsetKer = 0;
        auto it = particles.getDomainIterator();
        while (it.isNext()) {
            const T condVTOL = 1e2;

            // Get the points in the support of the DCPSE kernel and store the support for reuse
            Support support = supportBuilder.getSupport(it, requiredSupportSize,opt);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());

            // Vandermonde matrix computation
            Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                    vandermonde(support, monomialBasis, particles);
            vandermonde.getMatrix(V);

            T condV = conditionNumber(V, condVTOL);
            T eps = vandermonde.getEps();

            if (condV > condVTOL) {
                requiredSupportSize *= 2;
                std::cout
                        << "INFO: Increasing, requiredSupportSize = " << requiredSupportSize
                        << std::endl; // debug
                continue;
            } else {
                requiredSupportSize = monomialBasis.size();
            }

            auto key_o = particles.getOriginKey(it.get());

            localSupportRefs.get(key_o.getKey()) = support.getReferencePointKey();
            localSupportKeys.get(key_o.getKey()) = support.getKeys();
            kerOffsets.get(key_o.getKey()) = offsetKer;

            offsetKer += support.getKeys().size();
            ++it;
        }

        // move monomial basis to kernel
        auto& basis = monomialBasis.getBasis();
        openfpm::vector_custd<Monomial_gpu<dim>> basisTemp(basis.begin(), basis.end());
        basisTemp.template hostToDevice();
        MonomialBasis<dim, aggregate<Monomial_gpu<dim>>, openfpm::vector_custd_ker, memory_traits_inte> monomialBasisKernel(basisTemp.toKernel());

        size_t numMatrices = localSupportKeys.size();
        size_t monomialBasisSize = monomialBasis.size();
        size_t localSupportRefKeysTotalN = 0;
        size_t EMatN = 0, VMatN = 0;

        for (size_t i = 0; i < numMatrices; ++i) {
            size_t localSupportSize = localSupportKeys.get(i).size();
            localSupportRefKeysTotalN += localSupportSize;
            EMatN += localSupportSize * localSupportSize;
            VMatN += localSupportSize * monomialBasisSize;
        }

        openfpm::vector_custd<size_t> localSupportKeys1D(localSupportRefKeysTotalN+numMatrices+1);
        openfpm::vector_custd<size_t> EMatOffsets(numMatrices+1);
        openfpm::vector_custd<size_t> VMatOffsets(numMatrices+1);
        openfpm::vector_custd<T> EMat(EMatN);
        openfpm::vector_custd<T> VMat(VMatN);
        // B has the same dimensions as V
        openfpm::vector_custd<T> BMat(VMatN);

        // localSupportKeys1D is populated:
        // offsets of keys first, followed by keys 
        size_t offset = numMatrices+1;
        for (size_t i = 0; i < numMatrices; ++i) {
            localSupportKeys1D.get(i) = offset;
            offset += localSupportKeys.get(i).size();
        }
        localSupportKeys1D.get(numMatrices) = offset;

        offset = numMatrices+1;
        for (size_t i = 0; i < numMatrices; ++i)
            for (size_t j = 0; j < localSupportKeys.get(i).size(); ++j, ++offset)
                localSupportKeys1D.get(offset) = localSupportKeys.get(i).get(j);

        // fill the matrices' offsets for E, V (B has the same dim as V)
        size_t offsetE = 0, offsetV = 0;
        for (size_t i = 0; i < numMatrices; ++i) {
            EMatOffsets.get(i) = offsetE;
            VMatOffsets.get(i) = offsetV;
            offsetE += localSupportKeys.get(i).size() * localSupportKeys.get(i).size();
            offsetV += localSupportKeys.get(i).size() * monomialBasisSize;
        }

        // allocate device space for A, b on device
        T *d_A; cudaMalloc((void**)&d_A, sizeof(T)*numMatrices*monomialBasisSize*monomialBasisSize);
        T *d_b; cudaMalloc((void**)&d_b, sizeof(T)*numMatrices*monomialBasisSize);

        // create array of pointers to pass T** pointers to cublas subroutines
        T **h_A_pointers = (T**)malloc(numMatrices*sizeof(T*));
        for (int i = 0; i < numMatrices; i++) h_A_pointers[i] = d_A + i*monomialBasisSize*monomialBasisSize;

        T **h_b_pointers = (T**)malloc(numMatrices*sizeof(T*));
        for (int i = 0; i < numMatrices; i++) h_b_pointers[i] = d_b + i*monomialBasisSize;

        T **d_A_pointers; cudaMalloc((void**)&d_A_pointers, numMatrices*sizeof(T*));
        cudaMemcpy(d_A_pointers, h_A_pointers, numMatrices*sizeof(T*), cudaMemcpyHostToDevice);

        T **d_b_pointers; cudaMalloc((void**)&d_b_pointers, numMatrices*sizeof(T*));
        cudaMemcpy(d_b_pointers, h_b_pointers, numMatrices*sizeof(T*), cudaMemcpyHostToDevice);

        // assemble local matrices on GPU
        particles.hostToDevicePos();
        localSupportKeys1D.template hostToDevice();
        EMatOffsets.template hostToDevice();
        VMatOffsets.template hostToDevice();
        localSupportRefs.template hostToDevice();

        auto it2 = particles.getDomainIteratorGPU(512);
        assembleLocalMatrices_gpu<<<it2.wthr,it2.thr>>>(particles.toKernel(), differentialSignature, differentialOrder, monomialBasisKernel, localSupportRefs.toKernel(), localSupportKeys1D.toKernel(), 
            d_A_pointers, d_b_pointers, localEps.toKernel(), localEpsInvPow.toKernel(), EMatOffsets.toKernel(), EMat.toKernel(), VMatOffsets.toKernel(), VMat.toKernel(), BMat.toKernel());

        localEps.template deviceToHost();
        localEpsInvPow.template deviceToHost();

        //cublas lu solver
        int *d_infoArray; cudaMalloc((void**)&d_infoArray,  numMatrices*sizeof(int));
        int *h_infoArray = (int *)malloc(numMatrices*sizeof(int));

        cublasHandle_t cublas_handle; cublasCreate_v2(&cublas_handle);

        // if (std::is_same<T, float>::value)
            // cublasSgetrfBatched(cublas_handle, monomialBasisSize, d_A_pointers, monomialBasisSize, NULL, d_infoArray, numMatrices);
        // else
            cublasDgetrfBatched(cublas_handle, monomialBasisSize, d_A_pointers, monomialBasisSize, NULL, d_infoArray, numMatrices);
        cudaDeviceSynchronize();

        cudaMemcpy(h_infoArray, d_infoArray, numMatrices*sizeof(int), cudaMemcpyDeviceToHost);
        for (size_t i = 0; i < numMatrices; i++)
            if (h_infoArray[i] != 0) fprintf(stderr, "Factorization of matrix %d Failed: Matrix may be singular\n", i);

        const double alpha = 1.f;
        // if (std::is_same<T, float>::value) {
        //     cublasStrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
        //     cublasStrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
        // } else {
            cublasDtrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
            cublasDtrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
        // }
        cudaDeviceSynchronize();

        // populate the calcKernels on GPU
        calcKernels.resize(localSupportRefKeysTotalN);
        localEps.template hostToDevice();
        calcKernels_gpu<dim><<<it2.wthr,it2.thr>>>(particles.toKernel(), monomialBasisKernel, localSupportKeys1D.toKernel(), d_b_pointers, localEps.toKernel(), numMatrices, calcKernels.toKernel());
        calcKernels.template deviceToHost();

        // free the resources
        cublasDestroy_v2(cublas_handle);
        if (d_infoArray) cudaFree(d_infoArray);
        if (h_infoArray) free(h_infoArray);
        if (d_A) cudaFree(d_A);
        if (d_b) cudaFree(d_b);
        if (d_A_pointers) cudaFree(d_A_pointers);
        if (d_b_pointers) cudaFree(d_b_pointers);
        if (h_A_pointers) free(h_A_pointers);
        if (h_b_pointers) free(h_b_pointers);
    }


    void initializeStaticSize(vector_type &particles,
                              unsigned int convergenceOrder,
                              T rCut,
                              T supportSizeFactor) {
        this->rCut=rCut;
        this->supportSizeFactor=supportSizeFactor;
        this->convergenceOrder=convergenceOrder;
        SupportBuilder<vector_type>
                supportBuilder(particles, differentialSignature, rCut);
        unsigned int requiredSupportSize = monomialBasis.size() * supportSizeFactor;

        localSupportKeys.resize(particles.size_local_orig());
        localSupportRefs.resize(particles.size_local_orig());
        localEps.resize(particles.size_local_orig());
        localEpsInvPow.resize(particles.size_local_orig());
        kerOffsets.resize(particles.size_local_orig());
        kerOffsets.fill(-1);

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        size_t offsetKer = 0;
        auto it = particles.getDomainIterator();
        while (it.isNext()) {
            // Get the points in the support of the DCPSE kernel and store the support for reuse
            Support support = supportBuilder.getSupport(it, requiredSupportSize,opt);

            auto key_o = particles.getOriginKey(it.get());
            localSupportRefs.get(key_o.getKey()) = support.getReferencePointKey();
            localSupportKeys.get(key_o.getKey()) = support.getKeys();
            kerOffsets.get(key_o.getKey()) = offsetKer;

            offsetKer += support.getKeys().size();
            ++it;
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << "Support building took " << time_span2.count() * 1000. << " milliseconds." << std::endl;

        std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

        // move monomial basis to kernel
        auto& basis = monomialBasis.getBasis();
        openfpm::vector_custd<Monomial_gpu<dim>> basisTemp(basis.begin(), basis.end());
        basisTemp.template hostToDevice();
        MonomialBasis<dim, aggregate<Monomial_gpu<dim>>, openfpm::vector_custd_ker, memory_traits_inte> monomialBasisKernel(basisTemp.toKernel());

        size_t numMatrices = localSupportKeys.size();
        size_t monomialBasisSize = monomialBasis.size();
        size_t localSupportRefKeysTotalN = 0;
        size_t EMatN = 0, VMatN = 0;

        for (size_t i = 0; i < numMatrices; ++i) {
            size_t localSupportSize = localSupportKeys.get(i).size();
            localSupportRefKeysTotalN += localSupportSize;
            EMatN += localSupportSize * localSupportSize;
            VMatN += localSupportSize * monomialBasisSize;
        }

        openfpm::vector_custd<size_t> localSupportKeys1D(localSupportRefKeysTotalN+numMatrices+1);
        openfpm::vector_custd<size_t> EMatOffsets(numMatrices+1);
        openfpm::vector_custd<size_t> VMatOffsets(numMatrices+1);
        openfpm::vector_custd<T> EMat(EMatN);
        openfpm::vector_custd<T> VMat(VMatN);
        // B has the same dimensions as V
        openfpm::vector_custd<T> BMat(VMatN);

        // localSupportKeys1D is populated:
        // offsets of keys first, followed by keys 
        size_t offset = numMatrices+1;
        for (size_t i = 0; i < numMatrices; ++i) {
            localSupportKeys1D.get(i) = offset;
            offset += localSupportKeys.get(i).size();
        }
        localSupportKeys1D.get(numMatrices) = offset;

        offset = numMatrices+1;
        for (size_t i = 0; i < numMatrices; ++i)
            for (size_t j = 0; j < localSupportKeys.get(i).size(); ++j, ++offset)
                localSupportKeys1D.get(offset) = localSupportKeys.get(i).get(j);

        // fill the matrices' offsets for E, V (B has the same dim as V)
        size_t offsetE = 0, offsetV = 0;
        for (size_t i = 0; i < numMatrices; ++i) {
            EMatOffsets.get(i) = offsetE;
            VMatOffsets.get(i) = offsetV;
            offsetE += localSupportKeys.get(i).size() * localSupportKeys.get(i).size();
            offsetV += localSupportKeys.get(i).size() * monomialBasisSize;
        }

        // allocate device space for A, b on device
        T *d_A; cudaMalloc((void**)&d_A, sizeof(T)*numMatrices*monomialBasisSize*monomialBasisSize);
        T *d_b; cudaMalloc((void**)&d_b, sizeof(T)*numMatrices*monomialBasisSize);

        // create array of pointers to pass T** pointers to cublas subroutines
        T **h_A_pointers = (T**)malloc(numMatrices*sizeof(T*));
        for (int i = 0; i < numMatrices; i++) h_A_pointers[i] = d_A + i*monomialBasisSize*monomialBasisSize;

        T **h_b_pointers = (T**)malloc(numMatrices*sizeof(T*));
        for (int i = 0; i < numMatrices; i++) h_b_pointers[i] = d_b + i*monomialBasisSize;

        T **d_A_pointers; cudaMalloc((void**)&d_A_pointers, numMatrices*sizeof(T*));
        cudaMemcpy(d_A_pointers, h_A_pointers, numMatrices*sizeof(T*), cudaMemcpyHostToDevice);

        T **d_b_pointers; cudaMalloc((void**)&d_b_pointers, numMatrices*sizeof(T*));
        cudaMemcpy(d_b_pointers, h_b_pointers, numMatrices*sizeof(T*), cudaMemcpyHostToDevice);

        // assemble local matrices on GPU
        particles.hostToDevicePos();
        localSupportKeys1D.template hostToDevice();
        EMatOffsets.template hostToDevice();
        VMatOffsets.template hostToDevice();
        localSupportRefs.template hostToDevice();

        auto it2 = particles.getDomainIteratorGPU(512);
        assembleLocalMatrices_gpu<<<it2.wthr,it2.thr>>>(particles.toKernel(), differentialSignature, differentialOrder, monomialBasisKernel, localSupportRefs.toKernel(), localSupportKeys1D.toKernel(), 
            d_A_pointers, d_b_pointers, localEps.toKernel(), localEpsInvPow.toKernel(), EMatOffsets.toKernel(), EMat.toKernel(), VMatOffsets.toKernel(), VMat.toKernel(), BMat.toKernel());

        localEps.template deviceToHost();
        localEpsInvPow.template deviceToHost();

        //cublas lu solver
        int *d_infoArray; cudaMalloc((void**)&d_infoArray,  numMatrices*sizeof(int));
        int *h_infoArray = (int *)malloc(numMatrices*sizeof(int));

        cublasHandle_t cublas_handle; cublasCreate_v2(&cublas_handle);

        // if (std::is_same<T, float>::value)
            // cublasSgetrfBatched(cublas_handle, monomialBasisSize, d_A_pointers, monomialBasisSize, NULL, d_infoArray, numMatrices);
        // else
            cublasDgetrfBatched(cublas_handle, monomialBasisSize, d_A_pointers, monomialBasisSize, NULL, d_infoArray, numMatrices);
        cudaDeviceSynchronize();

        cudaMemcpy(h_infoArray, d_infoArray, numMatrices*sizeof(int), cudaMemcpyDeviceToHost);
        for (size_t i = 0; i < numMatrices; i++)
            if (h_infoArray[i] != 0) fprintf(stderr, "Factorization of matrix %d Failed: Matrix may be singular\n", i);

        const double alpha = 1.f;
        // if (std::is_same<T, float>::value) {
        //     cublasStrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
        //     cublasStrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
        // } else {
            cublasDtrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
            cublasDtrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, monomialBasisSize, 1, &alpha, d_A_pointers, monomialBasisSize, d_b_pointers, monomialBasisSize, numMatrices);
        // }
        cudaDeviceSynchronize();

        // populate the calcKernels on GPU
        calcKernels.resize(localSupportRefKeysTotalN);
        localEps.template hostToDevice();
        calcKernels_gpu<dim><<<it2.wthr,it2.thr>>>(particles.toKernel(), monomialBasisKernel, localSupportKeys1D.toKernel(), d_b_pointers, localEps.toKernel(), numMatrices, calcKernels.toKernel());
        calcKernels.template deviceToHost();

        // free the resources
        cublasDestroy_v2(cublas_handle);
        if (d_infoArray) cudaFree(d_infoArray);
        if (h_infoArray) free(h_infoArray);
        if (d_A) cudaFree(d_A);
        if (d_b) cudaFree(d_b);
        if (d_A_pointers) cudaFree(d_A_pointers);
        if (d_b_pointers) cudaFree(d_b_pointers);
        if (h_A_pointers) free(h_A_pointers);
        if (h_b_pointers) free(h_b_pointers);

        std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
        std::cout << "Matrices inverse took " << time_span.count() * 1000. << " milliseconds." << std::endl;
    }


    T computeKernel(Point<dim, T> x, EMatrix<T, Eigen::Dynamic, 1> & a) const {
        T res = 0;
        unsigned int counter = 0;
        T expFactor = exp(-norm2(x));

        size_t N = monomialBasis.getElements().size();
        for (size_t i = 0; i < N; ++i)
        {
            const Monomial<dim> &m = monomialBasis.getElement(i);

            T coeff = a(counter);
            T mbValue = m.evaluate(x);
            res += coeff * mbValue * expFactor;
            ++counter;
        }
        return res;
    }


    // template <unsigned int a_dim>
    // T computeKernel(Point<dim, T> x, const T (& a) [a_dim]) const {
    T computeKernel(Point<dim, T> x, const T* a) const {
        T res = 0;
        unsigned int counter = 0;
        T expFactor = exp(-norm2(x));

        size_t N = monomialBasis.getElements().size();
        for (size_t i = 0; i < N; ++i)
        {
            const Monomial<dim> &m = monomialBasis.getElement(i);

            T coeff = a[counter];
            T mbValue = m.evaluate(x);
            res += coeff * mbValue * expFactor;
            ++counter;
        }
        return res;
    }


    T conditionNumber(const EMatrix<T, -1, -1> &V, T condTOL) const {
        std::cout << "conditionNumber" << std::endl;
        Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(V);
        T cond = svd.singularValues()(0)
                 / svd.singularValues()(svd.singularValues().size() - 1);
        if (cond > condTOL) {
            std::cout
                    << "WARNING: cond(V) = " << cond
                    << " is greater than TOL = " << condTOL
                    << ",  numPoints(V) = " << V.rows()
                    << std::endl; // debug
        }
        return cond;
    }

};


template<unsigned int dim, typename T, typename particles_type, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename matrixOffset_type, typename matrix_type>
__global__ void assembleLocalMatrices_gpu(
        particles_type particles, Point<dim, unsigned int> differentialSignature, unsigned int differentialOrder, monomialBasis_type monomialBasis, 
        supportKey_type localSupportRefs, supportKey_type localSupportKeys1D, T** h_A, T** h_b, localEps_type localEps, localEps_type localEpsInvPow, 
        matrixOffset_type EMatOffsets, matrix_type EMat, matrixOffset_type VMatOffsets, matrix_type VMat, matrix_type BMat) 
    {
    auto p_key = GET_PARTICLE(particles);
    Point<dim, T> xa = particles.getPos(p_key);

    size_t  monomialBasisSize = monomialBasis.size();
    size_t  supportKeysSize = localSupportKeys1D.get(p_key+1)-localSupportKeys1D.get(p_key);
    size_t* supportKeys = &((size_t*)localSupportKeys1D.getPointer())[localSupportKeys1D.get(p_key)];
    size_t  xpK = localSupportRefs.get(p_key);

    // EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
    // EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());
    size_t VStartPos = VMatOffsets.get(p_key);
    size_t EStartPos = EMatOffsets.get(p_key);
    size_t BStartPos = VMatOffsets.get(p_key);
    T* V = &((T*)VMat.getPointer())[VStartPos];
    T* E = &((T*)EMat.getPointer())[EStartPos];
    T* B = &((T*)BMat.getPointer())[BStartPos];

    // Vandermonde matrix computation
    // Pointer to E is passed to reuse memory for offset construction inside Vandermonde. 
    // Memory needed is supportKeysSize*dim. In rare case where supportKeysSize < dim, memory is allocated with new[] 
    Vandermonde_gpu<dim, T, MonomialBasis<dim, aggregate<Monomial_gpu<dim>>, openfpm::vector_custd_ker, memory_traits_inte>>
            vandermonde(E, xpK, supportKeysSize, supportKeys, monomialBasis, particles);
    vandermonde.getMatrix(V);

    T eps = vandermonde.getEps(); localEps.get(p_key) = eps;
    localEpsInvPow.get(p_key) = 1.0 / pow(eps,differentialOrder);

    DcpseDiagonalScalingMatrix<dim, MonomialBasis<dim, aggregate<Monomial_gpu<dim>>, openfpm::vector_custd_ker, memory_traits_inte>> diagonalScalingMatrix(monomialBasis);
    diagonalScalingMatrix.buildMatrix(E, xpK, supportKeysSize, supportKeys, eps, particles);

    // EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E * V;
    T sum = 0.0;
    for (int i = 0; i < supportKeysSize; ++i)
        for (int j = 0; j < monomialBasisSize; ++j) {
            for (int k = 0; k < supportKeysSize; ++k)
                sum += E[i*supportKeysSize+k] * V[k*monomialBasisSize+j];

            B[i*monomialBasisSize+j] = sum; sum = 0.0;
        }

    // EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = B.transpose() * B;
    for (int i = 0; i < monomialBasisSize; ++i)
        for (int j = 0; j < monomialBasisSize; ++j) {
            for (int k = 0; k < supportKeysSize; ++k)
                sum += B[k*monomialBasisSize+i] * B[k*monomialBasisSize+j];

            h_A[p_key][i*monomialBasisSize+j] = sum; sum = 0.0;
        }

    // Compute RHS vector b
    DcpseRhs_gpu<dim, MonomialBasis<dim, aggregate<Monomial_gpu<dim>>, openfpm::vector_custd_ker, memory_traits_inte>> rhs(monomialBasis, differentialSignature);
    rhs.template getVector<T>(h_b[p_key]);
}

template<unsigned int dim, typename particles_type, typename T, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename calcKernels_type>
__global__ void calcKernels_gpu(particles_type particles, monomialBasis_type monomialBasis, supportKey_type localSupportKeys1D, 
        T** h_b, localEps_type localEps, size_t numMatrices, calcKernels_type calcKernels)
    {
    auto p_key = GET_PARTICLE(particles);
    Point<dim, T> xa = particles.getPos(p_key);

    size_t  monomialBasisSize = monomialBasis.size();
    size_t  supportKeysSize = localSupportKeys1D.get(p_key+1)-localSupportKeys1D.get(p_key);
    size_t* supportKeys = &((size_t*)localSupportKeys1D.getPointer())[localSupportKeys1D.get(p_key)];
    T* d_calcKernels = &((T*)calcKernels.getPointer())[localSupportKeys1D.get(p_key)-numMatrices-1];
    T eps = localEps.get(p_key);

    for (size_t j = 0; j < supportKeysSize; ++j)
    {
        size_t xqK = supportKeys[j];
        Point<dim, T> xq = particles.getPos(xqK);
        Point<dim, T> normalizedArg = (xa - xq) / eps;

        d_calcKernels[j] = computeKernel_gpu(normalizedArg, h_b[p_key], monomialBasis);
    }
}

template <unsigned int dim, typename T, typename monomialBasis_type>
__device__ T computeKernel_gpu(Point<dim, T>& x, const T* a, const monomialBasis_type& monomialBasis) {
    T res = 0;
    unsigned int counter = 0;
    T expFactor = exp(-norm2(x));

    const auto& basisElements = monomialBasis.getElements();

    size_t N = basisElements.size();
    for (size_t i = 0; i < N; ++i)
    {
        const Monomial_gpu<dim> &m = basisElements.get(i);

        T coeff = a[counter];
        T mbValue = m.evaluate(x);
        res += coeff * mbValue * expFactor;
        ++counter;
    }
    return res;
}


#endif
#endif //OPENFPM_PDATA_DCPSE_CUH

