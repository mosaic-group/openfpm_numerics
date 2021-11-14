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
#include "SupportBuilder.cuh"
#include "Support.hpp"
#include "Vandermonde.hpp"
#include "DcpseDiagonalScalingMatrix.hpp"
#include "DcpseRhs.hpp"

#include <chrono>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>


template<unsigned int dim, typename particles_type, typename T, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename calcKernels_type>
__global__ void calcKernels_gpu(particles_type, monomialBasis_type, supportKey_type, supportKey_type, T**, localEps_type, size_t, calcKernels_type);

template<unsigned int dim, typename T, typename particles_type, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename matrix_type>
__global__ void assembleLocalMatrices_gpu( particles_type, Point<dim, unsigned int>, unsigned int, monomialBasis_type, supportKey_type, supportKey_type, supportKey_type,
    T**, T**, localEps_type, localEps_type, matrix_type, size_t, size_t);


template<unsigned int dim, typename vector_type, class T = typename vector_type::stype>
class Dcpse_gpu {
    static_assert(std::is_floating_point<T>::value, "CUBLAS supports only float or double");

public:
    typedef typename vector_type::value_type part_type;
    typedef vector_type vtype;

    #ifdef SE_CLASS1
    int update_ctr=0;
    #endif
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

    // shared local support previosly built by another operator
    bool isSharedSupport = false;
    openfpm::vector_custd<size_t> supportRefs; // Each MPI rank has just access to the local ones
    openfpm::vector_custd<size_t> kerOffsets;
    openfpm::vector_custd<size_t> supportKeys1D;

    openfpm::vector_custd<T> localEps; // Each MPI rank has just access to the local ones
    openfpm::vector_custd<T> localEpsInvPow; // Each MPI rank has just access to the local ones
    openfpm::vector_custd<T> calcKernels;

    openfpm::vector<size_t> subsetKeyPid;

    vector_type & particles;
    double rCut;
    unsigned int convergenceOrder;
    double supportSizeFactor;

    size_t maxSupportSize;
    size_t supportKeysTotalN;

    support_options opt;

public:
#ifdef SE_CLASS1
    int getUpdateCtr() const
    {
        return update_ctr;
    }
#endif

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
            maxSupportSize(0),
            supportKeysTotalN(0),
            opt(opt)
    {
        particles.ghost_get_subset();

        if (supportSizeFactor < 1) 
            initializeAdaptive(particles, convergenceOrder, rCut);
        else 
            initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
    }

    Dcpse_gpu(vector_type &particles,
          const Dcpse_gpu<dim, vector_type, T>& other,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1,
          support_options opt = support_options::N_PARTICLES)
        :particles(particles), opt(opt),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            subsetKeyPid(other.subsetKeyPid),
            supportRefs(other.supportRefs),
            supportKeys1D(other.supportKeys1D),
            kerOffsets(other.kerOffsets),
            maxSupportSize(other.maxSupportSize),
            supportKeysTotalN(other.supportKeysTotalN),
            isSharedSupport(true)
    {
        particles.ghost_get_subset();

        if (supportSizeFactor < 1)
            initializeAdaptive(particles, convergenceOrder, rCut);
        else
            initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k)
    {
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);

        size_t  supportKeysSize = kerOffsets.get(k+1)-kerOffsets.get(k);
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(k)];

        for (int i = 0; i < supportKeysSize; i++)
        {
            size_t xqK = supportKeys[i];

            particles.template getProp<prp>(xqK) += calcKernels.get(kerOff+i);
        }
    }

    template<unsigned int prp>
    void DrawKernelNN(vector_type &particles, int k)
    {
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        size_t  supportKeysSize = kerOffsets.get(k+1)-kerOffsets.get(k);
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(k)];

        for (int i = 0; i < supportKeysSize; i++)
        {
            size_t xqK = supportKeys[i];

            particles.template getProp<prp>(xqK) = 1.0;
        }
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k, int i)
    {
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        size_t  supportKeysSize = kerOffsets.get(k+1)-kerOffsets.get(k);
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(k)];

        for (int i = 0; i < supportKeysSize; i++)
        {
            size_t xqK = supportKeys[i];

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

        size_t N = particles.size_local();
        for (size_t j = 0; j < N; ++j)
        {
            double eps = localEps.get(j);

            for (int i = 0; i < momenta.size(); i++)
            {
                momenta_accu.template get<0>(i) =  0.0;
            }

            size_t xpK = supportRefs.get(j);
            Point<dim, T> xp = particles.getPos(xpK);

            size_t kerOff = kerOffsets.get(xpK);
            size_t  supportKeysSize = kerOffsets.get(j+1)-kerOffsets.get(j);
            size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(j)];

            for (int i = 0; i < supportKeysSize; i++)
            {
                size_t xqK = supportKeys[i];
                Point<dim, T> xq = particles.getPos(xqK);
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

        size_t N = particles.size_local();
        for (size_t j = 0; j < N; ++j)
        {
            double epsInvPow = localEpsInvPow.get(j);

            T Dfxp = 0;
            size_t xpK = supportRefs.get(j);
            Point<dim, typename vector_type::stype> xp = particles.getPos(xpK);
            T fxp = sign * particles.template getProp<fValuePos>(xpK);
            size_t kerOff = kerOffsets.get(xpK);

            size_t  supportKeysSize = kerOffsets.get(j+1)-kerOffsets.get(j);
            size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(j)];

            for (int i = 0; i < supportKeysSize; i++)
            {
                size_t xqK = supportKeys[i];
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
        return kerOffsets.get(key.getKey()+1)-kerOffsets.get(key.getKey());
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
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(key.getKey())];
        return supportKeys[j];
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

        size_t localKey = subsetKeyPid.get(key.getKey());

        double eps = localEps.get(localKey);
        double epsInvPow = localEpsInvPow.get(localKey);

        auto &particles = o1.getVector();

#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif


        expr_type Dfxp = 0;
        size_t xpK = supportRefs.get(localKey);
        Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key);
        size_t kerOff = kerOffsets.get(xpK);

        size_t  supportKeysSize = kerOffsets.get(localKey+1)-kerOffsets.get(localKey);
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(localKey)];

        for (int i = 0; i < supportKeysSize; i++)
        {
            size_t xqK = supportKeys[i];
            expr_type fxq = o1.value(vect_dist_key_dx(xqK));
            Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+i);
        }
        Dfxp = Dfxp * epsInvPow;

        // T trueDfxp = particles.template getProp<2>(xpK);
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

        size_t localKey = subsetKeyPid.get(key.getKey());

        double eps = localEps.get(localKey);
        double epsInvPow = localEpsInvPow(localKey);

        auto &particles = o1.getVector();

#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif

        expr_type Dfxp = 0;
        size_t xpK = supportRefs.get(localKey);

        Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key)[i];
        size_t kerOff = kerOffsets.get(xpK);
        size_t  supportKeysSize = kerOffsets.get(localKey+1)-kerOffsets.get(localKey);
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(localKey)];

        for (int j = 0; j < supportKeysSize; j++)
        {
            size_t xqK = supportKeys[j];
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
#ifdef SE_CLASS1
        update_ctr=particles.getMapCtr();
#endif

        kerOffsets.clear();
        supportKeys1D.clear();
        supportRefs.clear();
        localEps.clear();
        localEpsInvPow.clear();
        calcKernels.clear();
        subsetKeyPid.clear();

        initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
    }

private:

    void initializeAdaptive(vector_type &particles,
                            unsigned int convergenceOrder,
                            double rCut) {
        // Still need to be tested
#ifdef SE_CLASS1
        this->update_ctr=particles.getMapCtr();
#endif
        SupportBuilder<vector_type> supportBuilder(particles, differentialSignature, rCut);
        unsigned int requiredSupportSize = monomialBasis.size();

        subsetKeyPid.resize(particles.size_local_orig());
        if (!isSharedSupport)
            supportRefs.resize(particles.size_local());
        localEps.resize(particles.size_local());
        localEpsInvPow.resize(particles.size_local());
        kerOffsets.resize(particles.size_local()+1);

        // need to resize supportKeys1D to yet unknown supportKeysTotalN
        // add() takes too long
        openfpm::vector<openfpm::vector<size_t>> tempSupportKeys(supportRefs.size());
        const T condVTOL = 1e2;

        auto it = particles.getDomainIterator();
        while (it.isNext()) {
            size_t localSupportSize;
            auto key_o = it.get(); subsetKeyPid.get(particles.getOriginKey(key_o).getKey()) = key_o.getKey();

            if (!isSharedSupport){
                auto key_o = particles.getOriginKey(it.get());

                Support support = supportBuilder.getSupport(it, requiredSupportSize, opt);
                supportRefs.get(key_o.getKey()) = key_o.getKey();
                tempSupportKeys.get(key_o.getKey()) = support.getKeys();
                kerOffsets.get(key_o.getKey()) = supportKeysTotalN;

                if (maxSupportSize < support.size()) maxSupportSize = support.size();
                supportKeysTotalN += support.size();


                EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());
                // Vandermonde matrix computation
                Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                        vandermonde(support, monomialBasis, particles);
                vandermonde.getMatrix(V);

                T condV = conditionNumber(V, condVTOL);
                T eps = vandermonde.getEps();

                if (condV > condVTOL) {
                    requiredSupportSize *= 2;
                    std::cout << "INFO: Increasing, requiredSupportSize = " << requiredSupportSize << std::endl; // debug
                    continue;
                } else requiredSupportSize = monomialBasis.size();
            }
            ++it;
        }

        if (!isSharedSupport){
            kerOffsets.get(supportRefs.size()) = supportKeysTotalN;
            supportKeys1D.resize(supportKeysTotalN);

            size_t offset = 0;
            for (size_t i = 0; i < tempSupportKeys.size(); ++i)
                for (size_t j = 0; j < tempSupportKeys.get(i).size(); ++j, ++offset)
                    supportKeys1D.get(offset) = tempSupportKeys.get(i).get(j);
        }

        kerOffsets.hostToDevice(); supportKeys1D.hostToDevice();
        assembleLocalMatrices(cublasDgetrfBatched, cublasDtrsmBatched);
    }

    void initializeAdaptive(vector_type &particles,
                            unsigned int convergenceOrder,
                            float rCut) {
        // Still need to be tested
#ifdef SE_CLASS1
        this->update_ctr=particles.getMapCtr();
#endif
        SupportBuilder<vector_type> supportBuilder(particles, differentialSignature, rCut);
        unsigned int requiredSupportSize = monomialBasis.size();

        subsetKeyPid.resize(particles.size_local_orig());
        if (!isSharedSupport)
            supportRefs.resize(particles.size_local());
        localEps.resize(particles.size_local());
        localEpsInvPow.resize(particles.size_local());
        kerOffsets.resize(particles.size_local()+1);

        // need to resize supportKeys1D to yet unknown supportKeysTotalN
        // add() takes too long
        openfpm::vector<openfpm::vector<size_t>> tempSupportKeys(supportRefs.size());
        const T condVTOL = 1e2;

        auto it = particles.getDomainIterator();
        while (it.isNext()) {
            size_t localSupportSize;
            auto key_o = it.get(); subsetKeyPid.get(particles.getOriginKey(key_o).getKey()) = key_o.getKey();

            if (!isSharedSupport){
                auto key_o = particles.getOriginKey(it.get());

                Support support = supportBuilder.getSupport(it, requiredSupportSize, opt);
                supportRefs.get(key_o.getKey()) = key_o.getKey();
                tempSupportKeys.get(key_o.getKey()) = support.getKeys();
                kerOffsets.get(key_o.getKey()) = supportKeysTotalN;

                if (maxSupportSize < support.size()) maxSupportSize = support.size();
                supportKeysTotalN += support.size();


                EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());
                // Vandermonde matrix computation
                Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                        vandermonde(support, monomialBasis, particles);
                vandermonde.getMatrix(V);

                T condV = conditionNumber(V, condVTOL);
                T eps = vandermonde.getEps();

                if (condV > condVTOL) {
                    requiredSupportSize *= 2;
                    std::cout << "INFO: Increasing, requiredSupportSize = " << requiredSupportSize << std::endl; // debug
                    continue;
                } else requiredSupportSize = monomialBasis.size();
            }
            ++it;
        }

        if (!isSharedSupport){
            kerOffsets.get(supportRefs.size()) = supportKeysTotalN;
            supportKeys1D.resize(supportKeysTotalN);

            size_t offset = 0;
            for (size_t i = 0; i < tempSupportKeys.size(); ++i)
                for (size_t j = 0; j < tempSupportKeys.get(i).size(); ++j, ++offset)
                    supportKeys1D.get(offset) = tempSupportKeys.get(i).get(j);
        }

        kerOffsets.hostToDevice(); supportKeys1D.hostToDevice();
        assembleLocalMatrices(cublasSgetrfBatched, cublasStrsmBatched);
    }

    void initializeStaticSize(vector_type &particles,
                              unsigned int convergenceOrder,
                              double rCut,
                              double supportSizeFactor) {
#ifdef SE_CLASS1
        this->update_ctr=particles.getMapCtr();
#endif
        this->rCut=rCut;
        this->supportSizeFactor=supportSizeFactor;
        this->convergenceOrder=convergenceOrder;

        subsetKeyPid.resize(particles.size_local_orig());
        if (!isSharedSupport)
            supportRefs.resize(particles.size_local());
        localEps.resize(particles.size_local());
        localEpsInvPow.resize(particles.size_local());

std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        auto it = particles.getDomainIterator();

        if (opt==support_options::RADIUS) {
            while (it.isNext()) {
                auto key_o = it.get(); subsetKeyPid.get(particles.getOriginKey(key_o).getKey()) = key_o.getKey();

                if (!isSharedSupport)
                    supportRefs.get(key_o.getKey()) = key_o.getKey();
                ++it;
            }

            SupportBuilderGPU<vector_type> supportBuilder(particles, differentialSignature, rCut);
            supportBuilder.getSupport(supportRefs.size(), kerOffsets, supportKeys1D, maxSupportSize, supportKeysTotalN);
        } else {
            size_t requiredSupportSize = monomialBasis.size() * supportSizeFactor;
            SupportBuilder<vector_type> supportBuilder(particles, differentialSignature, rCut);

            kerOffsets.resize(supportRefs.size()+1);
            // need to resize supportKeys1D to yet unknown supportKeysTotalN
            // add() takes too long
            openfpm::vector<openfpm::vector<size_t>> tempSupportKeys(supportRefs.size());

            while (it.isNext()) {
                if (!isSharedSupport){
                    auto key_o = it.get(); subsetKeyPid.get(particles.getOriginKey(key_o).getKey()) = key_o.getKey();

                    Support support = supportBuilder.getSupport(it, requiredSupportSize, opt);
                    supportRefs.get(key_o.getKey()) = key_o.getKey();
                    tempSupportKeys.get(key_o.getKey()) = support.getKeys();
                    kerOffsets.get(key_o.getKey()) = supportKeysTotalN;

                    if (maxSupportSize < support.size()) maxSupportSize = support.size();
                    supportKeysTotalN += support.size();
                }
                ++it;
            }

            if (!isSharedSupport){
                kerOffsets.get(supportRefs.size()) = supportKeysTotalN;
                supportKeys1D.resize(supportKeysTotalN);

                size_t offset = 0;
                for (size_t i = 0; i < tempSupportKeys.size(); ++i)
                    for (size_t j = 0; j < tempSupportKeys.get(i).size(); ++j, ++offset)
                        supportKeys1D.get(offset) = tempSupportKeys.get(i).get(j);
            }

            kerOffsets.hostToDevice(); supportKeys1D.hostToDevice();
        }

std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
std::cout << "Support building took " << time_span2.count() * 1000. << " milliseconds." << std::endl;

        assembleLocalMatrices(cublasDgetrfBatched, cublasDtrsmBatched);
    }

    // ad hoc solution to template specialization for float/double
    void initializeStaticSize(vector_type &particles,
                              unsigned int convergenceOrder,
                              float rCut,
                              float supportSizeFactor) {
#ifdef SE_CLASS1
        this->update_ctr=particles.getMapCtr();
#endif
        this->rCut=rCut;
        this->supportSizeFactor=supportSizeFactor;
        this->convergenceOrder=convergenceOrder;

        subsetKeyPid.resize(particles.size_local_orig());
        if (!isSharedSupport)
            supportRefs.resize(particles.size_local());
        localEps.resize(particles.size_local());
        localEpsInvPow.resize(particles.size_local());

std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        auto it = particles.getDomainIterator();

        if (opt==support_options::RADIUS) {
            while (it.isNext()) {
                auto key_o = it.get(); subsetKeyPid.get(particles.getOriginKey(key_o).getKey()) = key_o.getKey();

                if (!isSharedSupport)
                    supportRefs.get(key_o.getKey()) = key_o.getKey();
                ++it;
            }

            SupportBuilderGPU<vector_type> supportBuilder(particles, differentialSignature, rCut);
            supportBuilder.getSupport(supportRefs.size(), kerOffsets, supportKeys1D, maxSupportSize, supportKeysTotalN);

        } else {
            size_t requiredSupportSize = monomialBasis.size() * supportSizeFactor;
            SupportBuilder<vector_type> supportBuilder(particles, differentialSignature, rCut);

            kerOffsets.resize(supportRefs.size()+1);
            // need to resize supportKeys1D to yet unknown supportKeysTotalN
            // add() takes too long
            openfpm::vector<openfpm::vector<size_t>> tempSupportKeys(supportRefs.size());

            while (it.isNext()) {
                if (!isSharedSupport){
                    auto key_o = it.get(); subsetKeyPid.get(particles.getOriginKey(key_o).getKey()) = key_o.getKey();

                    Support support = supportBuilder.getSupport(it, requiredSupportSize, opt);
                    supportRefs.get(key_o.getKey()) = key_o.getKey();
                    tempSupportKeys.get(key_o.getKey()) = support.getKeys();
                    kerOffsets.get(key_o.getKey()) = supportKeysTotalN;

                    if (maxSupportSize < support.size()) maxSupportSize = support.size();
                    supportKeysTotalN += support.size();
                }
                ++it;
            }

            if (!isSharedSupport){
                kerOffsets.get(supportRefs.size()) = supportKeysTotalN;
                supportKeys1D.resize(supportKeysTotalN);

                size_t offset = 0;
                for (size_t i = 0; i < tempSupportKeys.size(); ++i)
                    for (size_t j = 0; j < tempSupportKeys.get(i).size(); ++j, ++offset)
                        supportKeys1D.get(offset) = tempSupportKeys.get(i).get(j);
            }

            kerOffsets.hostToDevice(); supportKeys1D.hostToDevice();
        }

std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
std::cout << "Support building took " << time_span2.count() * 1000. << " milliseconds." << std::endl;

        assembleLocalMatrices(cublasSgetrfBatched, cublasStrsmBatched);
    }

    template<typename cublasLUDec_type, typename cublasTriangSolve_type>
    void assembleLocalMatrices(cublasLUDec_type cublasLUDecFunc, cublasTriangSolve_type cublasTriangSolveFunc) {
        std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();

        // move monomial basis to kernel
        auto& basis = monomialBasis.getBasis();
        openfpm::vector_custd<Monomial_gpu<dim>> basisTemp(basis.begin(), basis.end());
        basisTemp.template hostToDevice();
        MonomialBasis<dim, aggregate<Monomial_gpu<dim>>, openfpm::vector_custd_ker, memory_traits_inte> monomialBasisKernel(basisTemp.toKernel());

        size_t numMatrices = supportRefs.size();
        size_t monomialBasisSize = monomialBasis.size();

        int numSMs, numSMsMult = 1;
        cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
        size_t numThreads = numSMs*numSMsMult*256;
        std::cout << "numThreads " << numThreads << " numMatrices " << numMatrices << std::endl;

        // B is an intermediate matrix
        openfpm::vector_custd<T> BMat(numThreads * maxSupportSize * monomialBasisSize);
        // allocate device space for A, b
        openfpm::vector_custd<T> AMat(numMatrices*monomialBasisSize*monomialBasisSize);
        openfpm::vector_custd<T> bVec(numMatrices*monomialBasisSize);

        // create array of pointers to pass T** pointers to cublas subroutines
        openfpm::vector_custd<T*> AMatPointers(numMatrices);
        openfpm::vector_custd<T*> bVecPointers(numMatrices);

        auto AMatKernel = AMat.toKernel(); T* AMatKernelPointer = (T*) AMatKernel.getPointer();
        for (size_t i = 0; i < numMatrices; i++) AMatPointers.get(i) = AMatKernelPointer + i*monomialBasisSize*monomialBasisSize;

        auto bVecKernel = bVec.toKernel(); T* bVecKernelPointer = (T*) bVecKernel.getPointer();
        for (size_t i = 0; i < numMatrices; i++) bVecPointers.get(i) = bVecKernelPointer + i*monomialBasisSize;

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span0 = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t3);
        std::cout << "Preallocation took " << time_span0.count() * 1000. << " milliseconds." << std::endl;

        // assemble local matrices on GPU
        std::chrono::high_resolution_clock::time_point t9 = std::chrono::high_resolution_clock::now();
        particles.hostToDevicePos();
        supportRefs.template hostToDevice();
        AMatPointers.template hostToDevice();
        bVecPointers.template hostToDevice();

        auto AMatPointersKernel = AMatPointers.toKernel(); T** AMatPointersKernelPointer = (T**) AMatPointersKernel.getPointer();
        auto bVecPointersKernel = bVecPointers.toKernel(); T** bVecPointersKernelPointer = (T**) bVecPointersKernel.getPointer();

        assembleLocalMatrices_gpu<<<numSMsMult*numSMs, 256>>>(particles.toKernel(), differentialSignature, differentialOrder, monomialBasisKernel, supportRefs.toKernel(), kerOffsets.toKernel(), supportKeys1D.toKernel(),
            AMatPointersKernelPointer, bVecPointersKernelPointer, localEps.toKernel(), localEpsInvPow.toKernel(), BMat.toKernel(), numMatrices, maxSupportSize);

        localEps.template deviceToHost();
        localEpsInvPow.template deviceToHost();

        std::chrono::high_resolution_clock::time_point t10 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span3 = std::chrono::duration_cast<std::chrono::duration<double>>(t10 - t9);
        std::cout << "assembleLocalMatrices_gpu took " << time_span3.count() * 1000. << " milliseconds." << std::endl;

        //cublas lu solver
        std::chrono::high_resolution_clock::time_point t7 = std::chrono::high_resolution_clock::now();
        cublasHandle_t cublas_handle; cublasCreate_v2(&cublas_handle);

        openfpm::vector_custd<int> infoArray(numMatrices); auto infoArrayKernel = infoArray.toKernel();
        cublasLUDecFunc(cublas_handle, monomialBasisSize, AMatPointersKernelPointer, monomialBasisSize, NULL, (int*) infoArrayKernel.getPointer(), numMatrices);
        cudaDeviceSynchronize();

        infoArray.template deviceToHost();
        for (size_t i = 0; i < numMatrices; i++)
            if (infoArray.get(i) != 0) fprintf(stderr, "Factorization of matrix %d Failed: Matrix may be singular\n", i);

        const double alpha = 1.f;
        cublasTriangSolveFunc(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, monomialBasisSize, 1, &alpha, AMatPointersKernelPointer, monomialBasisSize, bVecPointersKernelPointer, monomialBasisSize, numMatrices);
        cublasTriangSolveFunc(cublas_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, monomialBasisSize, 1, &alpha, AMatPointersKernelPointer, monomialBasisSize, bVecPointersKernelPointer, monomialBasisSize, numMatrices);
        cudaDeviceSynchronize();

        std::chrono::high_resolution_clock::time_point t8 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span4 = std::chrono::duration_cast<std::chrono::duration<double>>(t8 - t7);
        std::cout << "cublas took " << time_span4.count() * 1000. << " milliseconds." << std::endl;

        std::chrono::high_resolution_clock::time_point t5 = std::chrono::high_resolution_clock::now();
        // populate the calcKernels on GPU
        calcKernels.resize(supportKeysTotalN);
        localEps.template hostToDevice();
        auto it2 = particles.getDomainIteratorGPU(512);
        calcKernels_gpu<dim><<<it2.wthr,it2.thr>>>(particles.toKernel(), monomialBasisKernel, kerOffsets.toKernel(), supportKeys1D.toKernel(), bVecPointersKernelPointer, localEps.toKernel(), numMatrices, calcKernels.toKernel());
        calcKernels.template deviceToHost();

        std::chrono::high_resolution_clock::time_point t6 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span5 = std::chrono::duration_cast<std::chrono::duration<double>>(t6 - t5);
        std::cout << "calcKernels_gpu took " << time_span5.count() * 1000. << " milliseconds." << std::endl;

        // free the resources
        cublasDestroy_v2(cublas_handle);

        std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3);
        std::cout << "Matrices inverse took " << time_span.count() * 1000. << " milliseconds." << std::endl;
    }

    T computeKernel(Point<dim, T> x, EMatrix<T, Eigen::Dynamic, 1> & a) const {
        unsigned int counter = 0;
        T res = 0, expFactor = exp(-norm2(x));

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
        unsigned int counter = 0;
        T res = 0, expFactor = exp(-norm2(x));

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


template<unsigned int dim, typename T, typename particles_type, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename matrix_type>
__global__ void assembleLocalMatrices_gpu(
        particles_type particles, Point<dim, unsigned int> differentialSignature, unsigned int differentialOrder, monomialBasis_type monomialBasis, 
        supportKey_type supportRefs, supportKey_type kerOffsets, supportKey_type supportKeys1D, T** h_A, T** h_b, localEps_type localEps, localEps_type localEpsInvPow,
        matrix_type BMat, size_t numMatrices, size_t maxSupportSize)
    {
    auto p_key = GET_PARTICLE(particles);
    size_t monomialBasisSize = monomialBasis.size();
    size_t BStartPos = maxSupportSize * monomialBasisSize * p_key; T* B = &((T*)BMat.getPointer())[BStartPos];
    const auto& basisElements = monomialBasis.getElements();
    int rhsSign = (Monomial_gpu<dim>(differentialSignature).order() % 2 == 0) ? 1 : -1;

    for (; 
        p_key < numMatrices; 
        p_key += blockDim.x * gridDim.x) 
    {
        Point<dim, T> xa = particles.getPos(p_key);

        size_t  supportKeysSize = kerOffsets.get(p_key+1)-kerOffsets.get(p_key);
        size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(p_key)];
        size_t  xpK = supportRefs.get(p_key);

    assert(supportKeysSize >= monomialBasis.size());

        T FACTOR = 2, avgNeighbourSpacing = 0;

        for (int i = 0 ; i < supportKeysSize; i++) {
            Point<dim,T> off = xa; off -= particles.getPosOrig(supportKeys[i]);

            for (size_t j = 0; j < dim; ++j)
                avgNeighbourSpacing += fabs(off.value(j));
        }

        avgNeighbourSpacing /= supportKeysSize;
        T eps = FACTOR * avgNeighbourSpacing;

    assert(eps != 0);

        localEps.get(p_key) = eps;
        localEpsInvPow.get(p_key) = 1.0 / pow(eps,differentialOrder);

        // EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E * V;
        for (int i = 0; i < supportKeysSize; ++i)
            for (int j = 0; j < monomialBasisSize; ++j) {
                Point<dim,T> off = xa; off -= particles.getPosOrig(supportKeys[i]);
                const Monomial_gpu<dim>& m = basisElements.get(j);

                T V_ij = m.evaluate(off) / pow(eps, m.order());
                T E_ii = exp(- norm2(off) / (2.0 * eps * eps));
                B[i*monomialBasisSize+j] = E_ii * V_ij;
            }

        T sum = 0.0;
        // EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = B.transpose() * B;
        for (int i = 0; i < monomialBasisSize; ++i)
            for (int j = 0; j < monomialBasisSize; ++j) {
                for (int k = 0; k < supportKeysSize; ++k)
                    sum += B[k*monomialBasisSize+i] * B[k*monomialBasisSize+j];

                h_A[p_key][i*monomialBasisSize+j] = sum; sum = 0.0;
            }

        // Compute RHS vector b
        for (size_t i = 0; i < monomialBasisSize; ++i) {
            const Monomial_gpu<dim>& dm = basisElements.get(i).getDerivative(differentialSignature);
            h_b[p_key][i] = rhsSign * dm.evaluate(Point<dim, T>(0));
        }
    }
}

template<unsigned int dim, typename particles_type, typename T, typename monomialBasis_type, typename supportKey_type, typename localEps_type, typename calcKernels_type>
__global__ void calcKernels_gpu(particles_type particles, monomialBasis_type monomialBasis, supportKey_type kerOffsets, supportKey_type supportKeys1D,
        T** h_b, localEps_type localEps, size_t numMatrices, calcKernels_type calcKernels)
    {
    auto p_key = GET_PARTICLE(particles);
    Point<dim, T> xa = particles.getPos(p_key);

    size_t  monomialBasisSize = monomialBasis.size();
    const auto& basisElements = monomialBasis.getElements();
    size_t  supportKeysSize = kerOffsets.get(p_key+1)-kerOffsets.get(p_key);
    size_t* supportKeys = &((size_t*)supportKeys1D.getPointer())[kerOffsets.get(p_key)];

    T* calcKernelsLocal = &((T*)calcKernels.getPointer())[kerOffsets.get(p_key)];
    T eps = localEps.get(p_key);

    for (size_t j = 0; j < supportKeysSize; ++j)
    {
        size_t xqK = supportKeys[j];
        Point<dim, T> xq = particles.getPosOrig(xqK);
        Point<dim, T> offNorm = (xa - xq) / eps;
        T expFactor = exp(-norm2(offNorm));

        T res = 0;
        for (size_t i = 0; i < monomialBasisSize; ++i) {
            const Monomial_gpu<dim> &m = basisElements.get(i);
            T mbValue = m.evaluate(offNorm);
            T coeff = h_b[p_key][i];

            res += coeff * mbValue * expFactor;
        }
        calcKernelsLocal[j] = res;
    }
}

#endif
#endif //OPENFPM_PDATA_DCPSE_CUH

