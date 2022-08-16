//
// DCPSE Created by tommaso on 29/03/19.
// Modified, Updated and Maintained by Abhinav and Pietro
//Surface Operators by Abhinav Singh on 07/10/2021
#ifndef OPENFPM_PDATA_DCPSE_HPP
#define OPENFPM_PDATA_DCPSE_HPP

#ifdef HAVE_EIGEN

#include "Vector/vector_dist.hpp"
#include "MonomialBasis.hpp"
#include "DMatrix/EMatrix.hpp"
#include "SupportBuilder.hpp"
#include "Support.hpp"
#include "Vandermonde.hpp"
#include "DcpseDiagonalScalingMatrix.hpp"
#include "DcpseRhs.hpp"
#include "hash_map/hopscotch_map.h"

template<unsigned int N> struct value_t {};

template<bool cond>
struct is_scalar {
	template<typename op_type>
	static auto
	analyze(const vect_dist_key_dx &key, op_type &o1) -> typename std::remove_reference<decltype(o1.value(
			key))>::type {
		return o1.value(key);
	};
};

template<>
struct is_scalar<false> {
	template<typename op_type>
	static auto
	analyze(const vect_dist_key_dx &key, op_type &o1) -> typename std::remove_reference<decltype(o1.value(
			key))>::type {
		return o1.value(key);
	};
};

template<unsigned int dim, typename vector_type,typename vector_type2=vector_type>
class Dcpse {
public:

    typedef typename vector_type::stype T;
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
    ////c=HOverEpsilon. Note that the Eps value is computed by <h>/c (<h>=local average spacing for each particle and its support). This factor c is used in the Vandermonde.hpp.
    double HOverEpsilon=0.9;
private:
    const Point<dim, unsigned int> differentialSignature;
    const unsigned int differentialOrder;
    const MonomialBasis<dim> monomialBasis;

    bool isSharedLocalSupport = false;
    openfpm::vector<Support> localSupports; // Each MPI rank has just access to the local ones
    openfpm::vector<T> localEps; // Each MPI rank has just access to the local ones
    openfpm::vector<T> localEpsInvPow; // Each MPI rank has just access to the local ones

    openfpm::vector<size_t> kerOffsets,accKerOffsets;
    openfpm::vector<T> calcKernels;
    openfpm::vector<T> accCalcKernels;
    openfpm::vector<T> nSpacings;
    vector_type & particlesFrom;
    vector_type2 & particlesTo;
    double rCut,supportSizeFactor=1,nSpacing;
    unsigned int convergenceOrder,nCount;

    bool isSurfaceDerivative=false;
    size_t initialParticleSize;


    support_options opt;
public:
    template<unsigned int NORMAL_ID>
    void createNormalParticles(vector_type &particles)
    {
        particles.template ghost_get<NORMAL_ID>(SKIP_LABELLING);
        initialParticleSize=particles.size_local_with_ghost();
        auto it = particles.getDomainAndGhostIterator();
        while(it.isNext()){
            auto key=it.get();
            Point<dim,T> xp=particles.getPos(key), Normals=particles.template getProp<NORMAL_ID>(key);
            if(opt==support_options::ADAPTIVE_SURFACE)
            {
                nSpacing=nSpacings.get(key.getKey());
            }
            for(int i=1;i<=nCount;i++){
                particles.addAtEnd();
                for(size_t j=0;j<dim;j++)
                {particles.getLastPosEnd()[j]=xp[j]+i*nSpacing*Normals[j];}
                particles.addAtEnd();
                for(size_t j=0;j<dim;j++)
                {particles.getLastPosEnd()[j]=xp[j]-i*nSpacing*Normals[j];}
            }
            ++it;
        }
    }

    void accumulateAndDeleteNormalParticles(vector_type &particles)
    {
        tsl::hopscotch_map<size_t, size_t> nMap;
        auto it = particles.getDomainIterator();
        auto supportsIt = localSupports.begin();
        openfpm::vector_std<size_t> supportBuffer;
        accCalcKernels.clear();
        accKerOffsets.clear();
        accKerOffsets.resize(initialParticleSize);
        accKerOffsets.fill(-1);
        while(it.isNext()){
            supportBuffer.clear();
            nMap.clear();
            auto key=it.get();
            Support support = *supportsIt;
            size_t xpK = support.getReferencePointKey();
            size_t kerOff = kerOffsets.get(xpK);
            auto &keys = support.getKeys();
            accKerOffsets.get(xpK)=accCalcKernels.size();
            for (int i = 0 ; i < keys.size() ; i++)
            {
                size_t xqK = keys.get(i);
                int real_particle=(xqK-initialParticleSize)/(2.*nCount);
                if(real_particle<0)
                    {
                        real_particle=xqK;
                    }
                auto found=nMap.find(real_particle);
                if(found!=nMap.end()){
                    accCalcKernels.get(found->second)+=calcKernels.get(kerOff+i);
                }
                else{
                    supportBuffer.add();
                    supportBuffer.get(supportBuffer.size()-1)=real_particle;
                    accCalcKernels.add();
                    accCalcKernels.get(accCalcKernels.size()-1)=calcKernels.get(kerOff+i);
                    nMap[real_particle]=accCalcKernels.size()-1;
                }
            }
            keys.swap(supportBuffer);
            localSupports.get(xpK) = support;
            ++supportsIt;
            ++it;
        }
        particles.resizeAtEnd(initialParticleSize);
        localEps.resize(initialParticleSize);
        localEpsInvPow.resize(initialParticleSize);
        localSupports.resize(initialParticleSize);
        calcKernels.swap(accCalcKernels);
        kerOffsets.swap(accKerOffsets);
    }

#ifdef SE_CLASS1
    int getUpdateCtr() const
    {
        return update_ctr;
    }
#endif

    // Here we require the first element of the aggregate to be:
    // 1) the value of the function f on the point
    Dcpse(vector_type &particles,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1,                               //Maybe change this to epsilon/h or h/epsilon = c 0.9. Benchmark
          support_options opt = support_options::RADIUS)
		:particlesFrom(particles),
         particlesTo(particles),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            opt(opt)
    {
        particles.ghost_get_subset();         // This communicates which ghost particles to be excluded from support
        if (supportSizeFactor < 1) 
        {
            initializeAdaptive(particles, particles, convergenceOrder, rCut);
        } 
        else 
        {
            initializeStaticSize(particles, particles, convergenceOrder, rCut, supportSizeFactor);
        }
    }

    //Surface DCPSE Constructor
    template<unsigned int NORMAL_ID>
    Dcpse(vector_type &particles,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T nSpacing,
          value_t< NORMAL_ID >,
          support_options opt = support_options::RADIUS)
		:particlesFrom(particles),
         particlesTo(particles),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            opt(opt),isSurfaceDerivative(true),nSpacing(nSpacing),nCount(floor(rCut/nSpacing))
    {
        particles.ghost_get_subset();         // This communicates which ghost particles to be excluded from support

         if(opt==support_options::ADAPTIVE_SURFACE) {
             supportSizeFactor=nSpacing;
             if(dim==2){
                 nCount=3;
             }
             else{
                 nCount=2;
             }
              SupportBuilder<vector_type,vector_type2>
                supportBuilder(particlesFrom,particlesTo, differentialSignature, rCut, differentialOrder == 0);

                auto it = particlesTo.getDomainIterator();
                while (it.isNext()) {
                    auto key_o = particlesTo.getOriginKey(it.get());
                    Support support = supportBuilder.getSupport(it,1,opt);
                    nSpacings.add(supportBuilder.getLastAvgspacing());
                    ++it;
                  }

         }
        createNormalParticles<NORMAL_ID>(particles);
#ifdef SE_CLASS1
        particles.write("WithNormalParticlesQC");
#endif
        initializeStaticSize(particles, particles, convergenceOrder, rCut, supportSizeFactor);
        accumulateAndDeleteNormalParticles(particles);
    }

    Dcpse(vector_type &particles,
          const Dcpse<dim, vector_type>& other,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1,
          support_options opt = support_options::RADIUS)
        :particlesFrom(particles), particlesTo(particles), opt(opt),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            localSupports(other.localSupports),
            isSharedLocalSupport(true)
    {
        particles.ghost_get_subset();
        if (supportSizeFactor < 1)
            initializeAdaptive(particles, particles, convergenceOrder, rCut);
        else
            initializeStaticSize(particles, particles, convergenceOrder, rCut, supportSizeFactor);
    }
    Dcpse(vector_type &particlesFrom,vector_type2 &particlesTo,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1,
          support_options opt = support_options::RADIUS)
            :particlesFrom(particlesFrom),particlesTo(particlesTo),
             differentialSignature(differentialSignature),
             differentialOrder(Monomial<dim>(differentialSignature).order()),
             monomialBasis(differentialSignature.asArray(), convergenceOrder),
             opt(opt)
    {
        particlesFrom.ghost_get_subset();
        if (supportSizeFactor < 1)
            initializeAdaptive(particlesFrom,particlesTo,convergenceOrder, rCut);
        else
            initializeStaticSize(particlesFrom,particlesTo,convergenceOrder, rCut, supportSizeFactor);
    }


    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k)
    {
        Support support = localSupports.get(k);
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        auto & keys = support.getKeys();
        for (int i = 0 ; i < keys.size() ; i++)
        {
            size_t xqK = keys.get(i);
            particles.template getProp<prp>(xqK) += calcKernels.get(kerOff+i);
        }
    }

    template<unsigned int prp>
    void DrawKernelNN(vector_type &particles, int k)
    {
        Support support = localSupports.get(k);
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        auto & keys = support.getKeys();
        for (int i = 0 ; i < keys.size() ; i++)
        {
            size_t xqK = keys.get(i);
            particles.template getProp<prp>(xqK) = 1.0;
        }
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k, int i)
    {

        Support support = localSupports.get(k);
        size_t xpK = k;
        size_t kerOff = kerOffsets.get(k);
        auto & keys = support.getKeys();
        for (int i = 0 ; i < keys.size() ; i++)
        {
            size_t xqK = keys.get(i);
            particles.template getProp<prp>(xqK)[i] += calcKernels.get(kerOff+i);
        }
    }
    /*
     * breif Particle to Particle Interpolation Evaluation
     */
    template<unsigned int prp1,unsigned int prp2>
    void p2p()
    {
        typedef typename std::remove_reference<decltype(particlesTo.template getProp<prp2>(0))>::type T2;

        auto it = particlesTo.getDomainIterator();
        auto supportsIt = localSupports.begin();
        auto epsItInvPow = localEpsInvPow.begin();
        while (it.isNext()){
            double epsInvPow = *epsItInvPow;
            T2 Dfxp = 0;
            Support support = *supportsIt;
            size_t xpK = support.getReferencePointKey();
            //Point<dim, typename vector_type::stype> xp = particlesTo.getPos(xpK);
            //T fxp = sign * particlesTo.template getProp<fValuePos>(xpK);
            size_t kerOff = kerOffsets.get(xpK);
            auto & keys = support.getKeys();
            for (int i = 0 ; i < keys.size() ; i++)
            {
                size_t xqK = keys.get(i);
                T2 fxq = particlesFrom.template getProp<prp1>(xqK);
                Dfxp += fxq * calcKernels.get(kerOff+i);
            }
            Dfxp = epsInvPow*Dfxp;
            //
            //T trueDfxp = particles.template getProp<2>(xpK);
            // Store Dfxp in the right position
            particlesTo.template getProp<prp2>(xpK) = Dfxp;
            //
            ++it;
            ++supportsIt;
            ++epsItInvPow;
        }
    }
    /*! \brief Save the DCPSE computations
     *
     */
    void save(const std::string &file){
        auto & v_cl=create_vcluster();
        size_t req = 0;

		Packer<decltype(localSupports),HeapMemory>::packRequest(localSupports,req);
        Packer<decltype(localEps),HeapMemory>::packRequest(localEps,req);
        Packer<decltype(localEpsInvPow),HeapMemory>::packRequest(localEpsInvPow,req);
        Packer<decltype(calcKernels),HeapMemory>::packRequest(calcKernels,req);
        Packer<decltype(kerOffsets),HeapMemory>::packRequest(kerOffsets,req);

		// allocate the memory
		HeapMemory pmem;
		//pmem.allocate(req);
		ExtPreAlloc<HeapMemory> mem(req,pmem);

		//Packing
		Pack_stat sts;
		Packer<decltype(localSupports),HeapMemory>::pack(mem,localSupports,sts);
        Packer<decltype(localEps),HeapMemory>::pack(mem,localEps,sts);
        Packer<decltype(localEpsInvPow),HeapMemory>::pack(mem,localEpsInvPow,sts);
        Packer<decltype(calcKernels),HeapMemory>::pack(mem,calcKernels,sts);
        Packer<decltype(kerOffsets),HeapMemory>::pack(mem,kerOffsets,sts);

		// Save into a binary file
	    std::ofstream dump (file+"_"+std::to_string(v_cl.rank()), std::ios::out | std::ios::binary);
	    if (dump.is_open() == false)
        {   std::cerr << __FILE__ << ":" << __LINE__ <<" Unable to write since dump is open at rank "<<v_cl.rank()<<std::endl;
	    	return;
            }
	    dump.write ((const char *)pmem.getPointer(), pmem.size());
	    return;
    }
    /*! \brief Load the DCPSE computations
     *
     *
     */
    void load(const std::string & file)
	{
        auto & v_cl=create_vcluster();
	    std::ifstream fs (file+"_"+std::to_string(v_cl.rank()), std::ios::in | std::ios::binary | std::ios::ate );
	    if (fs.is_open() == false)
	    {
	    	std::cerr << __FILE__ << ":" << __LINE__ << " error, opening file: " << file << std::endl;
	    	return;
	    }

	    // take the size of the file
	    size_t sz = fs.tellg();

	    fs.close();

	    // reopen the file without ios::ate to read
	    std::ifstream input (file+"_"+std::to_string(v_cl.rank()), std::ios::in | std::ios::binary );
	    if (input.is_open() == false)
        {//some message here maybe
	    	return;}

	    // Create the HeapMemory and the ExtPreAlloc memory
	    size_t req = 0;
	    req += sz;
	    HeapMemory pmem;
		ExtPreAlloc<HeapMemory> mem(req,pmem);

		mem.allocate(pmem.size());

		// read
	    input.read((char *)pmem.getPointer(), sz);

	    //close the file
	    input.close();

		//Unpacking
		Unpack_stat ps;
	 	Unpacker<decltype(localSupports),HeapMemory>::unpack(mem,localSupports,ps);
        Unpacker<decltype(localEps),HeapMemory>::unpack(mem,localEps,ps);
        Unpacker<decltype(localEpsInvPow),HeapMemory>::unpack(mem,localEpsInvPow,ps);
        Unpacker<decltype(calcKernels),HeapMemory>::unpack(mem,calcKernels,ps);
        Unpacker<decltype(kerOffsets),HeapMemory>::unpack(mem,kerOffsets,ps);
	 	return;
	}


    void checkMomenta(vector_type &particles)
    {
        openfpm::vector<aggregate<double,double>> momenta;
        openfpm::vector<aggregate<double,double>> momenta_accu;

        momenta.resize(monomialBasis.size());
        momenta_accu.resize(monomialBasis.size());

        for (int i = 0 ; i < momenta.size() ; i++)
        {
            momenta.template get<0>(i) =  3000000000.0;
            momenta.template get<1>(i) = -3000000000.0;
        }

        auto it = particles.getDomainIterator();
        auto supportsIt = localSupports.begin();
        auto epsIt = localEps.begin();
        while (it.isNext())
        {
            double eps = *epsIt;

            for (int i = 0 ; i < momenta.size() ; i++)
            {
                momenta_accu.template get<0>(i) =  0.0;
            }

            Support support = *supportsIt;
            size_t xpK = support.getReferencePointKey();
            Point<dim, T> xp = particles.getPos(support.getReferencePointKey());
            size_t kerOff = kerOffsets.get(xpK);
            auto & keys = support.getKeys();
            for (int i = 0 ; i < keys.size() ; i++)
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

            for (int i = 0 ; i < momenta.size() ; i++)
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

            //
            ++it;
            ++supportsIt;
            ++epsIt;
        }

        for (size_t i = 0 ; i < momenta.size() ; i++)
        {
            std::cout << "MOMENTA: " << monomialBasis.getElement(i) << "Min: " << momenta.template get<0>(i) << "  " << "Max: " << momenta.template get<1>(i) << std::endl;
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

        auto it = particles.getDomainIterator();
        auto supportsIt = localSupports.begin();
        auto epsItInvPow = localEpsInvPow.begin();
        while (it.isNext()) {
            double epsInvPow = *epsItInvPow;

            T Dfxp = 0;
            Support support = *supportsIt;
            size_t xpK = support.getReferencePointKey();
            //Point<dim, typename vector_type::stype> xp = particles.getPos(support.getReferencePointKey());
            T fxp = sign * particles.template getProp<fValuePos>(xpK);
            size_t kerOff = kerOffsets.get(xpK);
            auto & keys = support.getKeys();
            for (int i = 0 ; i < keys.size() ; i++)
            {
                size_t xqK = keys.get(i);
                T fxq = particles.template getProp<fValuePos>(xqK);

                Dfxp += (fxq + fxp) * calcKernels.get(kerOff+i);
            }
            Dfxp *= epsInvPow;
            //
            //T trueDfxp = particles.template getProp<2>(xpK);
            // Store Dfxp in the right position
            particles.template getProp<DfValuePos>(xpK) = Dfxp;
            //
            ++it;
            ++supportsIt;
            ++epsItInvPow;
        }
    }


    /*! \brief Get the number of neighbours
     *
     * \return the number of neighbours
     *
     */
    inline int getNumNN(const vect_dist_key_dx &key)
    {
        return localSupports.get(key.getKey()).size();
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
        return localSupports.get(key.getKey()).getKeys().get(j);
    }


    inline T getSign()
    {
        T sign = 1.0;
        if (differentialOrder % 2 == 0 && differentialOrder!=0) {
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

#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif

        expr_type Dfxp = 0;
        Support support = localSupports.get(key.getKey());
        size_t xpK = support.getReferencePointKey();
        //Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key);
        size_t kerOff = kerOffsets.get(xpK);
        auto & keys = support.getKeys();
        for (int i = 0 ; i < keys.size() ; i++)
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

        //typedef typename decltype(o1.value(key))::blabla blabla;

        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        double eps = localEps.get(key.getKey());
        double epsInvPow = localEpsInvPow.get(key.getKey());

        auto &particles = o1.getVector();
#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif

        expr_type Dfxp = 0;
        Support support = localSupports.get(key.getKey());
        size_t xpK = support.getReferencePointKey();
        //Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key)[i];
        size_t kerOff = kerOffsets.get(xpK);
        auto & keys = support.getKeys();
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


    /**
     * Computes the value of the Surface differential operator for one particle for o1 representing a scalar
     *
     * \param key particle
     * \param o1 source property
     * \return the selected derivative
     *
     */
    template<typename op_type>
    auto computeSurfaceDifferentialOperator(const vect_dist_key_dx &key,
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

#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif

        expr_type Dfxp = 0;
        Support support = localSupports.get(key.getKey());
        size_t xpK = support.getReferencePointKey();
        //Point<dim, T> xp = particles.getPos(xpK);
        expr_type fxp = sign * o1.value(key);
        size_t kerOff = kerOffsets.get(xpK);
        auto & keys = support.getKeys();
        for (int i = 0 ; i < keys.size() ; i++)
        {
            size_t xqK = keys.get(i);
            expr_type fxq = o1.value(vect_dist_key_dx(xqK));
            Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+i);
        }
        //additional contribution of particles normal to reference Particle
        Dfxp = Dfxp + (o1.value(key)+fxp) * calcKernels.get(kerOff+keys.size());
        Dfxp = Dfxp * epsInvPow;
        //
        //T trueDfxp = particles.template getProp<2>(xpK);
        // Store Dfxp in the right position
        return Dfxp;
    }


    void initializeUpdate(vector_type &particlesFrom,vector_type2 &particlesTo)
    {
#ifdef SE_CLASS1
        update_ctr=particlesFrom.getMapCtr();
#endif

        localSupports.clear();
        localEps.clear();
        localEpsInvPow.clear();
        calcKernels.clear();
        kerOffsets.clear();
        initializeStaticSize(particlesFrom,particlesTo, convergenceOrder, rCut, supportSizeFactor);
    }

    void initializeUpdate(vector_type &particles)
    {
#ifdef SE_CLASS1
        update_ctr=particles.getMapCtr();
#endif

        localSupports.clear();
        localEps.clear();
        localEpsInvPow.clear();
        calcKernels.clear();
        kerOffsets.clear();

        initializeStaticSize(particles,particles, convergenceOrder, rCut, supportSizeFactor);
    }

private:

    void initializeAdaptive(vector_type &particlesFrom, 
                            vector_type2 &particlesTo,
                            unsigned int convergenceOrder,
                            T rCut) {
        SupportBuilder<vector_type,vector_type2>
                supportBuilder(particlesFrom, particlesTo, differentialSignature, rCut, differentialOrder == 0);
        unsigned int requiredSupportSize = monomialBasis.size();

        if (!isSharedLocalSupport)
            localSupports.resize(particlesTo.size_local_orig());
        localEps.resize(particlesTo.size_local_orig());
        localEpsInvPow.resize(particlesTo.size_local_orig());
        kerOffsets.resize(particlesTo.size_local_orig());
        kerOffsets.fill(-1);

        auto it = particlesTo.getDomainIterator();
        while (it.isNext()) {
            const T condVTOL = 1e2;
            auto key_o = particlesTo.getOriginKey(it.get());

            if (!isSharedLocalSupport)
                localSupports.get(key_o.getKey()) = supportBuilder.getSupport(it, requiredSupportSize,opt);

            Support& support = localSupports.get(key_o.getKey());

            // Get the points in the support of the DCPSE kernel and store the support for reuse
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());

            // Vandermonde matrix computation
            Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                    vandermonde(support, monomialBasis,particlesFrom, particlesTo,HOverEpsilon);
            vandermonde.getMatrix(V);

            T eps = vandermonde.getEps();

            if (!isSharedLocalSupport) {
                T condV = conditionNumber(V, condVTOL);

                if (condV > condVTOL) {
                    requiredSupportSize *= 2;
                    std::cout
                            << "INFO: Increasing, requiredSupportSize = " << requiredSupportSize
                            << std::endl; // debug
                    continue;
                } else
                    requiredSupportSize = monomialBasis.size();
            }

            localSupports.get(key_o.getKey()) = support;
            localEps.get(key_o.getKey()) = eps;
            localEpsInvPow.get(key_o.getKey()) = 1.0 / openfpm::math::intpowlog(eps,differentialOrder);
            // Compute the diagonal matrix E
            DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
            diagonalScalingMatrix.buildMatrix(E, support, eps, particlesFrom, particlesTo);
            // Compute intermediate matrix B
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E * V;
            // Compute matrix A
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = B.transpose() * B;
            // Compute RHS vector b
            DcpseRhs<dim> rhs(monomialBasis, differentialSignature);
            EMatrix<T, Eigen::Dynamic, 1> b(monomialBasis.size(), 1);
            rhs.template getVector<T>(b);
            // Get the vector where to store the coefficients...
            EMatrix<T, Eigen::Dynamic, 1> a(monomialBasis.size(), 1);
            // ...solve the linear system...
            a = A.colPivHouseholderQr().solve(b);
            // ...and store the solution for later reuse
            kerOffsets.get(key_o.getKey()) = calcKernels.size();

            Point<dim, T> xp = particlesTo.getPosOrig(key_o);

            const auto& support_keys = support.getKeys();
            size_t N = support_keys.size();
            for (size_t i = 0; i < N; ++i)
            {
                const auto& xqK = support_keys.get(i);
                Point<dim, T> xq = particlesFrom.getPosOrig(xqK);
                Point<dim, T> normalizedArg = (xp - xq) / eps;

                calcKernels.add(computeKernel(normalizedArg, a));
            }
            //
            ++it;
        }
    }


    void initializeStaticSize(vector_type &particlesFrom,vector_type2 &particlesTo,
                              unsigned int convergenceOrder,
                              T rCut,
                              T supportSizeFactor) {
#ifdef SE_CLASS1
        this->update_ctr=particlesFrom.getMapCtr();
#endif
        this->rCut=rCut;
        this->supportSizeFactor=supportSizeFactor;
        this->convergenceOrder=convergenceOrder;
        auto & v_cl=create_vcluster();
        if(this->opt==LOAD){
            if(v_cl.rank()==0)
            {std::cout<<"Warning: Creating empty DC-PSE operator! Please use update or load to get kernels."<<std::endl;}
            return;
        }
        SupportBuilder<vector_type,vector_type2>
                supportBuilder(particlesFrom,particlesTo, differentialSignature, rCut, differentialOrder == 0);
        unsigned int requiredSupportSize = monomialBasis.size() * supportSizeFactor;

        if (!isSharedLocalSupport)
            localSupports.resize(particlesTo.size_local_orig());
        localEps.resize(particlesTo.size_local_orig());
        localEpsInvPow.resize(particlesTo.size_local_orig());
        kerOffsets.resize(particlesTo.size_local_orig());
        kerOffsets.fill(-1);
        T avgSpacingGlobal=0;
        size_t Counter=0;
        auto it = particlesTo.getDomainIterator();
        while (it.isNext()) {
            // Get the points in the support of the DCPSE kernel and store the support for reuse
            auto key_o = particlesTo.getOriginKey(it.get());

            if (!isSharedLocalSupport)
                localSupports.get(key_o.getKey()) = supportBuilder.getSupport(it, requiredSupportSize,opt);

            Support& support = localSupports.get(key_o.getKey());

            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());

            // Vandermonde matrix computation
            Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                    vandermonde(support, monomialBasis,particlesFrom,particlesTo,HOverEpsilon);
            vandermonde.getMatrix(V);

            T eps = vandermonde.getEps();
            avgSpacingGlobal+=eps;

            localEps.get(key_o.getKey()) = eps;
            localEpsInvPow.get(key_o.getKey()) = 1.0 / openfpm::math::intpowlog(eps,differentialOrder);
            // Compute the diagonal matrix E
            DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
            diagonalScalingMatrix.buildMatrix(E, support, eps, particlesFrom, particlesTo);
            // Compute intermediate matrix B
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> B = E * V;
            // Compute matrix A
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = B.transpose() * B;

            // Compute RHS vector b
            DcpseRhs<dim> rhs(monomialBasis, differentialSignature);
            EMatrix<T, Eigen::Dynamic, 1> b(monomialBasis.size(), 1);
            rhs.template getVector<T>(b);
            // Get the vector where to store the coefficients...
            EMatrix<T, Eigen::Dynamic, 1> a(monomialBasis.size(), 1);
            // ...solve the linear system...
            a = A.colPivHouseholderQr().solve(b);
            // ...and store the solution for later reuse
            kerOffsets.get(key_o.getKey()) = calcKernels.size();

            Point<dim, T> xp = particlesTo.getPosOrig(key_o);

            const auto& support_keys = support.getKeys();
            size_t N = support_keys.size();
            for (size_t i = 0; i < N; ++i)
            {
                const auto& xqK = support_keys.get(i);
                Point<dim, T> xq = particlesFrom.getPosOrig(xqK);
                Point<dim, T> normalizedArg = (xp - xq) / eps;
                calcKernels.add(computeKernel(normalizedArg, a));
            }
            //
            ++it;
            ++Counter;
        }

        v_cl.sum(avgSpacingGlobal);
        v_cl.sum(Counter);
        v_cl.execute();
        if(v_cl.rank()==0)
        {std::cout<<"DCPSE Operator Construction Complete. The average spacing <h> is: "<<HOverEpsilon*avgSpacingGlobal/(T(Counter))<<". c="<<HOverEpsilon<<std::endl;}


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


    T conditionNumber(const EMatrix<T, -1, -1> &V, T condTOL) const {
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

#endif
#endif //OPENFPM_PDATA_DCPSE_HPP

