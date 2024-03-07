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
    double rCut,supportSizeFactor=1,nSpacing,AdapFac;
    unsigned int convergenceOrder,nCount;

    bool isSurfaceDerivative=false;
    size_t initialParticleSize;


    support_options opt;
public:

  /*!\fn createNormalParticless(vector_type &particlesFrom)
   * \brief Normal extension
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.
   * \param particlesFrom Particle set on which the normal extension is done.
   *
   * Creates the particles along the normal of each surface particle.
   */
    template<unsigned int NORMAL_ID>
    void createNormalParticles(vector_type &particlesFrom)
    {
        particlesFrom.template ghost_get<NORMAL_ID>(SKIP_LABELLING);
        initialParticleSize=particlesFrom.size_local_with_ghost();
        auto it = particlesFrom.getDomainAndGhostIterator();
        while(it.isNext()){
            auto key=it.get();
            Point<dim,T> xp=particlesFrom.getPos(key), Normals=particlesFrom.template getProp<NORMAL_ID>(key);
            if(opt==support_options::ADAPTIVE)
            {
                nSpacing=nSpacings.get(key.getKey());
            }
            for(int i=1;i<=nCount;i++){
                particlesFrom.addAtEnd();
                for(size_t j=0;j<dim;j++)
                {particlesFrom.getLastPosEnd()[j]=xp[j]+i*nSpacing*Normals[j];}
                particlesFrom.addAtEnd();
                for(size_t j=0;j<dim;j++)
                {particlesFrom.getLastPosEnd()[j]=xp[j]-i*nSpacing*Normals[j];}
            }
            ++it;
        }
    }

  /*!\fn accumulateAndDeleteNormalParticless(vector_type &particlesFrom, vector_type2 &particlesTo)
   * \brief Creates the Surface DCPSE operator by accumulating the kernel of support particles
   *
   * \param particlesFrom Particle set from which the support is constructed
   * \param particlesTo Particle set on which the support is constructed
   *
   * This function creates the surface DC-PSE kernels as explained in the original paper.
   * The support of each particle in particlesTo (which would only have surface particles -- surf_To-- )
   * is built from the surface particles in particlesFrom -- surf_From -- and the normally extended particles
   * in particlesFrom -- normal_From -- (these were created in the createNormalParticles() function.
   * For all particles in particlesTo we iterate through its support and accumulate/sum all the kernels of the normal_From particles
   * on the corresponding surf_From particle.
   * In the end, the normal_From particles are deleted: the vector is resized to its original size.
   */
    void accumulateAndDeleteNormalParticles(vector_type &particlesFrom, vector_type2 &particlesTo)
    {
        tsl::hopscotch_map<size_t, size_t> nMap;
        auto it = particlesTo.getDomainIterator();
        auto particlesTo_supportsIt = localSupports.begin(); // Iterator through the supports of all particles in particleTo
        openfpm::vector_std<size_t> supportBuffer;
        accCalcKernels.clear();
        accKerOffsets.clear();
        accKerOffsets.resize(particlesTo.size_local_orig());
	accKerOffsets.fill(-1);
        while(it.isNext()){
	    supportBuffer.clear();
            nMap.clear();
            auto key=it.get();
            Support particlesTo_support = *particlesTo_supportsIt; // Support of the particleTo with 'key'
            size_t xpK = particlesTo_support.getReferencePointKey();
            size_t kerOff = kerOffsets.get(xpK);
            auto &keys = particlesTo_support.getKeys();            
	    accKerOffsets.get(xpK)=accCalcKernels.size(); // accumulate kernel offsets of each particle in particlesTo

	    // Loop through all particles in the support of the current particle: normal_From and surf_From
            for (int i = 0 ; i < keys.size() ; i++)
            {
                size_t xqK = keys.get(i);		
		int difference = static_cast<int>(xqK) - static_cast<int>(initialParticleSize); // find out whether particle is a real particle (surf_From) or a virtual normal particle (normal_From)
		int real_particle;
		if (std::signbit(difference)) { // it's a real particle
		    real_particle = xqK;
		} else {			// it's not
		    real_particle = difference / (2 * nCount);
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
            keys.swap(supportBuffer); // store keys of the surface support (only surface particles)
            localSupports.get(xpK) = particlesTo_support;
            ++particlesTo_supportsIt;
            ++it;
        }
        particlesFrom.resizeAtEnd(initialParticleSize); // Delete all normal_From particles in the particlesFrom
        localEps.resize(particlesTo.size_local_orig());
        localEpsInvPow.resize(particlesTo.size_local_orig());
        localSupports.resize(particlesTo.size_local_orig());	
        calcKernels.swap(accCalcKernels); // store accumulated kernels (including normal_From) into surf_From particles
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
        initializeStaticSize(particles, particles, convergenceOrder, rCut, supportSizeFactor);
    }

  /*!\fn Dcpse
   * \brief Surface DC-PSE constructor
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.
   * \param particlesFrom Particle set from which the support is constructed. The two sets are different only when used for interpolation.
   * \param particlesTo Particle set on which the operator is constructed. The two sets are different only when used for interpolation.
   * \param differentialSignature Vector that contains the information on the order of the derivative. For example, df/(dxdy) on a 3D domain: [1,1,0].
   * \param convergenceOrder Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include sufficient enough particles to create the support.
   * \param nSpacing Spacing of the particlesTo (on the surface).
   * \param NORMAL_ID Property ID for the normal field of the particle set.
   * \param opt Type of support.
   *
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.
   *
   * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as rCut*nSpacing.
   *
   * \note When NOT used for particle to particle interpolation, particlesTo is set to be the same as particlesFrom.
   */
    template<unsigned int NORMAL_ID>
    Dcpse(vector_type &particlesFrom,vector_type2 &particlesTo,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T nSpacing,
          value_t< NORMAL_ID >,
          support_options opt = support_options::RADIUS)
		:particlesFrom(particlesFrom),
         particlesTo(particlesTo),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            opt(opt),isSurfaceDerivative(true),nSpacing(nSpacing),nCount(floor(rCut/nSpacing))
    {
	    particlesFrom.ghost_get_subset();         // This communicates which ghost particles to be excluded from support

         if(opt==support_options::ADAPTIVE) {
             this->AdapFac=nSpacing;
             if(dim==2){
                 nCount=3;
             }
             else{
                 nCount=2;
             }
             SupportBuilder<vector_type,vector_type2>
                supportBuilder(particlesFrom,particlesTo, differentialSignature, rCut, differentialOrder == 0);
                supportBuilder.setAdapFac(nSpacing);
                auto it = particlesTo.getDomainAndGhostIterator();
                while (it.isNext()) {
                    auto key_o = particlesTo.getOriginKey(it.get());
                    Support support = supportBuilder.getSupport(it,monomialBasis.size(),opt);
                    nSpacings.add(supportBuilder.getLastMinspacing());
                    ++it;
                  }

         }
         if(opt!=support_options::LOAD) {
             createNormalParticles<NORMAL_ID>(particlesFrom);
#ifdef SE_CLASS1
             particlesFrom.write("WithNormalParticlesQC");
#endif
         }
        initializeStaticSize(particlesFrom, particlesTo, convergenceOrder, rCut, supportSizeFactor);
         if(opt!=support_options::LOAD) {
             accumulateAndDeleteNormalParticles(particlesFrom, particlesTo);
         }
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
        supportBuilder.setAdapFac(AdapFac);

        if (!isSharedLocalSupport)
            localSupports.resize(particlesTo.size_local_orig());
        localEps.resize(particlesTo.size_local_orig());
        localEpsInvPow.resize(particlesTo.size_local_orig());
        kerOffsets.resize(particlesTo.size_local_orig());
        kerOffsets.fill(-1);
        T avgSpacingGlobal=0,avgSpacingGlobal2=0,maxSpacingGlobal=0,minSpacingGlobal=std::numeric_limits<T>::max();
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
            T tSpacing = vandermonde.getMinSpacing();
            avgSpacingGlobal2+=tSpacing;
            if(tSpacing>maxSpacingGlobal)
            {
                maxSpacingGlobal=tSpacing;
            }
            if(tSpacing<minSpacingGlobal)
            {
                minSpacingGlobal=tSpacing;
            }

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
        v_cl.sum(avgSpacingGlobal2);
        v_cl.max(maxSpacingGlobal);
        v_cl.min(minSpacingGlobal);
        v_cl.sum(Counter);
        v_cl.execute();
        if(v_cl.rank()==0)
        {std::cout<<"DCPSE Operator Construction Complete. The global avg spacing in the support <h> is: "<<HOverEpsilon*avgSpacingGlobal/(T(Counter))<<" (c="<<HOverEpsilon<<"). Avg:"<<avgSpacingGlobal2/(T(Counter))<<" Range:["<<minSpacingGlobal<<","<<maxSpacingGlobal<<"]."<<std::endl;}
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

