//
// DCPSE Created by tommaso on 29/03/19.
// Modified, Updated and Maintained by Abhinav and Pietro
// Surface Operators by Abhinav Singh on 07/10/2021 and update to use verlet list Hackathon24.
#ifndef OPENFPM_PDATA_DCPSE_HPP
#define OPENFPM_PDATA_DCPSE_HPP

#ifdef HAVE_EIGEN

#include "Vector/vector_dist.hpp"
#include "MonomialBasis.hpp"
#include "DMatrix/EMatrix.hpp"
#include "DcpseRhs.hpp"
#include "hash_map/hopscotch_map.h"

template<unsigned int N> struct value_t {};
enum support_option
{
    CONSTRUCT,
    LOAD
};

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
    typedef decltype(std::declval<vector_type>().getVerlet(0.0)) list_type;
    typedef decltype(std::declval<vector_type>().getVerlet(0.0).getNNIterator(0)) nbd_it_type;
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

    openfpm::vector<T> localEps; // Each MPI rank has just access to the local ones
    openfpm::vector<T> localEpsInvPow; // Each MPI rank has just access to the local ones

    openfpm::vector<size_t> kerOffsets,accKerOffsets;
    openfpm::vector<T> calcKernels;
    openfpm::vector<T> accCalcKernels;
    openfpm::vector<T> nSpacings;
    vector_type & particlesFrom;
    vector_type2 & particlesTo;
    list_type &verletList;
    double nSpacing;
    unsigned int convergenceOrder,nCount;
    bool isSurfaceDerivative=false;
    size_t initialParticleSize;
    support_option opt;

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
            //if(opt==support_option::ADAPTIVE)
            //{
            //    nSpacing=nSpacings.get(key.getKey());
            //}
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
        openfpm::vector_std<size_t> supportBuffer;
        accCalcKernels.clear();
        accKerOffsets.clear();
        accKerOffsets.resize(initialParticleSize);
        accKerOffsets.fill(-1);
        while(it.isNext()){
            supportBuffer.clear();
            nMap.clear();
            auto key=it.get();
            auto xpK = it.get();
            size_t kerOff = kerOffsets.get(xpK.getKey());
            accKerOffsets.get(xpK.getKey())=accCalcKernels.size();
            auto itNN = verletList.getNNIterator(xpK.getKey());
            size_t i=0;
            while (itNN.isNext()) {
                auto xqK = itNN.get();
                if(xpK.getKey()!=xqK) {
                    int difference = static_cast<int>(xqK) - static_cast<int>(initialParticleSize);
                    int real_particle;
                    if (difference<0) {
                        real_particle = xqK;
                    } else {
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
                            ++i;
                }
                ++itNN;
            }
            verletList.replace_p(xpK.getKey(), supportBuffer);
            ++it;
        }
        particles.resizeAtEnd(initialParticleSize);
        localEps.resize(initialParticleSize);
        localEpsInvPow.resize(initialParticleSize);
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
          list_type &verletList,
          support_option opt = support_option::CONSTRUCT)
		:particlesFrom(particles),
         particlesTo(particles),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            verletList(verletList),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            opt(opt)
    {
        particles.ghost_get_subset();         // This communicates which ghost particles to be excluded from support
        initializeStaticSize(particles, particles);
    }

    //Surface DCPSE Constructor
    template<unsigned int NORMAL_ID>
    Dcpse(vector_type &particles,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          list_type &verletList,
          T nSpacing, T nCount,
          value_t< NORMAL_ID >,
          support_option opt = support_option::CONSTRUCT)
		:particlesFrom(particles),
         particlesTo(particles),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            verletList(verletList),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            opt(opt),isSurfaceDerivative(true),nSpacing(nSpacing),nCount(nCount)
    {
        particles.ghost_get_subset();         // This communicates which ghost particles to be excluded from support

         /*if(opt==support_option::ADAPTIVE) {
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
        */
         if(opt!=support_option::LOAD) {
             createNormalParticles<NORMAL_ID>(particles);
#ifdef SE_CLASS1
             particles.write("WithNormalParticlesQC");
#endif
         }
        initializeStaticSize(particles, particles);
         if(opt!=support_option::LOAD) {
             accumulateAndDeleteNormalParticles(particles);
         }
    }
    /* breif Special constructor for vectorial expression
     *
     */
    Dcpse(vector_type &particles,
          const Dcpse<dim, vector_type>& other,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          list_type &verletList,
          support_option opt = support_option::CONSTRUCT)
        :particlesFrom(particles), particlesTo(particles), opt(opt),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            verletList(verletList),
            monomialBasis(differentialSignature.asArray(), convergenceOrder)
    {
        particles.ghost_get_subset();
        initializeStaticSize(particles, particles);
    }

    Dcpse(vector_type &particlesFrom,vector_type2 &particlesTo,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          list_type &verletList,
          support_option opt = support_option::CONSTRUCT)
            :particlesFrom(particlesFrom),particlesTo(particlesTo),
             differentialSignature(differentialSignature),
             differentialOrder(Monomial<dim>(differentialSignature).order()),
             verletList(verletList),
             monomialBasis(differentialSignature.asArray(), convergenceOrder),
             opt(opt)
    {
        particlesFrom.ghost_get_subset();
        initializeStaticSize(particlesFrom,particlesTo,differentialSignature==0);
    }


    template<unsigned int prp>
    void DrawKernel(vector_type &particles, const vect_dist_key_dx &xpK)
    {
        size_t kerOff = kerOffsets.get(xpK.getKey());
        auto itNN = verletList.getNNIterator(xpK.getKey());
            size_t i=0;
            while (itNN.isNext()) {
                auto xqK = itNN.get();
                if(xpK.getKey()!=xqK) {
                    particles.template getProp<prp>(xqK) += calcKernels.get(kerOff+i);
                    ++i;
                }
                ++itNN;
            }
    }

    template<unsigned int prp>
    void DrawKernelNN(vector_type &particles, const vect_dist_key_dx &xpK)
    {
        size_t kerOff = kerOffsets.get(xpK.getKey());
        auto itNN = verletList.getNNIterator(xpK.getKey());
            size_t i=0;
            while (itNN.isNext()) {
                auto xqK = itNN.get();
                if(xpK.getKey()!=xqK) {
                    particles.template getProp<prp>(xqK) = 1.0;
                    ++i;
                }
                ++itNN;
        }
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles,const vect_dist_key_dx &xpK, const size_t &i)
    {
        size_t kerOff = kerOffsets.get(xpK.getKey());
        auto itNN = verletList.getNNIterator(xpK.getKey());
        size_t j=0;
        while(itNN.isNext()) {
                auto xqK = itNN.get();
                if(xpK.getKey()!=xqK) {
                     particles.template getProp<prp>(xqK)[i] += calcKernels.get(kerOff+j);
                     ++j;
                }
                ++itNN;
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
        auto epsItInvPow = localEpsInvPow.begin();
        while (it.isNext()){
            double epsInvPow = *epsItInvPow;
            T2 Dfxp = 0;
            auto xpK = it.get();
            size_t kerOff = kerOffsets.get(xpK.getKey());
            auto itNN = verletList.getNNIterator(xpK.getKey());
            size_t i=0;
            while(itNN.isNext()) {
                auto xqK = itNN.get();
                T2 fxq = particlesFrom.template getProp<prp1>(xqK);
                Dfxp += fxq * calcKernels.get(kerOff + i);
                ++i;
                ++itNN;
            }
            Dfxp = epsInvPow*Dfxp;
            particlesTo.template getProp<prp2>(xpK) = Dfxp;
            ++it;
            ++epsItInvPow;
        }
    }
    /*! \brief Save the DCPSE computations
     *
     */
    void save(const std::string &file){
        auto & v_cl=create_vcluster();
        size_t req = 0;
        Packer<typename std::remove_reference<decltype(verletList)>::type,HeapMemory>::packRequest(verletList,req);
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
        Packer<typename std::remove_reference<decltype(verletList)>::type,HeapMemory>::pack(mem,verletList,sts);
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
        Unpacker<typename std::remove_reference<decltype(verletList)>::type,HeapMemory>::unpack(mem,verletList,ps);
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
        auto epsIt = localEps.begin();
        while (it.isNext())
        {
            double eps = *epsIt;
            for (int i = 0 ; i < momenta.size() ; i++)
            {
                momenta_accu.template get<0>(i) =  0.0;
            }
            auto key=it.get();
            size_t xpK = key.getKey();
            Point<dim, T> xp = particles.getPos(key);
            size_t kerOff = kerOffsets.get(xpK);
            //NN iterator
            auto itNN = verletList.getNNIterator(xpK);
            size_t i=0;
            while (itNN.isNext()) {
                auto xqK = itNN.get();
                if(xpK !=xqK) {
                    Point<dim, T> xq = particles.getPos(xqK);
                    Point<dim, T> normalizedArg = (xp - xq) / eps;
                    auto ker = calcKernels.get(kerOff+i);
                    int counter = 0;
                    size_t N = monomialBasis.getElements().size();
                    for (size_t j = 0; j < N; ++j) {
                        const Monomial<dim> &m = monomialBasis.getElement(j);
                        T mbValue = m.evaluate(normalizedArg);
                        momenta_accu.template get<0>(counter) += mbValue * ker;
                        ++counter;
                    }
                    ++i;
                }
                ++itNN;
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
            ++it;
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
        auto epsItInvPow = localEpsInvPow.begin();
        while (it.isNext()) {
            auto xpK=it.get();
            double epsInvPow = *epsItInvPow;
            T Dfxp = 0;
            T fxp = sign * particles.template getProp<fValuePos>(xpK);
            size_t kerOff = kerOffsets.get(xpK.getKey());
            auto itNN = verletList.getNNIterator(xpK.getKey());
            size_t i=0;
            while (itNN.isNext()) {
                auto xqK = itNN.get();
                if(xpK.getKey()!=xqK) {
                    T fxq = particles.template getProp<fValuePos>(xqK);
                    Dfxp += (fxq + fxp) * calcKernels.get(kerOff+i);
                    ++i;
                }
                ++itNN;
            }
            Dfxp *= epsInvPow;
            particles.template getProp<DfValuePos>(xpK) = Dfxp;
            //
            ++it;
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
        return verletList.getNNPart(key.getKey())-1;
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
    inline size_t getIndexNN(const vect_dist_key_dx &i,const vect_dist_key_dx &j) const
    {

        return verletList.getNeighborId(i.getKey(),j.getKey());
        //localSupports.get(key.getKey()).getKeys().get(j);
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
    auto computeDifferentialOperator(const vect_dist_key_dx &xpK,
                                     op_type &o1) -> decltype(is_scalar<std::is_fundamental<decltype(o1.value(
            xpK))>::value>::analyze(xpK, o1)) {

        typedef decltype(is_scalar<std::is_fundamental<decltype(o1.value(xpK))>::value>::analyze(xpK, o1)) expr_type;

        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        double eps = localEps.get(xpK.getKey());
        double epsInvPow = localEpsInvPow.get(xpK.getKey());

        auto &particles = o1.getVector();

#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif

        expr_type Dfxp = 0;
        expr_type fxp = sign * o1.value(xpK);
        size_t kerOff = kerOffsets.get(xpK.getKey());
        auto itNN = verletList.getNNIterator(xpK.getKey());
        size_t i=0;
        while (itNN.isNext()) {
            auto xqK = itNN.get();
            if(xpK.getKey()!=xqK) {
                expr_type fxq = o1.value(vect_dist_key_dx(xqK));
                Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+i);
                ++i;
            }
            ++itNN;
        }
        Dfxp = Dfxp * epsInvPow;
        return Dfxp;
    }

    /**
     * Computes the value of the differential operator for one particle for o1 representing a vector with components
     *
     * \param key particle
     * \param o1 source property
     * \param i component
     * \return the selected derivative
     *
     */
    template<typename op_type>
    auto computeDifferentialOperator(const vect_dist_key_dx &xpK,
                                     op_type &o1,
                                     int i) -> typename decltype(is_scalar<std::is_fundamental<decltype(o1.value(
            xpK))>::value>::analyze(xpK, o1))::coord_type {

        typedef typename decltype(is_scalar<std::is_fundamental<decltype(o1.value(xpK))>::value>::analyze(xpK, o1))::coord_type expr_type;
        //typedef typename decltype(o1.value(xpK))::blabla blabla; // This is used for understanding and debugging templated code.
        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }
        double eps = localEps.get(xpK.getKey());
        double epsInvPow = localEpsInvPow.get(xpK.getKey());
        auto &particles = o1.getVector();
#ifdef SE_CLASS1
        if(particles.getMapCtr()!=this->getUpdateCtr())
        {
            std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
        }
#endif
        expr_type Dfxp = 0;
        expr_type fxp = sign * o1.value(xpK)[i];
        size_t kerOff = kerOffsets.get(xpK.getKey());
        auto itNN = verletList.getNNIterator(xpK.getKey());
        size_t j=0;
        while (itNN.isNext()) {
            auto xqK = itNN.get();
            if(xpK.getKey()!=xqK) {
                expr_type fxq = o1.value(vect_dist_key_dx(xqK))[i];
                Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+j);
                ++j;
            }
            ++itNN;
        }
        Dfxp = Dfxp * epsInvPow;
        return Dfxp;
    }

    void initializeUpdate(vector_type &particlesFrom,vector_type2 &particlesTo)
    {
#ifdef SE_CLASS1
        update_ctr=particlesFrom.getMapCtr();
#endif

        localEps.clear();
        localEpsInvPow.clear();
        calcKernels.clear();
        kerOffsets.clear();
        initializeStaticSize(particlesFrom,particlesTo, convergenceOrder);
    }

    void initializeUpdate(vector_type &particles)
    {
#ifdef SE_CLASS1
        update_ctr=particles.getMapCtr();
#endif
        localEps.clear();
        localEpsInvPow.clear();
        calcKernels.clear();
        kerOffsets.clear();

        initializeStaticSize(particles,particles);
    }

private:
    void initializeStaticSize(vector_type &particlesFrom,vector_type2 &particlesTo,
                              const bool &isInterpolation=0) {
#ifdef SE_CLASS1
        this->update_ctr=particlesFrom.getMapCtr();
#endif
        auto & v_cl=create_vcluster();
        if(this->opt==LOAD){
            if(v_cl.rank()==0)
            {std::cout<<"Warning: Creating empty DC-PSE operator! Please use update or load to get kernels."<<std::endl;}
            return;
        }

        //SupportBuilder<vector_type,vector_type2>
        //        supportBuilder(particlesFrom,particlesTo, differentialSignature, rCut, differentialOrder == 0);
        unsigned int requiredSupportSize = monomialBasis.size();

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

            size_t nnbs = verletList.getNNPart(key_o.getKey()); //no. of neighbours minus center particle
            if(!isInterpolation)
                {nnbs-=1;}
            // First check that the number of points given is enough for building the Vandermonde matrix
            if (nnbs < monomialBasis.size()) {
                ACTION_ON_ERROR(
                        std::length_error("Not enough neighbour points passed for Vandermonde matrix construction!"));
            }

            //First compute eps and min spacing.
            T avgNeighbourSpacing = 0;
            T minSpacing = std::numeric_limits<T>::max();
            {
                auto itNN = verletList.getNNIterator(key_o.getKey());
                while (itNN.isNext()) {
                    auto key = itNN.get();
                    if(key_o.getKey()==key && !isInterpolation){
                        ++itNN;
                        continue;
                    }
                    Point<dim, T> xp = particlesFrom.getPosOrig(key_o);
                    Point<dim, T> xq = particlesTo.getPosOrig(key);
                    Point<dim, T> Arg = xp-xq; //Possible OpenFPM bug here. The subtraction is not working properly if directly using Point. particlesFrom.getPosOrig(key_o)-particlesTo.getPosOrig(key)
                    double dist = norm(Arg);
                    avgNeighbourSpacing += Arg.norm1();
                    if (minSpacing > dist) { minSpacing = dist; }
                    ++itNN;
                }
            }
            avgNeighbourSpacing/=T(nnbs);
            T eps = avgNeighbourSpacing/HOverEpsilon;
            if (eps==0) {
                ACTION_ON_ERROR(std::length_error("Average neighbour spacing is for particle: "));
                //ACTION_ON_ERROR(std::cout<<key_o.getKey()<<std::endl);
            }

            //Now we build vandermonde V and Diagonal Scaling matrix E on the fly.
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(nnbs, monomialBasis.size());
            Eigen::DiagonalMatrix<T, Eigen::Dynamic> E(nnbs);
#ifdef SE_CLASS1
            //Extra Checks in SE_CLASS1 to see if Neighbors are enough
            assert(nnbs >= monomialBasis.size());
#endif
            {
                //here we precompute a smart dictionary for avoiding recomputations of pow.
                T epsSq = (2.0 * eps * eps);
                auto& basisElements = monomialBasis.getElements();
                auto itNN = verletList.getNNIterator(key_o.getKey());
                size_t i = 0;
                while (itNN.isNext()) {
                    auto key = itNN.get();
                    if(key_o.getKey()==key && !isInterpolation){
                        ++itNN;
                        continue;
                    }
                    Point<dim, T> xp = particlesFrom.getPosOrig(key_o);
                    Point<dim, T> xq = particlesTo.getPosOrig(key);
                    Point<dim, T> Arg = xp-xq;
                    double dist = norm(Arg);
                    for (size_t j = 0; j < basisElements.size(); ++j) {
                        const Monomial<dim> &m =  basisElements.get(j);
                        //double temp=m.evaluate(Arg);
                        V(i, j) = m.evaluate(Arg);
                    }
                    //Eigen access diagonal matrix entry i
                    E.diagonal()[i]=exp(- norm2(Arg) / epsSq);
                    ++i;
                    ++itNN;
                }
                //computing columnwise prefactor for as many basis elements
                for (size_t j = 0; j < basisElements.size(); ++j) {
                    const Monomial<dim> &m =  basisElements.get(j);
                    V.col(j) /= openfpm::math::intpowlog(eps, m.order());
                }
            }
            //Compute Spacing Statistics for convergence and output.
            avgNeighbourSpacing/=T(nnbs);
            avgSpacingGlobal+=eps;
            T tSpacing = minSpacing;
            avgSpacingGlobal2+=tSpacing;
            if(tSpacing>maxSpacingGlobal)
            {
                maxSpacingGlobal=tSpacing;
            }
            if(tSpacing<minSpacingGlobal)
            {
                minSpacingGlobal=tSpacing;
            }
            //Store the computed eps and its inverse power for later use.
            localEps.get(key_o.getKey()) = eps;
            localEpsInvPow.get(key_o.getKey()) = 1.0 / openfpm::math::intpowlog(eps,differentialOrder);
            // Compute matrix A, Note that in dcpse, intermediate B = E * V.
            //EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> A = (E*V).transpose()*(E*V);
            // Compute RHS vector b
            EMatrix<T, Eigen::Dynamic, 1> b(monomialBasis.size(), 1), a(monomialBasis.size(), 1);
            DcpseRhs<dim> rhs(monomialBasis, differentialSignature);
            rhs.template getVector<T>(b);
            // ...solve the linear system...
            //a = ((E*V).transpose()*(E*V)).colPivHouseholderQr().solve(b);
            a = ((E*V).transpose()*(E*V)).lu().solve(b);
            // ...and store the solution for later reuse
            kerOffsets.get(key_o.getKey()) = calcKernels.size();

            Point<dim, T> xp = particlesTo.getPosOrig(key_o);
            {
                auto itNN = verletList.getNNIterator(key_o.getKey());
                while (itNN.isNext()) {
                    auto xqK = itNN.get();
                    if(key_o.getKey()==xqK && !isInterpolation){
                        ++itNN;
                        continue;
                    }
                    Point<dim, T> xq = particlesTo.getPosOrig(xqK);
                    Point<dim, T> normalizedArg = (xp-xq)/ eps;
                    calcKernels.add(computeKernel(normalizedArg, a));
                    ++itNN;
                }
            }
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

