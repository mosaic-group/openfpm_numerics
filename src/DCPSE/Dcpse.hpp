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
#include "Vandermonde.hpp"
#include "DcpseDiagonalScalingMatrix.hpp"
#include "DcpseRhs.hpp"
#include "hash_map/hopscotch_map.h"
#include <type_traits>

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


template<unsigned int dim, typename VerletList_type, typename vector_type, typename vector_type2=vector_type>
class Dcpse {
public:
	typedef typename vector_type::stype T;
	typedef typename vector_type::value_type part_type;
	typedef vector_type vtype;

	#ifdef SE_CLASS1
	int update_ctr=0;
	#endif

protected:
	const Point<dim, unsigned int> differentialSignature;
	const unsigned int differentialOrder;
	const MonomialBasis<dim> monomialBasis;

	openfpm::vector<T> localEps; // Each MPI rank has just access to the local ones
	openfpm::vector<T> localEpsInvPow; // Each MPI rank has just access to the local ones

	openfpm::vector<size_t> kerOffsets;
	openfpm::vector<T> calcKernels;
	VerletList_type & verletList;
	vector_type & particlesSupport;
	vector_type2 & particlesDomain;
	T rCut;
	unsigned int convergenceOrder;

	support_options opt;

public:
	// This works in this way:
	// 1) User constructs this by giving a domain of points (where one of the properties is the value of our f),
	//    the signature of the differential operator and the error order bound.
	// 2) The machinery for assembling and solving the linear system for coefficients starts...
	// 3) The user can then call an evaluate(point) method to get the evaluation of the differential operator
	//    on the given point.
	////c=HOverEpsilon. Note that the Eps value is computed by <h>/c (<h>=local average spacing for each particle and its support). This factor c is used in the Vandermonde.hpp.
	T HOverEpsilon=0.9;

#ifdef SE_CLASS1
	int getUpdateCtr() const
	{
		return update_ctr;
	}
#endif

	// Here we require the first element of the aggregate to be:
	// 1) the value of the function f on the point
	Dcpse(
		vector_type& particles,
		VerletList_type& verletList,
		Point<dim, unsigned int> differentialSignature,
		unsigned int convergenceOrder,
		T rCut,
		support_options opt = support_options::RADIUS
	):
		particlesSupport(particles),
		particlesDomain(particles),
		verletList(verletList),
		differentialSignature(differentialSignature),
		differentialOrder(Monomial<dim>(differentialSignature).order()),
		monomialBasis(differentialSignature.asArray(), convergenceOrder),
		opt(opt)
	{
		particles.ghost_get_subset();         // This communicates which ghost particles to be excluded from support
		initializeStaticSize(particles, particles, convergenceOrder, rCut);
	}

	Dcpse(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		Point<dim, unsigned int> differentialSignature,
		unsigned int convergenceOrder,
		T rCut,
		support_options opt = support_options::RADIUS
	):
		particlesSupport(particlesSupport),
		particlesDomain(particlesDomain),
		verletList(verletList),
		differentialSignature(differentialSignature),
		differentialOrder(Monomial<dim>(differentialSignature).order()),
		monomialBasis(differentialSignature.asArray(), convergenceOrder),
		opt(opt)
	{
		particlesSupport.ghost_get_subset();
		initializeStaticSize(particlesSupport,particlesDomain,convergenceOrder, rCut);
	}

	// Default constructor to call from SurfaceDcpse
	// to initialize protected members
	Dcpse(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		Point<dim, unsigned int> differentialSignature,
		unsigned int convergenceOrder,
		support_options opt
	):
		particlesSupport(particlesSupport),
		particlesDomain(particlesDomain),
		verletList(verletList),
		differentialSignature(differentialSignature),
		differentialOrder(Monomial<dim>(differentialSignature).order()),
		monomialBasis(differentialSignature.asArray(), convergenceOrder),
		opt(opt)
	{}

	template<unsigned int prp>
	void DrawKernel(vector_type &particles, int p)
	{
		size_t kerOff = kerOffsets.get(p);
		auto verletIt = this->verletList.getNNIterator(p);
		int i = 0;
		while (verletIt.isNext())
		{
			size_t q = verletIt.get();
			particles.template getProp<prp>(q) += calcKernels.get(kerOff+i);
			++verletIt; ++i;
		}
	}

	template<unsigned int prp>
	void DrawKernelNN(vector_type &particles, int p)
	{
		size_t kerOff = kerOffsets.get(p);
		auto verletIt = this->verletList.getNNIterator(p);
		int i = 0;
		while (verletIt.isNext())
		{
			size_t q = verletIt.get();

			particles.template getProp<prp>(q) = 1.0;
			++verletIt; ++i;
		}
	}

	template<unsigned int prp>
	void DrawKernel(vector_type &particles, int p, int index)
	{

		size_t kerOff = kerOffsets.get(p);
		auto verletIt = this->verletList.getNNIterator(p);
		int i = 0;
		while (verletIt.isNext())
		{
			size_t q = verletIt.get();
			particles.template getProp<prp>(q)[index] += calcKernels.get(kerOff+i);
			++verletIt; ++i;
		}
	}

	/*
	 * breif Particle to Particle Interpolation Evaluation
	 */
	template<unsigned int prp1,unsigned int prp2>
	void p2p()
	{
		typedef typename std::remove_reference<decltype(particlesDomain.template getProp<prp2>(0))>::type T2;

		auto it = particlesDomain.getDomainIterator();
		auto epsItInvPow = localEpsInvPow.begin();
		while (it.isNext()){
			T epsInvPow = *epsItInvPow;
			T2 Dfxp = 0;

			size_t p = it.get();
			auto verletIt = this->verletList.getNNIterator(p);

			size_t kerOff = kerOffsets.get(p);
			int i = 0;
			while (verletIt.isNext())
			{
				size_t q = verletIt.get();
				T2 fxq = particlesSupport.template getProp<prp1>(q);
				Dfxp += fxq * calcKernels.get(kerOff+i);
				++verletIt; ++i;
			}
			Dfxp = epsInvPow*Dfxp;
			particlesDomain.template getProp<prp2>(p) = Dfxp;
			++it;
			++epsItInvPow;
		}
	}

	// foggia 16.09.24
	/*!\fn p2p()
	*
	* \brief Method to perform the particle to particle interpolation of VECTOR fields using DC-PSE kernels.
	*  
	* \tparam propFrom Property ID for the property to interpolate from (vector property, e.g. double[3]).
	* \tparam propTo Property ID for the property to interpolate to (vector property, e.g. double[3]).
	* \tparam N1 Number of elements in the vector property (e.g., for double[3], N1=3).
	*
	*/
	template<unsigned int prp1,unsigned int prp2, size_t N1>
	void p2p()
	{
		typedef typename std::remove_reference<decltype(particlesDomain.template getProp<prp2>(0)[0])>::type T2;

		// Using this one could probably get rid of the N1 tparam. It requires some thought.
		// size_t extent_prp2{std::extent<typename std::remove_reference<decltype(particlesDomain.template getProp<prp2>(0))>::type,0>::value};
		// std::cout << "extent: " << extent_prp2 << std::endl;

		auto it = particlesDomain.getDomainIterator();
		auto epsItInvPow = localEpsInvPow.begin();
		while (it.isNext()){
			double epsInvPow = *epsItInvPow;
			T2 Dfxp[N1];

			for (size_t i1 = 0; i1 < N1; ++i1)
			{
				Dfxp[i1] = 0;
			}

			size_t p = it.get();
			auto verletIt = this->verletList.getNNIterator(p);
			//Point<dim, typename vector_type::stype> xp = particlesDomain.getPos(p);
			//T fxp = sign * particlesDomain.template getProp<fValuePos>(p);
			size_t kerOff = kerOffsets.get(p);

			int i = 0;
			while (verletIt.isNext())
			{
				size_t q = verletIt.get();
				T2 fxq[N1];

				for (size_t i1 = 0; i1 < N1; ++i1)
				{
					fxq[i1] = particlesSupport.template getProp<prp1>(q)[i1];
					Dfxp[i1] += fxq[i1] * calcKernels.get(kerOff+i);
				}

				++verletIt; ++i;
			}

			for (size_t i1 = 0; i1 < N1; ++i1)
			{
				Dfxp[i1] = epsInvPow*Dfxp[i1];
			}

			//
			//T trueDfxp = particles.template getProp<2>(p);
			// Store Dfxp in the right position

			for (size_t i1 = 0; i1 < N1; ++i1)
			{
				particlesDomain.template getProp<prp2>(p)[i1] = Dfxp[i1];
			}

			//
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
		Unpacker<decltype(localEps),HeapMemory>::unpack(mem,localEps,ps);
		Unpacker<decltype(localEpsInvPow),HeapMemory>::unpack(mem,localEpsInvPow,ps);
		Unpacker<decltype(calcKernels),HeapMemory>::unpack(mem,calcKernels,ps);
		Unpacker<decltype(kerOffsets),HeapMemory>::unpack(mem,kerOffsets,ps);
		return;
	}

	void checkMomenta(vector_type &particles)
	{
		openfpm::vector<aggregate<T,T>> momenta;
		openfpm::vector<aggregate<T,T>> momenta_accu;

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
			T eps = *epsIt;

			for (int i = 0 ; i < momenta.size() ; i++)
			{
				momenta_accu.template get<0>(i) =  0.0;
			}

			size_t p = it.get();
			auto verletIt = this->verletList.getNNIterator(p);

			Point<dim, T> xp = particles.getPos(p);
			size_t kerOff = kerOffsets.get(p);
			int i = 0;
			while (verletIt.isNext())
			{
				size_t q = verletIt.get();

				Point<dim, T> xq = particles.getPos(q);
				Point<dim, T> normalizedArg = (xp - xq) / eps;

				auto ker = calcKernels.get(kerOff+i);

				int counter = 0;
				size_t N = monomialBasis.getElements().size();

				for (size_t j = 0; j < N; ++j)
				{
					const Monomial<dim> &m = monomialBasis.getElement(j);

					T mbValue = m.evaluate(normalizedArg);
					momenta_accu.template get<0>(counter) += mbValue * ker;

					++counter;
				}
				++verletIt; ++i;
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
			T epsInvPow = *epsItInvPow;

			T Dfxp = 0;
			size_t p = it.get();
			auto verletIt = this->verletList.getNNIterator(p);
			T fxp = sign * particles.template getProp<fValuePos>(p);
			size_t kerOff = kerOffsets.get(p);
			int i = 0;
			while (verletIt.isNext())
			{
				size_t q = verletIt.get();
				T fxq = particles.template getProp<fValuePos>(q);

				Dfxp += (fxq + fxp) * calcKernels.get(kerOff+i);
				++verletIt; ++i;
			}
			Dfxp *= epsInvPow;
			particles.template getProp<DfValuePos>(p) = Dfxp;

			++it;
			++epsItInvPow;
		}
	}


	/*! \brief Get the number of neighbours
	 *
	 * \return the number of neighbours
	 *
	 */
	inline int getNumNN(const vect_dist_key_dx &p)
	{
		return verletList.getNNPart(p);
	}

	/*! \brief Get the coefficent j (Neighbour) of the particle p
	 *
	 * \param p particle
	 * \param j neighbour
	 *
	 * \return the coefficent
	 *
	 */
	inline T getCoeffNN(const vect_dist_key_dx &p, int j)
	{
		size_t base = kerOffsets.get(p.getKey());
		return calcKernels.get(base + j);
	}

	/*! \brief Get the number of neighbours
 *
 * \return the number of neighbours
 *
 */
	inline size_t getIndexNN(const vect_dist_key_dx &p, int q)
	{
		return verletList.get(p, q);
	}


	inline T getSign()
	{
		T sign = 1.0;
		if (differentialOrder % 2 == 0 && differentialOrder!=0) {
			sign = -1;
		}

		return sign;
	}

	T getEpsilonInvPrefactor(const vect_dist_key_dx &p)
	{
		return localEpsInvPow.get(p.getKey());
	}

	/**
	 * Computes the value of the differential operator for one particle for o1 representing a scalar
	 *
	 * \param p particle
	 * \param o1 source property
	 * \return the selected derivative
	 *
	 */
	template<typename op_type>
	auto computeDifferentialOperator(const vect_dist_key_dx &p,
									 op_type &o1) -> decltype(is_scalar<std::is_fundamental<decltype(o1.value(
			p))>::value>::analyze(p, o1)) {

		typedef decltype(is_scalar<std::is_fundamental<decltype(o1.value(p))>::value>::analyze(p, o1)) expr_type;

		T sign = 1.0;
		if (differentialOrder % 2 == 0) {
			sign = -1;
		}

		T eps = localEps.get(p.getKey());
		T epsInvPow = localEpsInvPow.get(p.getKey());

		auto &particles = o1.getVector();

#ifdef SE_CLASS1
		if(particles.getMapCtr()!=this->getUpdateCtr())
		{
			std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
		}
#endif

		expr_type Dfxp = 0;
		auto verletIt = this->verletList.getNNIterator(p);

		expr_type fxp = sign * o1.value(p);
		size_t kerOff = kerOffsets.get(p);

		int i = 0;
		while (verletIt.isNext())
		{
			size_t q = verletIt.get();

			expr_type fxq = o1.value(vect_dist_key_dx(q));
			Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+i);
			++verletIt; ++i;
		}

		Dfxp = Dfxp * epsInvPow;

		return Dfxp;
	}

	/**
	 * Computes the value of the differential operator for one particle for o1 representing a vector
	 *
	 * \param p particle
	 * \param o1 source property
	 * \param i component
	 * \return the selected derivative
	 *
	 */
	template<typename op_type>
	auto computeDifferentialOperator(
		const vect_dist_key_dx &p,
		op_type &o1,
		int i
	) -> typename decltype(is_scalar<std::is_fundamental<decltype(o1.value(p))>::value>::analyze(p, o1))::coord_type
	{
		typedef typename decltype(is_scalar<std::is_fundamental<decltype(o1.value(p))>::value>::analyze(p, o1))::coord_type expr_type;

		T sign = 1.0;
		if (differentialOrder % 2 == 0) {
			sign = -1;
		}

		T eps = localEps.get(p.getKey());
		T epsInvPow = localEpsInvPow.get(p.getKey());

		auto &particles = o1.getVector();
#ifdef SE_CLASS1
		if(particles.getMapCtr()!=this->getUpdateCtr())
		{
			std::cerr<<__FILE__<<":"<<__LINE__<<" Error: You forgot a DCPSE operator update after map."<<std::endl;
		}
#endif

		expr_type Dfxp = 0;

		expr_type fxp = sign * o1.value(p)[i];
		size_t kerOff = kerOffsets.get(p);

		auto verletIt = this->verletList.getNNIterator(p);
		int j = 0;

		while (verletIt.isNext())
		{
			size_t q = verletIt.get();

			expr_type fxq = o1.value(vect_dist_key_dx(q))[i];
			Dfxp = Dfxp + (fxq + fxp) * calcKernels.get(kerOff+j);
			++verletIt; ++j;
		}

		Dfxp = Dfxp * epsInvPow;

		return Dfxp;
	}

	void initializeUpdate(vector_type &particlesSupport,vector_type2 &particlesDomain)
	{
#ifdef SE_CLASS1
		update_ctr=particlesSupport.getMapCtr();
#endif

		localEps.clear();
		localEpsInvPow.clear();
		calcKernels.clear();
		kerOffsets.clear();

		initializeStaticSize(particlesSupport, particlesDomain, convergenceOrder, rCut);
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

		initializeStaticSize(particles, particles, convergenceOrder, rCut);
	}

protected:
	void initializeStaticSize(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		unsigned int convergenceOrder,
		T rCut
	) {
#ifdef SE_CLASS1
		this->update_ctr=particlesSupport.getMapCtr();
#endif
		this->rCut=rCut;
		this->convergenceOrder=convergenceOrder;
		auto & v_cl=create_vcluster();
#ifdef DCPSE_VERBOSE
		if(this->opt==LOAD){
			if(v_cl.rank()==0)
			{std::cout<<"Warning: Creating empty DC-PSE operator! Please use update or load to get kernels."<<std::endl;}
			return;
		}
#endif
		localEps.resize(particlesDomain.size_local());
		localEpsInvPow.resize(particlesDomain.size_local());
		kerOffsets.resize(particlesDomain.size_local());
		kerOffsets.fill(-1);
		T avgSpacingGlobal=0,avgSpacingGlobal2=0,maxSpacingGlobal=0,minSpacingGlobal=std::numeric_limits<T>::max();
		size_t Counter=0;

		auto it = particlesDomain.getDomainIterator();
		while (it.isNext()) {
			// Get the points in the support of the DCPSE kernel and store the support for reuse
			auto p = it.get();

			auto verletIt = verletList.getNNIterator(p);
			// Vandermonde matrix computation
			Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
				vandermonde(p, verletIt, monomialBasis,particlesSupport,particlesDomain,HOverEpsilon);

			EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(verletList.getNNPart(p), monomialBasis.size());
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

			localEps.get(p.getKey()) = eps;
			localEpsInvPow.get(p.getKey()) = 1.0 / openfpm::math::intpowlog(eps,differentialOrder);
			// Compute the diagonal matrix E
			DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
			EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(verletList.getNNPart(p), verletList.getNNPart(p));

			assert(verletList.getNNPart(p) >= monomialBasis.size());
			assert(E.rows() == verletList.getNNPart(p));
			assert(E.cols() == verletList.getNNPart(p));

			verletIt.reset();
			diagonalScalingMatrix.buildMatrix(E, p, verletIt, eps, particlesSupport, particlesDomain);
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
			kerOffsets.get(p.getKey()) = calcKernels.size();

			Point<dim, T> xp = particlesDomain.getPos(p);

			verletIt.reset();
			unsigned matVRow = 0;
			while (verletIt.isNext())
			{
				size_t q = verletIt.get();
				Point<dim, T> xq = particlesSupport.getPos(q);
				Point<dim, T> x_pqNorm = (xp - xq) / eps;
				calcKernels.add(computeKernel(x_pqNorm, a, V, matVRow, monomialBasis.getElements().size()));
				++verletIt;
				++matVRow;
			}
			++it;
			++Counter;
		}
#ifdef DCPSE_VERBOSE
		v_cl.sum(avgSpacingGlobal);
		v_cl.sum(avgSpacingGlobal2);
		v_cl.max(maxSpacingGlobal);
		v_cl.min(minSpacingGlobal);
		v_cl.sum(Counter);
		v_cl.execute();
		if(v_cl.rank()==0)
		{std::cout<<"DCPSE Operator Construction Complete. The global avg spacing in the support <h> is: "<<HOverEpsilon*avgSpacingGlobal/(T(Counter))<<" (c="<<HOverEpsilon<<"). Avg:"<<avgSpacingGlobal2/(T(Counter))<<" Range:["<<minSpacingGlobal<<","<<maxSpacingGlobal<<"]."<<std::endl;}
#endif
	}

	T computeKernel(Point<dim, T> x, EMatrix<T, Eigen::Dynamic, 1> & a, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>& V, size_t matVRow, size_t monomialBasisSize) const {
		T res = 0;
		T expFactor = exp(-norm2(x));

		for (size_t i = 0; i < monomialBasisSize; ++i)
		{
			T coeff = a(i);
			T mbValue = V(matVRow, i);
			res += coeff * mbValue * expFactor;
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


template<unsigned int dim, typename VerletList_type, typename vector_type, typename vector_type2=vector_type>
class SurfaceDcpse : Dcpse<dim, VerletList_type, vector_type, vector_type2> {
public:
	typedef typename vector_type::stype T;

protected:
	openfpm::vector<size_t> accKerOffsets;
	openfpm::vector<T> accCalcKernels;
	openfpm::vector<T> nSpacings;
	T nSpacing;
	unsigned int nCount;

	bool isSurfaceDerivative=false;
	size_t initialParticleSize;


  /*!\fn createNormalParticless(vector_type &particlesSupport)
   * \brief Normal extension
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.
   * \param particlesSupport Particle set on which the normal extension is done.
   *
   * Creates the particles along the normal of each surface particle.
   */
	template<unsigned int NORMAL_ID>
	void createNormalParticles()
	{
		this->particlesSupport.template ghost_get<NORMAL_ID>(SKIP_LABELLING);
		initialParticleSize=this->particlesSupport.size_local_with_ghost();
		T nSpacing_p = nSpacing;

		auto it = this->particlesSupport.getDomainAndGhostIterator();
		while (it.isNext()) {
			size_t p = it.get();

			Point<dim,T> xp=this->particlesSupport.getPos(p), Normals=this->particlesSupport.template getProp<NORMAL_ID>(p);

			if (this->opt == support_options::ADAPTIVE)
				nSpacing_p = nSpacings.get(p);

			for(int i=1;i<=nCount;i++) {
			  this->particlesSupport.appendLocal();
				for(size_t j=0;j<dim;j++)
					this->particlesSupport.getLastPosEnd()[j] = xp[j]+i*nSpacing_p*Normals[j];

				this->particlesSupport.appendLocal();
				for(size_t j=0;j<dim;j++)
					this->particlesSupport.getLastPosEnd()[j] = xp[j]-i*nSpacing_p*Normals[j];
			}
			++it;
		}

		if (this->opt==support_options::ADAPTIVE)
		  {
		    // Get all rCuts
		    openfpm::vector<T> rCuts;
		    auto it = this->particlesDomain.getDomainIterator();
		    while (it.isNext()) {
		      size_t p = it.get();
		      rCuts.add();
		      rCuts.get(rCuts.size()-1) = this->verletList.getRCuts(p);
		      this->verletList.clear(p); // clear the Verlet list before refilling it
		      ++it;
		    }
#ifdef SE_CLASS1
		    if (rCuts.size() != this->particlesDomain.size_local())
		      {
			std::cerr << __FILE__ << ":" << __LINE__
				  << " ERROR: when updating adaptive cut-off Verlet list in createNormalParticles, rCuts.size() != particlesDomain.size_local(), ["
				  << rCuts.size() << "!=" << this->particlesDomain.size_local() << "]" << std::endl;
			std::runtime_error("Runtime adaptive cut-off Verlet list error");
		      }
#endif
		    // Fill out the Verlet list
		    auto domainIt = this->particlesDomain.getDomainIterator();
		    this->verletList.fillNonSymmAdaptiveIterator(domainIt,
								 this->particlesDomain.getPosVector(),
								 this->particlesSupport.getPosVector(),
								 rCuts,
								 this->particlesDomain.size_local()
								 );
		  }
		else
		{
			this->verletList.initCl(
				this->verletList.getInternalCellList(),
				this->particlesSupport.getPosVector(),
				this->particlesSupport.size_local()
			);

			auto domainIt = this->particlesDomain.getDomainIterator();
			this->verletList.Initialize(
				this->verletList.getInternalCellList(),
				this->verletList.getRCut(),
				domainIt,
				this->particlesDomain.getPosVector(),
				this->particlesSupport.getPosVector(),
				this->particlesDomain.size_local()
			);
		}
	}

  /*!\fn accumulateAndDeleteNormalParticless(vector_type &particlesSupport, vector_type2 &particlesDomain)
   * \brief Creates the Surface DCPSE operator by accumulating the kernel of support particles
   *
   * \param particlesSupport Particle set from which the support is constructed
   * \param particlesDomain Particle set on which the support is constructed
   *
   * This function creates the surface DC-PSE kernels as explained in the original paper.
   * The support of each particle in particlesDomain (which would only have surface particles -- surfDomain-- )
   * is built from the surface particles in particlesSupport -- surfSupport -- and the normally extended particles
   * in particlesSupport -- normalSupport -- (these were created in the createNormalParticles() function.
   * For all particles in particlesDomain we iterate through its support and accumulate/sum all the kernels of the normalSupport particles
   * on the corresponding surfSupport particle.
   * In the end, the normalSupport particles are deleted: the vector is resized to its original size.
   */
	void accumulateAndDeleteNormalParticles()
	{
		accCalcKernels.clear();
		accKerOffsets.clear();
		accKerOffsets.resize(this->particlesDomain.size_local());
		accKerOffsets.fill(-1);

		tsl::hopscotch_map<size_t, size_t> nMap;
		openfpm::vector_std<size_t> supportBuffer; // list of real (surface) particles

		auto it = this->particlesDomain.getDomainIterator();
		while (it.isNext()) {
			supportBuffer.clear();
			nMap.clear();

			size_t p=it.get();
			size_t kerOff = this->kerOffsets.get(p);
			// accumulate kernel offsets of each particle in particlesDomain
			accKerOffsets.get(p)=accCalcKernels.size();
			auto verletIt = this->verletList.getNNIterator(p);
			int i = 0;

			while (verletIt.isNext())
			{
				size_t q = verletIt.get();

				int difference = static_cast<int>(q) - static_cast<int>(initialParticleSize);
				int real_particle;

				// error here, last element is overflown
				// find out whether particle is a real particle (surfSupport) or a virtual normal particle (normalSupport)
				if (std::signbit(difference))
					// it's a real (surface) particle
					real_particle = q;
				else
					// it's not (it's a particle along the normal)
					real_particle = difference / (2 * nCount);

				auto found=nMap.find(real_particle);

				if (found != nMap.end())
					accCalcKernels.get(found->second) += this->calcKernels.get(kerOff+i);
				else
				{
					supportBuffer.add();
					supportBuffer.get(supportBuffer.size()-1) = real_particle;
					accCalcKernels.add();
					accCalcKernels.get(accCalcKernels.size()-1) = this->calcKernels.get(kerOff+i);
					nMap[real_particle]=accCalcKernels.size()-1;
				}

				++verletIt; ++i;
			}

			this->verletList.clear(p);
			
			// Verlet list ends up with only one particle per particle (themselves)
			for (int i = 0; i < supportBuffer.size(); ++i)
				this->verletList.addPart(p, supportBuffer.get(i));

			++it;
		}

		// Delete all normal particles in the particlesSupport
		this->particlesSupport.discardLocalAppend(initialParticleSize);
		// this->localEps.resize(initialParticleSize);
		// this->localEpsInvPow.resize(initialParticleSize);
		// store accumulated kernels (including normalSupport) into surfSupport particles
		this->calcKernels.swap(accCalcKernels);
		this->kerOffsets.swap(accKerOffsets);
	}

public:
  /*!\fn Dcpse
   * \brief Surface DC-PSE constructor
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.
   * \param particlesSupport Particle set from which the support is constructed. The two sets are different only when used for interpolation.
   * \param particlesDomain Particle set on which the operator is constructed. The two sets are different only when used for interpolation.
   * \param verletList Verlet list for the particle set.
   * \param differentialSignature Vector that contains the information on the order of the derivative. For example, df/(dxdy) on a 3D domain: [1,1,0].
   * \param convergenceOrder Convergence order of the numerical operator.
   * \param rCut DEPRECARTED - Size of the support/argument for cell list construction. It has to include sufficient enough particles to create the support.
   * \param nSpacing Spacing of the particlesDomain (on the surface). Used for the default mode, where normal spacing is the same for all surface particles.
   * \param nCount Number of particles along the normal on each side of the surface.
   * \param NORMAL_ID Property ID for the normal field of the particle set.
   * \param opt Type of support. Default: RADIUS, where normal spacing is the same for all surface particles; ADAPTIVE, where the normal spacing is obtaned from the adaptive Verlet list: nSpacing_p = rCut_p/nCount.
   *
   * \note When NOT used for particle to particle interpolation, particlesDomain is set to be the same as particlesSupport.
   */
	template<unsigned int NORMAL_ID>
	SurfaceDcpse(
		vector_type& particlesSupport,
		vector_type2& particlesDomain,
		VerletList_type& verletList,
		Point<dim, unsigned int> differentialSignature,
		unsigned int convergenceOrder,
		T rCut, // TODO: delete this parameter, it is not used
		T nSpacing,
		unsigned int nCount,
		value_t< NORMAL_ID >,
		support_options opt = support_options::RADIUS)
	:
		Dcpse<dim, VerletList_type, vector_type, vector_type2>(particlesSupport, particlesDomain, verletList, differentialSignature, convergenceOrder, opt),
		isSurfaceDerivative(true),
		nSpacing(nSpacing),
		nCount(nCount)
	{
	        particlesSupport.ghost_get_subset(); // TODO: Delete -- This does nothing as that function definition is empty
		this->rCut = rCut;

		if(opt==support_options::ADAPTIVE) {
			// Get the normal spacing for each particle
			nSpacings.clear();
			auto it = particlesDomain.getDomainIterator();
		  	while (it.isNext()) {
		    		size_t p = it.get();
		    		nSpacings.add(verletList.getRCuts(p)/nCount);
				++it;
			}
		}

		if(opt!=support_options::LOAD) {
			createNormalParticles<NORMAL_ID>();
#ifdef SE_CLASS1
			particlesSupport.write("WithNormalParticlesQC"); // TODO: this gives an error in ParaView: the properties of the particles are not there I think.
#endif
		}

		this->initializeStaticSize(
			particlesSupport,
			particlesDomain,
			convergenceOrder,
			rCut
		);

		if(opt!=support_options::LOAD) {
			accumulateAndDeleteNormalParticles();
		}
	}
};

#endif
#endif //OPENFPM_PDATA_DCPSE_HPP

