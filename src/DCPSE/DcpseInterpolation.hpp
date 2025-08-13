//
// Created by Abhinav Singh on 03.11.21.
// Surface interpolation: Lennt Schulze and Ale Foggia on 06.03.24
//

#ifndef OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
#define OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
#include "DCPSE/Dcpse.hpp"

/*!\class PPInterpolation 
 * \brief Class to perform particle to particle interpolation using DC-PSE kernels.
 *
 * \tparam particlesSupport_type Type of the particle set from which to interpolate.
 * \tparam particlesDomain_type Type of the particle set to which to interpolate.
 * \tparam VerletList_type Type of the Verlet List of the particle support set (particlesSupport)
 * \tparam NORMAL_ID Property ID for the normal field of the particle set.
 * 
 * \param particlesFrom Particle set from which to interpolate.
 * \param particlesTo Particle set to which to interpolate.
 * \param ord Convergence order of the numerical operator.
 * \param rCut Size of the support/argument for cell list construction. It has to include sufficient enough particles to create the support.
 * \param isSurfaceInterpolation If not passed as true, interpolation is performed on the bulk.
 * \param oversampling_factor Multiplier to the minimum no. of particles required by the operator in support.
 * \param support_options default:RADIUS (selects all particles inside rCut, overrides oversampling).
 *
 * The interpolation is performed using the (Surface) DC-PSE operators corresponding to the zeroth order derivative.
 * Inside the constructor, the differential signature vector is set to zero, and a Dcpse object is created.
 * Interpolation is then performed when calling the p2p method passing the property ID of the two sets <prop_From,prop_To>.
 */
template<typename particlesSupport_type, typename particlesDomain_type, typename VerletList_type, size_t NORMAL_ID = INT_MAX>
class PPInterpolation 
{

    void *dcpse;

    particlesSupport_type & particlesSupport;
    particlesDomain_type & particlesDomain;


public:
    /*! \brief Constructor for Creating the DCPSE Operator Dx and objects and computes DCPSE Kernels.
     *
     * \param particlesSupport Particle set from which to interpolate.
     * \param particlesDomain Particle set to which to interpolate.
     * \param VerletList Verlet List of the particle support set (particlesSupport)
     * \param ord Convergence order of the numerical operator.
     * \param rCut Size of the support/argument for cell list construction. It has to include sufficient enough particles to create the support.
     * \param support_options default:RADIUS (selects all particles inside rCut, overrides oversampling).
     *
     * \return Operator F which is a function on Vector_dist_Expressions
     */
    PPInterpolation(
        particlesSupport_type &particlesSupport,
        particlesDomain_type &particlesDomain,
        VerletList_type& verletList,
        unsigned int ord,
        typename particlesSupport_type::stype rCut,
        support_options opt = support_options::RADIUS
    ):
        particlesSupport(particlesSupport),
        particlesDomain(particlesDomain)
    {
        Point<particlesSupport_type::dims, unsigned int> p; p.zero();

        dcpse = new Dcpse<
            particlesSupport_type::dims,
            VerletList_type,
            particlesSupport_type,
            particlesDomain_type>
        (
            particlesSupport,
            particlesDomain,
            verletList,
            p,
            ord,
            rCut,
            opt
        );
    }

    /*!\fn PPInterpolation
    * \tparam NORMAL_ID Enables the constructor for surface interpolation operator
    * 
    * \param particlesSupport Particle set from which to interpolate.
    * \param particlesDomain Particle set to which to interpolate.
    * \param VerletList Verlet List of the particle support set (particlesSupport)
    * \param ord Convergence order of the numerical operator.
    * \param rCut Size of the support/argument for cell list construction. It has to include sufficient enough particles to create the support.
    * \param support_options default:RADIUS (selects all particles inside rCut, overrides oversampling).
    *
    * \return Operator F which is a function on Vector_dist_Expressions
    *
    * \brief Constructor for the surface particle to particle interpolation. Only enabled when the property ID of the normal to the surface is passed as the third template parameter. 
    */
    PPInterpolation(
        particlesSupport_type &particlesSupport,
        particlesDomain_type &particlesDomain,
        VerletList_type& verletList,
        unsigned int ord,
        typename particlesSupport_type::stype rCut,
        typename particlesSupport_type::stype nSpacing,
        bool isSurfaceInterpolation,
        support_options opt = support_options::RADIUS
    ):
        particlesSupport(particlesSupport),
        particlesDomain(particlesDomain)
    {
        Point<particlesSupport_type::dims, unsigned int> p;
        p.zero();
        dcpse = new SurfaceDcpse<particlesSupport_type::dims,VerletList_type,particlesSupport_type,particlesDomain_type>(particlesSupport, particlesDomain, verletList, p, ord, rCut, nSpacing, static_cast<unsigned int>(rCut/nSpacing), value_t<NORMAL_ID>(), opt);
    }

    void deallocate() {
        delete (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type> *) dcpse;
    }

   /* template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse_type<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;
        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }*/

  /*!\fn p2p()
   *
   * \brief Method to perform the particle to particle interpolation of SCALAR fields using DC-PSE kernels.
   *  
   * \tparam propSupport Property ID for the property to interpolate from (scalar property, e.g. double).
   * \tparam propDomain Property ID for the property to interpolate to (scalar property, e.g. double).
   *
   */
   template<unsigned int propSupport,unsigned int propDomain>
   void p2p() {
       auto dcpse_temp = (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type>*) dcpse;
       dcpse_temp->template p2p<propSupport,propDomain>();
   }

  // foggia 16.09.24
   /*!\fn p2p()
   *
   * \brief Method to perform the particle to particle interpolation of VECTOR fields using DC-PSE kernels.
   *  
   * \tparam propSupport Property ID for the property to interpolate from (vector property, e.g. double[3]).
   * \tparam propDomain Property ID for the property to interpolate to (vector property, e.g. double[3]).
   * \tparam N1 Number of elements in the vector property (e.g., for double[3], N1=3).
   *
   */
  template<unsigned int propSupport,unsigned int propDomain, size_t N1>
  void p2p() {
    auto dcpse_temp = (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type>*) dcpse;
    dcpse_temp->template p2p<propSupport,propDomain,N1>();
  }

    // template<unsigned int prp, typename particles_type>
    // void DrawKernel(particles_type &particles, int k) {
    //     auto dcpse_temp = (Dcpse_type<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type> *) dcpse;
    //     dcpse_temp->template DrawKernel<prp>(particles, k);

    // }

    // template<unsigned int prp, typename particles_type>
    // void DrawKernelNN(particles_type &particles, int k) {
    //     auto dcpse_temp = (Dcpse_type<particlesSupport_type::dims, particlesSupport_type,particlesDomain_type> *) dcpse;
    //     dcpse_temp->template DrawKernelNN<prp>(particles, k);

    // }

    // template<typename particles_type>
    // void checkMomenta(particles_type &particles) {
    //     auto dcpse_temp = (Dcpse_type<particles_type::dims, particlesSupport_type, particlesDomain_type> *) dcpse;
    //     dcpse_temp->checkMomenta(particles);

    // }

    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     */
    void update() {
        auto dcpse_temp = (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type> *) dcpse;
        dcpse_temp->initializeUpdate(particlesSupport,particlesDomain);
    }
};



#endif //OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
