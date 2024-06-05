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
 * \tparam particlesFrom_type Type of the particle set from which to interpolate.
 * \tparam particlesTo_type Type of the particle set to which to interpolate.
 * \tparam NORMAL_ID Property ID for the normal field of the particle set. If not passed, interpolation is performed on the bulk.
 * 
 * \param particlesFrom Particle set from which to interpolate.
 * \param particlesTo Particle set to which to interpolate.
 * \param ord Convergence order of the numerical operator.
 * \param rCut Size of the support/argument for cell list construction. It has to include sufficient enough particles to create the support.
 * \param oversampling_factor Multiplier to the minimum no. of particles required by the operator in support.
 * \param support_options default:RADIUS (selects all particles inside rCut, overrides oversampling).
 *
 * The interpolation is performed using the (Surface) DC-PSE operators corresponding to the zeroth order derivative.
 * Inside the constructor, the differential signature vector is set to zero, and a Dcpse object is created.
 * Interpolation is then performed when calling the p2p method passing the property ID of the two sets <prop_From,prop_To>.
 */
template<typename particlesFrom_type, typename particlesTo_type, size_t NORMAL_ID = INT_MAX>
class PPInterpolation 
{

    void *dcpse;

    particlesFrom_type & particlesFrom;
    particlesTo_type & particlesTo;
    bool isSurfaceInterpolation=false;
public:
    /*!\fn PPInterpolation
     *
     * \brief Constructor for the bulk particle to particle interpolation.
     */
    PPInterpolation(particlesFrom_type &particlesFrom,particlesTo_type &particlesTo, unsigned int ord, typename particlesFrom_type::stype rCut,
                      double oversampling_factor = dcpse_oversampling_factor,
                      support_options opt = support_options::RADIUS)
    :particlesFrom(particlesFrom),particlesTo(particlesTo)
    {
        Point<particlesFrom_type::dims, unsigned int> p;
        p.zero();
        dcpse = new Dcpse<particlesFrom_type::dims, particlesFrom_type,particlesTo_type>(particlesFrom,particlesTo, p, ord, rCut, oversampling_factor, opt);
    }

  /*!\fn PPInterpolation
   *
   * \brief Constructor for the surface particle to particle interpolation. Only enabled when the property ID of the normal to the surface is passed as the third template parameter. 
   */
  //template<std::enable_if_t< (NORMAL_ID < INT_MAX),int> =0>
  PPInterpolation(particlesFrom_type &particlesFrom,particlesTo_type &particlesTo, unsigned int ord, typename particlesFrom_type::stype rCut,
		  typename particlesFrom_type::stype nSpacing,
		  support_options opt = support_options::RADIUS,
		  bool isSurfaceInterpolation = true)
    :particlesFrom(particlesFrom),particlesTo(particlesTo),isSurfaceInterpolation(true)
  {
    Point<particlesFrom_type::dims, unsigned int> p;
    p.zero();
    dcpse = new Dcpse<particlesFrom_type::dims,particlesFrom_type,particlesTo_type>(particlesFrom,particlesTo, p, ord, rCut, nSpacing,value_t<NORMAL_ID>(), opt);
  }

    void deallocate() {
        delete (Dcpse<particlesFrom_type::dims, particlesFrom_type, particlesTo_type> *) dcpse;
    }

   /* template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse_type<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;
        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }*/

  /*!\fn p2p()
   *
   * \brief Method to perform the particle to particle interpolation using DC-PSE kernels.
   *  
   * \tparam propFrom Property ID for the property to interpolate from.
   * \tparam propTo Property ID for the property to interpolate to.
   *
   */
   template<unsigned int propFrom,unsigned int propTo>
   void p2p() {
       auto dcpse_temp = (Dcpse<particlesFrom_type::dims, particlesFrom_type, particlesTo_type>*) dcpse;
       dcpse_temp->template p2p<propFrom,propTo>();

   }

    // template<unsigned int prp, typename particles_type>
    // void DrawKernel(particles_type &particles, int k) {
    //     auto dcpse_temp = (Dcpse_type<particlesFrom_type::dims, particlesFrom_type, particlesTo_type> *) dcpse;
    //     dcpse_temp->template DrawKernel<prp>(particles, k);

    // }

    // template<unsigned int prp, typename particles_type>
    // void DrawKernelNN(particles_type &particles, int k) {
    //     auto dcpse_temp = (Dcpse_type<particlesFrom_type::dims, particlesFrom_type,particlesTo_type> *) dcpse;
    //     dcpse_temp->template DrawKernelNN<prp>(particles, k);

    // }

    // template<typename particles_type>
    // void checkMomenta(particles_type &particles) {
    //     auto dcpse_temp = (Dcpse_type<particles_type::dims, particlesFrom_type, particlesTo_type> *) dcpse;
    //     dcpse_temp->checkMomenta(particles);

    // }

    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    void update() {
        auto dcpse_temp = (Dcpse<particlesFrom_type::dims, particlesFrom_type, particlesTo_type> *) dcpse;
        dcpse_temp->initializeUpdate(particlesFrom,particlesTo);

    }

};



#endif //OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
