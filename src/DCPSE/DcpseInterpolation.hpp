//
// Created by Abhinav Singh on 03.11.21.
//

#ifndef OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
#define OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
#include "DCPSE/Dcpse.hpp"

/*! \brief Class for Creating the DCPSE Operator For the function approximation objects and computes DCPSE Kernels.
 *
 *
 * \param parts particle set
 * \param ord order of convergence of the operator
 * \param rCut Argument for cell list construction
 * \param support_options default: RADIUS to select all particles inside rCut*
 * \return Operator Dx which is a function on Vector_dist_Expressions
 *
 */
template<typename particlesSupport_type, typename particlesDomain_type, typename VerletList_type>
class PPInterpolation 
{

    void *dcpse;

    particlesSupport_type & particlesSupport;
    particlesDomain_type & particlesDomain;

public:
    /*! \brief Constructor for Creating the DCPSE Operator Dx and objects and computes DCPSE Kernels.
     *
     *
     * \param parts particle set
     * \param ord order of convergence of the operator
     * \param rCut Argument for cell list construction
     * \param support_options default: RADIUS to select all particles inside rCut
     *
     * \return Operator F which is a function on Vector_dist_Expressions
     *
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

    void deallocate() {
        delete (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type> *) dcpse;
    }

   /* template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse_type<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;
        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }*/

   template<unsigned int prp1,unsigned int prp2>
   void p2p() {
       auto dcpse_temp = (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type>*) dcpse;
       dcpse_temp->template p2p<prp1,prp2>();

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
     *
     * \param parts particle set
     */
    void update() {
        auto dcpse_temp = (Dcpse<particlesSupport_type::dims, VerletList_type, particlesSupport_type, particlesDomain_type> *) dcpse;
        dcpse_temp->initializeUpdate(particlesSupport,particlesDomain);

    }

};



#endif //OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
