//
// Created by Abhinav Singh on 03.11.21.
//

#ifndef OPENFPM_PDATA_DCPSEINTERPOLATION_HPP
#define OPENFPM_PDATA_DCPSEINTERPOLATION_HPP

/*! \brief Class for Creating the DCPSE Operator For the function approximation objects and computes DCPSE Kernels.
 *
 *
 * \param parts particle set
 * \param ord order of convergence of the operator
 * \param rCut Argument for cell list construction
 * \param oversampling_factor multiplier to the minimum no. of particles required by the operator in support
 * \param support_options default:N_particles, Radius can be used to select all particles inside rCut. Overrides oversampling.
 *
 * \return Operator Dx which is a function on Vector_dist_Expressions
 *
 */
template<typename particlesFrom_type, typename particlesTo_type>
class PPInterpolation 
{

    void *dcpse;

    particlesFrom_type & particlesFrom;
    particlesTo_type & particlesTo;

public:
    /*! \brief Constructor for Creating the DCPSE Operator Dx and objects and computes DCPSE Kernels.
     *
     *
     * \param parts particle set
     * \param ord order of convergence of the operator
     * \param rCut Argument for cell list construction
     * \param oversampling_factor multiplier to the minimum no. of particles required by the operator in support
     * \param support_options default:N_particles, Radius can be used to select all particles inside rCut. Overrides oversampling.
     *
     * \return Operator F which is a function on Vector_dist_Expressions
     *
     */
    PPInterpolation(particlesFrom_type &partsFrom,particlesTo_type &partsTo, unsigned int ord, typename particlesFrom_type::stype rCut,
                      double oversampling_factor = dcpse_oversampling_factor,
                      support_options opt = support_options::RADIUS)
    :particlesFrom(particlesFrom),particlesTo(particlesTo)
    {
        Point<particlesFrom_type::dims, unsigned int> p;
        p.zero();
        dcpse = new Dcpse<particlesFrom_type::dims, particlesFrom_type,particlesTo_type>(partsFrom,partsTo, p, ord, rCut, oversampling_factor, opt);
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

   template<unsigned int prp1,unsigned int prp2>
   void p2p() {
       auto dcpse_temp = (Dcpse<particlesFrom_type::dims, particlesFrom_type, particlesTo_type>*) dcpse;
       dcpse_temp->template p2p<prp1,prp2>();

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
