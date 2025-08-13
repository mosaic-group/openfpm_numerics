//
// Created by Abhinav Singh on 15.11.21.
//

#ifndef OPENFPM_PDATA_DCPSE_SURFACE_OP_HPP
#define OPENFPM_PDATA_DCPSE_SURFACE_OP_HPP
#ifdef HAVE_EIGEN

#include "DCPSE/DCPSE_op/DCPSE_op.hpp"

template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_x {

    void *dcpse;

public:

  /*!\class SurfaceDerivative_x                 
   * \brief Class to create the surface derivative on the x direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
    template<typename particles_type>
    SurfaceDerivative_x(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 1;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

     /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }

    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);
    }
};

  /*!\class SurfaceDerivative_y
   * \brief Class to create the surface derivative on the y direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_y {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_y(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(1) = 1;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }


    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};

  /*!\class SurfaceDerivative_z
   * \brief Class to create the surface derivative on the z direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_z {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_z(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(2) = 1;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};


  /*!\class LaplaceBeltrami
   * \brief Class to create the surface Laplace Beltrami operator (for a scalar field).
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class Laplace_Beltrami {

    void *dcpse;

public:
    template<typename particles_type>
    Laplace_Beltrami(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 2;
        p.get(1) = 2;
        p.get(2) = 2;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);
    }
};

  /*!\class SurfaceDerivative_xx
   * \brief Class to create the second-order surface derivative on the x direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_xx {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_xx(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 2;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);
    }
};

  /*!\class SurfaceDerivative_yy
   * \brief Class to create the second-order surface derivative on the y direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_yy {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_yy(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(1) = 2;
        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};

  /*!\class SurfaceDerivative_zz
   * \brief Class to create the second-order surface derivative on the z direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_zz {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_zz(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(2) = 2;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        particles.write("With Normal");
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};

template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_xy {

    void *dcpse;

public:
    /*! \brief Class for Creating the DCPSE Operator Dxx and objects and computs DCPSE Kernels.
     *
     *
     * \param parts particle set
     * \param ord order of convergence of the operator
     * \param rCut Argument for cell list construction
     * \param oversampling_factor multiplier to the minimum no. of particles required by the operator in support
     * \param support_options default:N_particles, Radius can be used to select all particles inside rCut. Overrides oversampling.
     *
     * \return Operator Dxx which is a function on Vector_dist_Expressions
     *
     */
    template<typename particles_type>
    SurfaceDerivative_xy(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 1;
        p.get(1) = 1;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};

  /*!\class SurfaceDerivative_yz
   * \brief Class to create the second-order surface derivative on the yz direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_yz {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_yz(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(1) = 1;
        p.get(2) = 1;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};

  /*!\class SurfaceDerivative_xz
   * \brief Class to create the second-order surface derivative on the xz direction.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).                                                                                                                
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_xz {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_xz(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        support_options opt = support_options::RADIUS
    ) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 1;
        p.get(2) = 1;

        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};

  /*!\class SurfaceDerivative_G
   * \brief General class to create a surface derivative.
   *
   * \tparam NORMAL_ID Property ID for the normal field of the particle set.                                                                                                     
   * \param parts particle set                       
   * \param ord Convergence order of the numerical operator.
   * \param rCut Size of the support/argument for cell list construction. It has to include enough particles to create the support.                                   
   * \param nSpacing Spacing of the particlesTo (on the surface).
   * \param p Vector that contains the information on the order of the derivative. For example, df/(dxdy) on a 3D domain: [1,1,0].
   * \param opt Type of support.                                                                                                                                                 
   *                                                                                                                                                                             
   * \note The number of particles along the normal is determined as nCount = floor(rCut/nSpacing). The ghost layer has to be at at least as big as rCut.                        
   *                                                                                                                                                                                * \attention If opt = support_options::ADAPTIVE, the meaning of the rCut and nSpacing parameters changes. rCut is a slightly bigger number than the distance between two particlesTo (on the surface). nSpacing is a factor by which rCut is multiplied in order to determine the size of the support.  In this case, the algorithm takes the rCut as a suggestion in order to find the minimum distance in each particle's neighbourhood. Then, to construct the support for the operator, it multiplies that minimum distance by the nSpacing factor. In this case, the number of particles along the normal is set to be (hardcoded) 2 for a 3D problem and 3 for a 2D problem. The ghost layer has to be at least as big as Cut*nSpacing.                                                                                                                                                                     
   *                                                                                                                                                                             
   */
template<unsigned int NORMAL_ID, typename VerletList_type>
class SurfaceDerivative_G {

    void *dcpse;

public:
    template<typename particles_type>
    SurfaceDerivative_G(
        particles_type &parts,
        VerletList_type& verletList,
        unsigned int ord,
        typename particles_type::stype rCut,
        typename particles_type::stype nSpacing,
        unsigned int nCount,
        const Point<particles_type::dims, unsigned int> &p,
        support_options opt = support_options::RADIUS
    ) {
        dcpse = new SurfaceDcpse<particles_type::dims, VerletList_type, particles_type>(parts, parts, verletList, p, ord, rCut,nSpacing, nCount, value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, VerletList_type, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles,particles);

    }
};
#endif //Eigen
#endif //OPENFPM_PDATA_DCPSE_SURFACE_OP_HPP
