//
// Created by Abhinav Singh on 15.11.21.
//

#ifndef OPENFPM_PDATA_DCPSE_SURFACE_OP_HPP
#define OPENFPM_PDATA_DCPSE_SURFACE_OP_HPP

#include "DCPSE/DCPSE_op/DCPSE_op.hpp"

template<unsigned int NORMAL_ID>
class SurfaceDerivative_x {

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
    SurfaceDerivative_x(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 1;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

        /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }

    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);
    }
};

template<unsigned int NORMAL_ID>
class SurfaceDerivative_y {

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
    SurfaceDerivative_y(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(1) = 1;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }


    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};

template<unsigned int NORMAL_ID>
class SurfaceDerivative_z {

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
    SurfaceDerivative_z(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(2) = 1;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};

template<unsigned int NORMAL_ID>
class Laplace_Beltrami {

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
    Laplace_Beltrami(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 2;
        p.get(1) = 2;
        p.get(2) = 2;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);
    }
};

template<unsigned int NORMAL_ID>
class SurfaceDerivative_xx {

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
    SurfaceDerivative_xx(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 2;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);
    }
};


template<unsigned int NORMAL_ID>
class SurfaceDerivative_yy {

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
    SurfaceDerivative_yy(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(1) = 2;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};

template<unsigned int NORMAL_ID>
class SurfaceDerivative_zz {

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
    SurfaceDerivative_zz(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(2) = 2;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template createNormalParticles<NORMAL_ID>(particles);
        particles.write("With Normal");
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};

template<unsigned int NORMAL_ID>
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
    SurfaceDerivative_xy(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 1;
        p.get(1) = 1;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};


template<unsigned int NORMAL_ID>
class SurfaceDerivative_yz {

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
    SurfaceDerivative_yz(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(1) = 1;
        p.get(2) = 1;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};


template<unsigned int NORMAL_ID>
class SurfaceDerivative_xz {

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
    SurfaceDerivative_xz(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    support_options opt = support_options::RADIUS) {
        Point<particles_type::dims, unsigned int> p;
        p.zero();
        p.get(0) = 1;
        p.get(2) = 1;

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};

template<unsigned int NORMAL_ID>
class SurfaceDerivative_G {

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
    SurfaceDerivative_G(particles_type &parts, unsigned int ord, typename particles_type::stype rCut,typename particles_type::stype nSpacing,
                    const Point<particles_type::dims, unsigned int> &p,support_options opt = support_options::RADIUS) {

        dcpse = new Dcpse<particles_type::dims, particles_type>(parts, p, ord, rCut,nSpacing,value_t<NORMAL_ID>(), opt);
    }

    template<typename particles_type>
    void deallocate(particles_type &parts) {
        delete (Dcpse<particles_type::dims, particles_type> *) dcpse;
    }

    template<typename operand_type>
    vector_dist_expression_op<operand_type, Dcpse<operand_type::vtype::dims, typename operand_type::vtype>, VECT_DCPSE>
    operator()(operand_type arg) {
        typedef Dcpse<operand_type::vtype::dims, typename operand_type::vtype> dcpse_type;

        return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
    }

    template<typename particles_type>
    void checkMomenta(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->checkMomenta(particles);

    }

    template<unsigned int prp, typename particles_type>
    void DrawKernel(particles_type &particles, int k) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->template DrawKernel<prp>(particles, k);

    }

    /*! \brief Method for Saving the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be saved.
     */
    template<typename particles_type>
    void save(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->save(file);
    }
    /*! \brief Method for Loading the DCPSE Operator.
     *
     * \param parts particle set
     * \param file name for data to be loaded from.
     */
    template<typename particles_type>
    void load(particles_type &particles, const std::string &file) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->load(file);
    }
    /*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
     *
     *
     * \param parts particle set
     */
    template<typename particles_type>
    void update(particles_type &particles) {
        auto dcpse_temp = (Dcpse<particles_type::dims, particles_type> *) dcpse;
        dcpse_temp->createNormalParticles<NORMAL_ID>(particles);
        dcpse_temp->initializeUpdate(particles);
        dcpse_temp->accumulateAndDeleteNormalParticles(particles);

    }
};

#endif //OPENFPM_PDATA_DCPSE_SURFACE_OP_HPP
