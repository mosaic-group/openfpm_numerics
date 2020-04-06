//
// Created by tommaso on 29/03/19.
//
#ifndef OPENFPM_PDATA_DCPSE_HPP
#define OPENFPM_PDATA_DCPSE_HPP

#ifdef HAVE_EIGEN

#include "Vector/vector_dist.hpp"
#include "MonomialBasis.hpp"
#include "../../openfpm_numerics/src/DMatrix/EMatrix.hpp"
#include "SupportBuilder.hpp"
#include "Support.hpp"
#include "Vandermonde.hpp"
#include "DcpseDiagonalScalingMatrix.hpp"
#include "DcpseRhs.hpp"

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

template<unsigned int dim, typename vector_type>
class Dcpse {
public:

    typedef typename vector_type::stype T;
    typedef typename vector_type::value_type part_type;
    typedef vector_type vtype;

    // This works in this way:
    // 1) User constructs this by giving a domain of points (where one of the properties is the value of our f),
    //    the signature of the differential operator and the error order bound.
    // 2) The machinery for assembling and solving the linear system for coefficients starts...
    // 3) The user can then call an evaluate(point) method to get the evaluation of the differential operator
    //    on the given point.
private:
    const Point<dim, unsigned int> differentialSignature;
    const unsigned int differentialOrder;
    const MonomialBasis<dim> monomialBasis;
    std::vector<EMatrix<T, Eigen::Dynamic, 1>> localCoefficients; // Each MPI rank has just access to the local ones
    std::vector<Support<dim, T, typename vector_type::value_type>> localSupports; // Each MPI rank has just access to the local ones
    std::vector<T> localEps; // Each MPI rank has just access to the local ones

    vector_type & particles;

    typename vector_type::stype rcut_verlet = -1.0;

public:

    // Here we require the first element of the aggregate to be:
    // 1) the value of the function f on the point
    Dcpse(vector_type &particles,
          Point<dim, unsigned int> differentialSignature,
          unsigned int convergenceOrder,
          T rCut,
          T supportSizeFactor = 1,
          T rcut_verlet = -1.0) :
            particles(particles),
            differentialSignature(differentialSignature),
            differentialOrder(Monomial<dim>(differentialSignature).order()),
            monomialBasis(differentialSignature.asArray(), convergenceOrder),
            rcut_verlet(rcut_verlet)
            {
        if (supportSizeFactor < 1) {
            initializeAdaptive(particles, convergenceOrder, rCut);
        } else {
            initializeStaticSize(particles, convergenceOrder, rCut, supportSizeFactor);
        }
    }



    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k)
    {
        EMatrix<T, Eigen::Dynamic, 1> &a = localCoefficients[k];
        Support<dim, T, part_type> support = localSupports[k];
        auto eps = localEps[k];
        int DrawKernelKounter=0;
        size_t xpK = k;
        Point<dim, T> xp = support.getReferencePoint();
        for (auto &xqK : support.getKeys())
        {
            Point<dim, T> xq = particles.getPos(xqK);
            Point<dim, T> normalizedArg = (xp - xq) / eps;

            particles.template getProp<prp>(xqK) += computeKernel(normalizedArg, a);
            DrawKernelKounter++;
        }
        //std::cout<<"Number of Neighbours: "<<DrawKernelKounter<<std::endl;
    }

    template<unsigned int prp>
    void DrawKernel(vector_type &particles, int k, int i)
    {
        EMatrix<T, Eigen::Dynamic, 1> &a = localCoefficients[k];
        Support<dim, T, part_type> support = localSupports[k];
        auto eps = localEps[k];

        size_t xpK = k;
        Point<dim, T> xp = support.getReferencePoint();
        for (auto &xqK : support.getKeys())
        {
            Point<dim, T> xq = particles.getPos(xqK);
            Point<dim, T> normalizedArg = (xp - xq) / eps;

            particles.template getProp<prp>(xqK)[i] += computeKernel(normalizedArg, a);
        }
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
        auto coefficientsIt = localCoefficients.begin();
        auto supportsIt = localSupports.begin();
        auto epsIt = localEps.begin();
        while (it.isNext())
        {
            double eps = *epsIt;

            for (int i = 0 ; i < momenta.size() ; i++)
            {
                momenta_accu.template get<0>(i) =  0.0;
            }

            Support<dim, T, part_type> support = *supportsIt;
            size_t xpK = support.getReferencePointKey();
            Point<dim, T> xp = support.getReferencePoint();
            for (auto &xqK : support.getKeys())
            {
                Point<dim, T> xq = particles.getPos(xqK);
                Point<dim, T> normalizedArg = (xp - xq) / eps;
                EMatrix<T, Eigen::Dynamic, 1> &a = *coefficientsIt;

                int counter = 0;
                for (const Monomial<dim> &m : monomialBasis.getElements())
                {
                    T mbValue = m.evaluate(normalizedArg);


                    momenta_accu.template get<0>(counter) += mbValue * computeKernel(normalizedArg, a);

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
            ++coefficientsIt;
            ++supportsIt;
            ++epsIt;
        }

        for (int i = 0 ; i < momenta.size() ; i++)
        {
            std::cout << "MOMENTA: " << monomialBasis.getElements()[i] << "Min: " << momenta.template get<0>(i) << "  " << "Max: " << momenta.template get<1>(i) << std::endl;
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
        auto coefficientsIt = localCoefficients.begin();
        auto supportsIt = localSupports.begin();
        auto epsIt = localEps.begin();
        while (it.isNext()) {
            double eps = *epsIt;

            T Dfxp = 0;
            Support<dim, T, part_type> support = *supportsIt;
            size_t xpK = support.getReferencePointKey();
            Point<dim, T> xp = support.getReferencePoint();
            T fxp = sign * particles.template getProp<fValuePos>(xpK);
            for (auto &xqK : support.getKeys()) {
                Point<dim, T> xq = particles.getPos(xqK);
                T fxq = particles.template getProp<fValuePos>(xqK);
                Point<dim, T> normalizedArg = (xp - xq) / eps;
                EMatrix<T, Eigen::Dynamic, 1> &a = *coefficientsIt;
                Dfxp += (fxq + fxp) * computeKernel(normalizedArg, a);
            }
            Dfxp /= pow(eps, differentialOrder);
            //
            T trueDfxp = particles.template getProp<2>(xpK);
            // Store Dfxp in the right position
            particles.template getProp<DfValuePos>(xpK) = Dfxp;
            //
            ++it;
            ++coefficientsIt;
            ++supportsIt;
            ++epsIt;
        }
    }


    /*! \brief Get the number of neighbours
     *
     * \return the number of neighbours
     *
     */
    int getNumNN(const vect_dist_key_dx &key)
    {
        return localSupports[key.getKey()].size();
    }

    /*! \brief Get the coefficent j (Neighbour) of the particle key
     *
     * \param key particle
     * \param j neighbour
     *
     * \return the coefficent
     *
     */
    T getCoeffNN(const vect_dist_key_dx &key, int j)
    {
        double eps = localEps[key.getKey()];

        auto xqK = getIndexNN(key,j);

        Point<dim, T> xp = particles.getPos(key);
        Point<dim, T> xq = particles.getPos(xqK);
        Point<dim, T> normalizedArg = (xp - xq) / eps;

        EMatrix<T, Eigen::Dynamic, 1> &a = localCoefficients[key.getKey()];
        return computeKernel(normalizedArg, a);
    }

    /*! \brief Get the number of neighbours
 *
 * \return the number of neighbours
 *
 */
    size_t getIndexNN(const vect_dist_key_dx &key, int j)
    {
        return localSupports[key.getKey()].getKeys()[j];
    }


    T getSign()
    {
        T sign = 1.0;
        if (differentialOrder % 2 == 0) {
            sign = -1;
        }

        return sign;
    }

    T getEpsilonPrefactor(const vect_dist_key_dx &key)
    {
        double eps = localEps[key.getKey()];

        return pow(eps, differentialOrder);
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

        double eps = localEps[key.getKey()];

        auto &particles = o1.getVector();

        expr_type Dfxp = 0;
        Support<dim, T, part_type> support = localSupports[key.getKey()];
        size_t xpK = support.getReferencePointKey();
        Point<dim, T> xp = support.getReferencePoint();
        expr_type fxp = sign * o1.value(key);
        for (auto &xqK : support.getKeys()) {
            Point<dim, T> xq = particles.getPos(xqK);
            expr_type fxq = o1.value(vect_dist_key_dx(xqK));
            Point<dim, T> normalizedArg = (xp - xq) / eps;
            EMatrix<T, Eigen::Dynamic, 1> &a = localCoefficients[key.getKey()];
            Dfxp = Dfxp + (fxq + fxp) * computeKernel(normalizedArg, a);
        }
        Dfxp = Dfxp / pow(eps, differentialOrder);
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

        double eps = localEps[key.getKey()];

        auto &particles = o1.getVector();

        expr_type Dfxp = 0;
        Support<dim, T, part_type> support = localSupports[key.getKey()];
        size_t xpK = support.getReferencePointKey();
        Point<dim, T> xp = support.getReferencePoint();
        expr_type fxp = sign * o1.value(key)[i];
        for (auto &xqK : support.getKeys()) {
            Point<dim, T> xq = particles.getPos(xqK);
            expr_type fxq = o1.value(vect_dist_key_dx(xqK))[i];
            Point<dim, T> normalizedArg = (xp - xq) / eps;
            EMatrix<T, Eigen::Dynamic, 1> &a = localCoefficients[key.getKey()];
            Dfxp = Dfxp + (fxq + fxp) * computeKernel(normalizedArg, a);
        }
        Dfxp = Dfxp / pow(eps, differentialOrder);
        //
        //T trueDfxp = particles.template getProp<2>(xpK);
        // Store Dfxp in the right position
        return Dfxp;
    }

private:

    void initializeAdaptive(vector_type &particles,
                            unsigned int convergenceOrder,
                            T rCut) {
        SupportBuilder<vector_type>
                supportBuilder(particles, differentialSignature, rCut);
        unsigned int requiredSupportSize = monomialBasis.size();

        auto it = particles.getDomainIterator();
        while (it.isNext()) {
            const T condVTOL = 1e2;

            // Get the points in the support of the DCPSE kernel and store the support for reuse
            Support<dim, T, part_type> support = supportBuilder.getSupport(it, requiredSupportSize,rcut_verlet);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());

            // Vandermonde matrix computation
            Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                    vandermonde(support, monomialBasis);
            vandermonde.getMatrix(V);

            T condV = conditionNumber(V, condVTOL);
            T eps = vandermonde.getEps();

            if (condV > condVTOL) {
                requiredSupportSize *= 2;
                std::cout
                        << "INFO: Increasing, requiredSupportSize = " << requiredSupportSize
                        << std::endl; // debug
                continue;
            } else {
                requiredSupportSize = monomialBasis.size();
            }

            localSupports.push_back(support);
            localEps.push_back(eps);
            // Compute the diagonal matrix E
            DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
            diagonalScalingMatrix.buildMatrix(E, support, eps);
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
            localCoefficients.push_back(a);
            //
            ++it;
        }
    }


    void initializeStaticSize(vector_type &particles,
                              unsigned int convergenceOrder,
                              T rCut,
                              T supportSizeFactor) {
        SupportBuilder<vector_type>
                supportBuilder(particles, differentialSignature, rCut);
        unsigned int requiredSupportSize = monomialBasis.size() * supportSizeFactor;

        auto it = particles.getDomainIterator();
        while (it.isNext()) {
            // Get the points in the support of the DCPSE kernel and store the support for reuse
            Support<dim, T, part_type> support = supportBuilder.getSupport(it, requiredSupportSize,rcut_verlet);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> V(support.size(), monomialBasis.size());
/* Some Debug code
            if (it.get().getKey() == 5564)
            {
                int debug = 0;
                debug++;
            }
*/
            // Vandermonde matrix computation
            Vandermonde<dim, T, EMatrix<T, Eigen::Dynamic, Eigen::Dynamic>>
                    vandermonde(support, monomialBasis);
            vandermonde.getMatrix(V);

            T eps = vandermonde.getEps();

            localSupports.push_back(support);
            localEps.push_back(eps);
            // Compute the diagonal matrix E
            DcpseDiagonalScalingMatrix<dim> diagonalScalingMatrix(monomialBasis);
            EMatrix<T, Eigen::Dynamic, Eigen::Dynamic> E(support.size(), support.size());
            diagonalScalingMatrix.buildMatrix(E, support, eps);
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
            localCoefficients.push_back(a);
            //
            ++it;
        }
    }


    T computeKernel(Point<dim, T> x, EMatrix<T, Eigen::Dynamic, 1> a) const {
        T res = 0;
        unsigned int counter = 0;
        for (const Monomial<dim> &m : monomialBasis.getElements()) {
            T coeff = a(counter);
            T mbValue = m.evaluate(x);
            T expFactor = exp(-norm2(x));
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

