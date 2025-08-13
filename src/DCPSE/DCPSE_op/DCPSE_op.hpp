/*
 * DCPSE_op.hpp
 *
 *  Created on: Jan 7, 2020
 *      Author: Abhinav Singh, Pietro Incardona
 */

#ifndef DCPSE_OP_HPP_
#define DCPSE_OP_HPP_
#ifdef HAVE_EIGEN

#include "Decomposition/CartDecomposition.hpp"
#include "DCPSE/Dcpse.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#if defined(__NVCC__)
#include "DCPSE/Dcpse.cuh"
#endif


const double dcpse_oversampling_factor = 1.9;
const double rcut_verlet = 3.1;

/*! \brief Subtraction operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template<typename exp1, typename DCPSE_type>
class vector_dist_expression_op<exp1, DCPSE_type, VECT_DCPSE> {
	//! expression 1
	const exp1 o1;

	DCPSE_type &dcp;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	typedef typename exp1::vtype vtype;

	//! Costruct a subtraction expression out of two expressions
	inline vector_dist_expression_op(const exp1 &o1, DCPSE_type &dcp)
			: o1(o1), dcp(dcp) {}

	/*! \brief This function must be called before value
	*
	* it initialize the expression if needed
	*
	 */
	inline void init() const {
		o1.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(vect_dist_key_dx()))>::type>
	inline r_type value(const vect_dist_key_dx &key) const {
		return dcp.computeDifferentialOperator(key, o1);
	}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		// for all NN of key
		for (int j = 0; j < dcp.getNumNN(key); j++) {
			auto coeff_dc = dcp.getCoeffNN(key, j);
			auto k = dcp.getIndexNN(key, j);

			auto k_coeff = coeff_dc * coeff * dcp.getEpsilonInvPrefactor(key);
			o1.template value_nz<Sys_eqs>(p_map, k, cols, k_coeff, comp);

			auto kk_coeff = dcp.getSign() * coeff_dc * coeff * dcp.getEpsilonInvPrefactor(key);
			o1.template value_nz<Sys_eqs>(p_map, key, cols, kk_coeff, comp);
		}
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	vtype &getVector() {
		return o1.getVector();
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype &getVector() const {
		return o1.getVector();
	}
};

template<typename exp1, typename DCPSE_type>
class vector_dist_expression_op<exp1, DCPSE_type, VECT_DCPSE_V> {
	//! expression 1
	const exp1 o1;

	DCPSE_type (&dcp)[DCPSE_type::vtype::dims];

	static const int dims = DCPSE_type::vtype::dims;
	typedef typename DCPSE_type::vtype::stype stype;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	typedef typename exp1::vtype vtype;

	//! Costruct a subtraction expression out of two expressions
	inline vector_dist_expression_op(const exp1 &o1, DCPSE_type (&dcp)[DCPSE_type::vtype::dims])
			: o1(o1), dcp(dcp) {}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const {
		o1.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=VectorS<dims, stype> >
	inline r_type value(const vect_dist_key_dx &key) const {
		VectorS<dims, stype> v_grad;

		for (int i = 0; i < dims; i++) {
			v_grad.get(i) = dcp[i].computeDifferentialOperator(key, o1);
		}

		return v_grad;
	}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		for (int i = 0; i < DCPSE_type::vtype::dims; i++) {
			// for all NN of key
			for (int j = 0; j < dcp[i].getNumNN(key); j++) {
				auto coeff_dc = dcp[i].getCoeffNN(key, j);
				auto k = dcp[i].getIndexNN(key, j);

				auto k_coeff = coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, k, cols, k_coeff, comp);

				auto kk_coeff = dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, key, cols, kk_coeff, comp);


				//cols[p_map. template getProp<0>(k)*Sys_eqs::nvar + comp] += coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);

				//cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
			}
		}
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	vtype &getVector() {
		return o1.getVector();
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype &getVector() const {
		return o1.getVector();
	}

};


template<typename exp1, typename DCPSE_type>
class vector_dist_expression_op<exp1, DCPSE_type, VECT_DCPSE_V_CURL2D> {
	//! expression 1
	const exp1 o1;

	DCPSE_type (&dcp)[DCPSE_type::vtype::dims];

	static const int dims = DCPSE_type::vtype::dims;
	typedef typename DCPSE_type::vtype::stype stype;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	typedef typename exp1::vtype vtype;

	//! Costruct a subtraction expression out of two expressions
	inline vector_dist_expression_op(const exp1 &o1, DCPSE_type (&dcp)[DCPSE_type::vtype::dims])
			: o1(o1), dcp(dcp) {}

	/*! \brief This function must be called before value
	 *
	 * it initialize the expression if needed
	 *
	 */
	inline void init() const {
		o1.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 * \return the result of the expression
	 *
	 */
	template<typename r_type=VectorS<dims, stype> >
	inline r_type value(const vect_dist_key_dx &key) const {
		VectorS<dims, stype> v_grad;
		v_grad.get(0) = dcp[0].computeDifferentialOperator(key, o1);
		v_grad.get(1) = -dcp[1].computeDifferentialOperator(key, o1);

		return v_grad;
	}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		for (int i = 0; i < DCPSE_type::vtype::dims; i++) {
			// for all NN of key
			for (int j = 0; j < dcp[i].getNumNN(key); j++) {
				auto coeff_dc = dcp[i].getCoeffNN(key, j);
				auto k = dcp[i].getIndexNN(key, j);

				auto k_coeff = coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, k, cols, k_coeff, comp);

				auto kk_coeff = dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, key, cols, kk_coeff, comp);


				//cols[p_map. template getProp<0>(k)*Sys_eqs::nvar + comp] += coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);

				//cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
			}
		}
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	vtype &getVector() {
		return o1.getVector();
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype &getVector() const {
		return o1.getVector();
	}

};


template<>
class vector_dist_expression_op<void, void, VECT_COPY_N_TO_N> {
	mutable int i;

	//! expression 1
	openfpm::vector<aggregate<int>> &l1;
	openfpm::vector<aggregate<int>> &l2;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	inline vector_dist_expression_op(openfpm::vector<aggregate<int>> &l1, openfpm::vector<aggregate<int>> &l2)
			: l1(l1), l2(l2) {}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		if (l1.template get<0>(i) != key.getKey()) {
			std::cout << "ERROR" << std::endl;
		}

		cols[p_map.template getProp<0>(key) * Sys_eqs::nvar + comp] += coeff;
		std::cout << "L2: " << l2.template get<0>(i) << std::endl;
		cols[p_map.template getProp<0>(l2.template get<0>(i)) * Sys_eqs::nvar + comp] -= coeff;

		i++;
	}
};


template<>
class vector_dist_expression_op<void, void, VECT_COPY_1_TO_N> {
	mutable int i = 0;

	//! expression 1
	openfpm::vector<aggregate<int>> &l1;
	int l2_key;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	inline vector_dist_expression_op(openfpm::vector<aggregate<int>> &l1, int l2_key)
			: l1(l1), l2_key(l2_key) {}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		if (l1.template get<0>(i) != key.getKey()) {
			std::cout << "ERROR" << std::endl;
		}

		cols[p_map.template getProp<0>(key) * Sys_eqs::nvar + comp] += coeff;
		cols[p_map.template getProp<0>(l2_key) * Sys_eqs::nvar + comp] -= coeff;
		i++;
	}
};


template<typename exp1, typename DCPSE_type>
class vector_dist_expression_op<exp1, DCPSE_type, VECT_DCPSE_V_SUM> {
	//! expression 1
	const exp1 o1;

	DCPSE_type (&dcp)[DCPSE_type::vtype::dims];

	static const int dims = DCPSE_type::vtype::dims;
	typedef typename DCPSE_type::vtype::stype stype;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	typedef typename exp1::vtype vtype;

	inline vector_dist_expression_op(const exp1 &o1, DCPSE_type (&dcp)[DCPSE_type::vtype::dims])
			: o1(o1), dcp(dcp) {}

	inline void init() const {
		o1.init();
	}

	template<typename r_type= typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)))>::type>
	inline r_type value(const vect_dist_key_dx &key) const {

		typename std::remove_reference<decltype(o1.value(key))>::type v_lap;
		v_lap = 0.0;


		for (int i = 0; i < dims; i++) {
			v_lap += dcp[i].computeDifferentialOperator(key, o1);
		}

		return v_lap;
	}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		for (int i = 0; i < DCPSE_type::vtype::dims; i++) {
			// for all NN of key
			for (int j = 0; j < dcp[i].getNumNN(key); j++) {
				auto coeff_dc = dcp[i].getCoeffNN(key, j);
				auto k = dcp[i].getIndexNN(key, j);

				auto coeff_k = coeff_dc * coeff * dcp[i].getEpsilonInvPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, k, cols, coeff_k, comp);

				auto coeff_kk = dcp[i].getSign() * coeff_dc * coeff * dcp[i].getEpsilonInvPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, key, cols, coeff_kk, comp);
			}
		}
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	vtype &getVector() {
		return o1.getVector();
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype &getVector() const {
		return o1.getVector();
	}
};


template<typename exp1, typename DCPSE_type>
class vector_dist_expression_op<exp1, DCPSE_type, VECT_DCPSE_V_DIV> {
	//! expression 1
	const exp1 o1;

	DCPSE_type (&dcp)[DCPSE_type::vtype::dims];

	static const int dims = DCPSE_type::vtype::dims;
	typedef typename DCPSE_type::vtype::stype stype;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	typedef typename exp1::vtype vtype;

	inline vector_dist_expression_op(const exp1 &o1, DCPSE_type (&dcp)[DCPSE_type::vtype::dims])
			: o1(o1), dcp(dcp) {}

	inline void init() const {
		o1.init();
	}

	template<typename r_type= typename std::remove_reference<decltype(o1.value(vect_dist_key_dx(0)))>::type::coord_type>
	inline r_type value(const vect_dist_key_dx &key) const {
		//typedef typename std::remove_reference<decltype(o1.value(key))>::type::blabla blabla;

		typename std::remove_reference<decltype(o1.value(key))>::type::coord_type v_div;
		v_div = 0.0;

		for (int i = 0; i < dims; i++) {
			v_div += dcp[i].computeDifferentialOperator(key, o1, i);
		}

		return v_div;
	}

	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		for (int i = 0; i < DCPSE_type::vtype::dims; i++) {
			// for all NN of key
			for (int j = 0; j < dcp[i].getNumNN(key); j++) {
				auto coeff_dc = dcp[i].getCoeffNN(key, j);
				auto k = dcp[i].getIndexNN(key, j);

				auto k_coeff = coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, k, cols, k_coeff, comp);

				auto kk_coeff = dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				o1.template value_nz<Sys_eqs>(p_map, key, cols, kk_coeff, comp);


				//cols[p_map. template getProp<0>(k)*Sys_eqs::nvar + comp] += coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);

				//cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
			}
		}
	}

	/*! \brief Return the vector on which is acting
	 *
	 * It return the vector used in getVExpr, to get this object
	 *
	 * \return the vector
	 *
	 */
	vtype &getVector() {
		return o1.getVector();
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype &getVector() const {
		return o1.getVector();
	}
};

template<typename exp1, typename exp2_pr>
class vector_dist_expression_op<exp1, exp2_pr, VECT_DCPSE_V_DOT> {
	typedef typename std::tuple_element<1, exp2_pr>::type DCPSE_type;
	typedef typename std::tuple_element<0, exp2_pr>::type exp2;

	//! expression 1
	const exp1 o1;
	const exp2 o2;

	DCPSE_type (&dcp)[DCPSE_type::vtype::dims];

	static const int dims = DCPSE_type::vtype::dims;
	typedef typename DCPSE_type::vtype::stype stype;

public:

	typedef std::false_type is_ker;

	typedef std::false_type NN_type;

	//! The type of the internal vector
	typedef typename first_or_second<has_vtype<exp1>::value, exp1, exp2>::vtype vtype;
	//typedef typename exp2::vtype vtype;


	inline vector_dist_expression_op(const exp1 &o1, const exp2 &o2, DCPSE_type (&dcp)[DCPSE_type::vtype::dims])
			: o1(o1), o2(o2), dcp(dcp) {}

	inline void init() const {
		o1.init();
		o2.init();
	}

	//template<typename r_type=VectorS<dims,stype> > inline r_type value(const vect_dist_key_dx & key) const
	template<typename r_type= typename std::remove_reference<decltype(o2.value(vect_dist_key_dx(0)))>::type>
	inline r_type value(const vect_dist_key_dx &key) const {
		//typedef typename std::remove_reference<decltype(o1.value(key))>::type::blabla blabla;
		typename std::remove_reference<decltype(o2.value(key))>::type adv;
		adv = 0.0;
		for (int i = 0; i < dims; i++) {
			adv += o1.value(key)[i] * dcp[i].computeDifferentialOperator(key, o2);

		}
		return adv;
	}


	template<typename Sys_eqs, typename pmap_type, typename unordered_map_type, typename coeff_type>
	inline void value_nz(pmap_type &p_map, const vect_dist_key_dx &key, unordered_map_type &cols, coeff_type &coeff,
						 unsigned int comp) const {
		for (int i = 0; i < DCPSE_type::vtype::dims; i++) {
			// for all NN of key
			for (int j = 0; j < dcp[i].getNumNN(key); j++) {
				auto coeff_dc = dcp[i].getCoeffNN(key, j);
				auto k = dcp[i].getIndexNN(key, j);

				auto k_coeff = o1.value(k) * coeff_dc * coeff * dcp[i].getEpsilonInvPrefactor(key);
				o2.template value_nz<Sys_eqs>(p_map, k, cols, k_coeff, comp);

				auto kk_coeff =
						o1.value(key) * dcp[i].getSign() * coeff_dc * coeff * dcp[i].getEpsilonInvPrefactor(key);
				o2.template value_nz<Sys_eqs>(p_map, key, cols, kk_coeff, comp);

				//cols[p_map. template getProp<0>(k)*Sys_eqs::nvar + comp] += o1.value(key)[i] * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
				//cols[p_map. template getProp<0>(key)*Sys_eqs::nvar + comp] += o1.value(key)[i] * dcp[i].getSign() * coeff_dc * coeff / dcp[i].getEpsilonPrefactor(key);
			}
		}
	}

	/*! \brief Return the vector on which is acting
 *
 * It return the vector used in getVExpr, to get this object
 *
 * \return the vector
 *
 */
	vtype &getVector() {
		return first_or_second<has_vtype<exp1>::value, exp1, exp2>::getVector(o1, o2);
	}

	/*! \brief Return the vector on which is acting
	*
	* It return the vector used in getVExpr, to get this object
	*
	* \return the vector
	*
	*/
	const vtype &getVector() const {
		return first_or_second<has_vtype<exp1>::value, exp1, exp2>::getVector(o1, o2);
	}

};

	/*
	template<typename operand_type>
	class Derivative_x_node
	{
		operand_type arg;

	public:
		typedef int it_is_a_node;
		Derivative_x_node(operand_type &arg)
		:arg(arg)
		{}
	};
	*/


/*! \brief Class for Creating the DCPSE Operator Dx and objects and computes DCPSE Kernels.
 *
 *
 * \param parts particle set
 * \param ord order of convergence of the operator
 * \param rCut Argument for cell list construction
 * \param support_options default: RADIUS to select all particles inside rCut
 *
 * \return Operator Dx which is a function on Vector_dist_Expressions
 *
 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_x_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dx and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dx which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_x_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dx and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dx which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_x_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernelNN(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernelNN<prp>(particles, k);

	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}

};

/*! \brief Class for Creating the DCPSE Operator Dy and objects and computes DCPSE Kernels.
 *
 *
 * \param parts particle set
 * \param ord order of convergence of the operator
 * \param rCut Argument for cell list construction
 * \param support_options default: RADIUS to select all particles inside rCut
 *
 * \return Operator Dy which is a function on Vector_dist_Expressions
 *
 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_y_T {

	void *dcpse;

public:

	/*! \brief Constructor for Creating the DCPSE Operator Dy and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dy which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_y_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(1) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dy and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dy which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_y_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(1) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>

	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse2 = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse2->template DrawKernel<prp>(particles, k);

	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

/*! \brief Class for Creating the DCPSE Operator Dz and objects and computes DCPSE Kernels.
 *
 *
 * \param parts particle set
 * \param ord order of convergence of the operator
 * \param rCut Argument for cell list construction
 * \param support_options default: RADIUS to select all particles inside rCut
 *
 * \return Operator Dz which is a function on Vector_dist_Expressions
 *
 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_z_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dz which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_z_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(2) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dz which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_z_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(2) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

};

/*! \brief Class for Creating the Gradient Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Grad which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Gradient_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the Gradient Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Grad which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Gradient_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<particles_type::dims, VerletList_type, particles_type> DCPSE_type;

		dcpse = new unsigned char[particles_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			Point<particles_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 1;

			new(&dcpse_ptr[i]) Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
		}
	}

	/*! \brief Constructor for Creating the Gradient Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Grad which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Gradient_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> DCPSE_type;

		dcpse = new unsigned char[vector_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *dcpse_ptr = (Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *) dcpse;

		for (int i = 0; i < vector_type::dims; i++) {
			Point<vector_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 1;

			new(&dcpse_ptr[i]) Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
				particlesSupport,
				particlesDomain,
				verletList,
				p,
				ord,
				rCut,
				opt
			);
		}
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		for (int i = 0; i < particles_type::dims; i++) {
			delete &(((Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse)[i]);
		}
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE_V>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE_V>(arg,
																				 *(dcpse_type(*)[operand_type::vtype::dims]) dcpse);
	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].template DrawKernel<prp>(particles, i, k);
		}

	}

	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].initializeUpdate(particles);
		}

	}


};
/*! \brief Class for Creating the DCPSE 2D Curl Operator
	 *
	 *  Creates object which work in 2 dimension and computes DCPSE Kernels for Dx_i as required.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Curl2D_T {

	void *dcpse;
public:
	/*! \brief Constructor for Creating the DCPSE 2D Curl Operator
	 *
	 *  Creates object which work in 2 dimension and computes DCPSE Kernels for Dx_i as required.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Curl2D_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<particles_type::dims, VerletList_type, particles_type> DCPSE_type;

		dcpse = new unsigned char[particles_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		Point<particles_type::dims, unsigned int> p;

		p.zero();
		p.get(1) = 1;
		new(dcpse_ptr) Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);

		p.zero();
		p.get(0) = 1;
		new(dcpse_ptr+1) Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);

	}

	/*! \brief Constructor for Creating the DCPSE 2D Curl Operator
	 *
	 *  Creates object which work in 2 dimension and computes DCPSE Kernels for Dx_i as required.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Curl2D_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> DCPSE_type;

		dcpse = new unsigned char[vector_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *dcpse_ptr = (Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *) dcpse;
		Point<vector_type::dims, unsigned int> p;

		p.zero();
		p.get(1) = 1;
		new(dcpse_ptr) Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);

		p.zero();
		p.get(0) = 1;
		new(dcpse_ptr+1) Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);

	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE_V_CURL2D>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE_V_CURL2D>(arg,
																						*(dcpse_type(*)[operand_type::vtype::dims]) dcpse);
	}
};
/*! \brief Class for Creating the DCPSE Laplacian Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Laplacian_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Laplacian Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Laplacian_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<particles_type::dims, VerletList_type, particles_type> DCPSE_type;
		dcpse = new unsigned char[particles_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			Point<particles_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 2;

			new(&dcpse_ptr[i]) Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
		}
	}

	/*! \brief Constructor for Creating the DCPSE Laplacian Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Laplacian_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> DCPSE_type;
		dcpse = new unsigned char[vector_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *dcpse_ptr = (Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *) dcpse;

		for (int i = 0; i < vector_type::dims; i++) {
			Point<vector_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 2;

			new(&dcpse_ptr[i]) Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
				particlesSupport,
				particlesDomain,
				verletList,
				p,
				ord,
				rCut,
				opt
			);
		}
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE_V_SUM>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE_V_SUM>(arg,
																					 *(dcpse_type(*)[operand_type::vtype::dims]) dcpse);
	}


	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].checkMomenta(particles);
		}

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].template DrawKernel<prp>(particles, k);
		}

	}
	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].initializeUpdate(particles);
		}

	}


};

/*! \brief Class for Creating the DCPSE Divergence Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions. Computes Divergence of Vectors
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Divergence_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Divergence Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions. Computes Divergence of Vectors
	 *
	 */
	template<typename particles_type>
	Divergence_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<particles_type::dims, VerletList_type, particles_type> DCPSE_type;

		dcpse = new unsigned char[particles_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			Point<particles_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 1;

			new(&dcpse_ptr[i]) Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
		}
	}

	/*! \brief Constructor for Creating the DCPSE Divergence Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions. Computes Divergence of Vectors
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Divergence_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> DCPSE_type;

		dcpse = new unsigned char[vector_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *dcpse_ptr = (Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *) dcpse;

		for (int i = 0; i < vector_type::dims; i++) {
			Point<vector_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 1;

			new(&dcpse_ptr[i]) Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
				particlesSupport,
				particlesDomain,
				verletList,
				p,
				ord,
				rCut,
				opt
			);
		}
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE_V_DIV>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE_V_DIV>(arg,
																					 *(dcpse_type(*)[operand_type::vtype::dims]) dcpse);
	}

	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].initializeUpdate(particles);
		}

	}

};

/*! \brief Class for Creating the DCPSE Advection Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions. Computes Advection of Vectors Adv(v,u) = v.Grad(u)
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Advection_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Advection Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions. Computes Advection of Vectors Adv(v,u) = v.Grad(u)
	 *
	 */
	template<typename particles_type>
	Advection_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<particles_type::dims, VerletList_type, particles_type> DCPSE_type;

		dcpse = new unsigned char[particles_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			Point<particles_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 1;

			new(&dcpse_ptr[i]) Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
		}


	}

	/*! \brief Constructor for Creating the DCPSE Advection Operator
	 *
	 *  Creates object which work on any dimension and computes DCPSE Kernels for Dx_i in each dimension.
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator which is a function on Vector_dist_Expressions. Computes Advection of Vectors Adv(v,u) = v.Grad(u)
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Advection_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		typedef Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> DCPSE_type;

		dcpse = new unsigned char[vector_type::dims * sizeof(DCPSE_type)];

		Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *dcpse_ptr = (Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2> *) dcpse;

		for (int i = 0; i < vector_type::dims; i++) {
			Point<vector_type::dims, unsigned int> p;
			p.zero();
			p.get(i) = 1;

			new(&dcpse_ptr[i]) Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
				particlesSupport,
				particlesDomain,
				verletList,
				p,
				ord,
				rCut,
				opt
			);
		}


	}

	template<typename operand_type1, typename operand_type2>
	vector_dist_expression_op<operand_type1, std::pair<operand_type2, Dcpse_type<operand_type2::vtype::dims, VerletList_type, typename operand_type2::vtype>>, VECT_DCPSE_V_DOT>
	operator()(operand_type1 arg, operand_type2 arg2) {
		typedef Dcpse_type<operand_type2::vtype::dims, VerletList_type, typename operand_type2::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type1, std::pair<operand_type2, dcpse_type>, VECT_DCPSE_V_DOT>(arg,
																												arg2,
																												*(dcpse_type(*)[operand_type2::vtype::dims]) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].checkMomenta(particles);
		}

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].template DrawKernel<prp>(particles, i, k);
		}

	}

	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		for (int i = 0; i < particles_type::dims; i++) {
			dcpse_ptr[i].initializeUpdate(particles);
		}

	}


};

/*! \brief Class for Creating the DCPSE Operator Dxy and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxy which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xy_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dxy and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxy which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_xy_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;
		p.get(1) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dxy and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxy which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xy_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;
		p.get(1) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		dcpse_ptr[0].template DrawKernel<prp>(particles, k);

	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};
/*! \brief Class for Creating the DCPSE Operator Dyz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dyz which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_yz_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dyz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dyz which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_yz_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(1) = 1;
		p.get(2) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dyz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dyz which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_yz_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(1) = 1;
		p.get(2) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		dcpse_ptr[0].template DrawKernel<prp>(particles, k);

	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};
/*! \brief Class for Creating the DCPSE Operator Dxz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxz which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xz_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dxz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxz which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_xz_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;
		p.get(2) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dxz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxz which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xz_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;
		p.get(2) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		Dcpse_type<particles_type::dims, VerletList_type, particles_type> *dcpse_ptr = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;

		dcpse_ptr[0].template DrawKernel<prp>(particles, k);

	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

/*! \brief Constructor for Creating the DCPSE Operator Dxx and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxx which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xx_T {

	void *dcpse;

public:
	/*! \brief Class for Creating the DCPSE Operator Dxx and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxx which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_xx_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 2;
		p.get(1) = 0;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Class for Creating the DCPSE Operator Dxx and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dxx which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xx_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 2;
		p.get(1) = 0;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

/*! \brief Class for Creating the DCPSE Operator Dyy and objects and computes DCPSE Kernels.
 *
 *
 * \param parts particle set
 * \param ord order of convergence of the operator
 * \param rCut Argument for cell list construction
 * \param support_options default: RADIUS to select all particles inside rCut
 *
 * \return Operator Dyy which is a function on Vector_dist_Expressions
 *
 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_yy_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dyy and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dyy which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_yy_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 0;
		p.get(1) = 2;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dyy and objects and computes DCPSE Kernels.
	 *
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dyy which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_yy_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 0;
		p.get(1) = 2;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};
/*! \brief Class for Creating the DCPSE Operator Dzz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dzz which is a function on Vector_dist_Expressions
	 *
	 */
template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_zz_T {

	void *dcpse;

public:
	/*! \brief Constructor for Creating the DCPSE Operator Dzz and objects and computes DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dzz which is a function on Vector_dist_Expressions
	 *
	 */
	template<typename particles_type>
	Derivative_zz_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(2) = 2;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	/*! \brief Constructor for Creating the DCPSE Operator Dzz and objects and computes DCPSE Kernels.
	 *
	 *
	 *
	 * \param parts set of support particles
	 * \param parts set of domain particles
	 * \param ord order of convergence of the operator
	 * \param rCut Argument for cell list construction
	 * \param support_options default: RADIUS to select all particles inside rCut
	 *
	 * \return Operator Dzz which is a function on Vector_dist_Expressions
	 *
	 */
	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_zz_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(2) = 2;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename particles_type>
	void deallocate(particles_type &parts) {
		delete (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xxx_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_xxx_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 3;
		p.get(1) = 0;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xxx_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 3;
		p.get(1) = 0;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>

	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xxy_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_xxy_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 2;
		p.get(1) = 1;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xxy_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 2;
		p.get(1) = 1;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>

	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_yyx_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_yyx_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;
		p.get(1) = 2;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_yyx_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 1;
		p.get(1) = 2;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>

	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_yyy_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_yyy_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 0;
		p.get(1) = 3;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_yyy_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 0;
		p.get(1) = 3;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>

	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};


template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xxxx_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_xxxx_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 4;
		p.get(1) = 0;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xxxx_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 4;
		p.get(1) = 0;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};


template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_yyyy_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_yyyy_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 0;
		p.get(1) = 4;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_yyyy_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 0;
		p.get(1) = 4;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_xxyy_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_xxyy_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<particles_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 2;
		p.get(1) = 2;

		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_xxyy_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		support_options opt = support_options::RADIUS
	) {
		Point<vector_type::dims, unsigned int> p;
		p.zero();
		p.get(0) = 2;
		p.get(1) = 2;

		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>

	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};


template<typename VerletList_type, template<unsigned int, typename, typename...> class Dcpse_type = Dcpse>
class Derivative_G_T {

	void *dcpse;

public:

	template<typename particles_type>
	Derivative_G_T(
		particles_type &parts,
		VerletList_type& verletList,
		unsigned int ord,
		typename particles_type::stype rCut,
		const Point<particles_type::dims, unsigned int> &p,
		support_options opt = support_options::RADIUS
	) {
		dcpse = new Dcpse_type<particles_type::dims, VerletList_type, particles_type>(parts, verletList, p, ord, rCut, opt);
	}

	template<
		typename vector_type,
		typename vector_type2,
		std::enable_if_t<std::is_same<int, typename vector_type2::yes_i_am_vector_subset>::value>* = nullptr>
	Derivative_G_T(
		vector_type &particlesSupport,
		vector_type2 &particlesDomain,
		VerletList_type& verletList,
		unsigned int ord,
		typename vector_type::stype rCut,
		const Point<vector_type::dims, unsigned int> &p,
		support_options opt = support_options::RADIUS
	) {
		dcpse = new Dcpse_type<vector_type::dims, VerletList_type, vector_type, vector_type2>(
			particlesSupport,
			particlesDomain,
			verletList,
			p,
			ord,
			rCut,
			opt
		);
	}

	template<typename operand_type>
	vector_dist_expression_op<operand_type, Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype>, VECT_DCPSE>
	operator()(operand_type arg) {
		typedef Dcpse_type<operand_type::vtype::dims, VerletList_type, typename operand_type::vtype> dcpse_type;

		return vector_dist_expression_op<operand_type, dcpse_type, VECT_DCPSE>(arg, *(dcpse_type *) dcpse);
	}

	template<typename particles_type>
	void checkMomenta(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->checkMomenta(particles);

	}

	template<unsigned int prp, typename particles_type>
	void DrawKernel(particles_type &particles, int k) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->template DrawKernel<prp>(particles, k);

	}

	/*! \brief Method for Saving the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be saved.
	 */
	template<typename particles_type>
	void save(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->save(file);
	}
	/*! \brief Method for Loading the DCPSE Operator.
	 *
	 * \param parts particle set
	 * \param file name for data to be loaded from.
	 */
	template<typename particles_type>
	void load(particles_type &particles, const std::string &file) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->load(file);
	}
	/*! \brief Method for Updating the DCPSE Operator by recomputing DCPSE Kernels.
	 *
	 *
	 * \param parts particle set
	 */
	template<typename particles_type>
	void update(particles_type &particles) {
		auto dcpse_temp = (Dcpse_type<particles_type::dims, VerletList_type, particles_type> *) dcpse;
		dcpse_temp->initializeUpdate(particles);

	}
};

//typedef PPInterpolation_T<Dcpse> PPInterpolation;
template <typename VerletList_type> using Derivative_x = Derivative_x_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_x = Derivative_x_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_y = Derivative_y_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_z = Derivative_z_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Gradient = Gradient_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Curl2D = Curl2D_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Laplacian = Laplacian_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Divergence = Divergence_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Advection = Advection_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xy = Derivative_xy_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_yz = Derivative_yz_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xz = Derivative_xz_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xx = Derivative_xx_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_yy = Derivative_yy_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_zz = Derivative_zz_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xxx = Derivative_xxx_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xxy = Derivative_xxy_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_yyx = Derivative_yyx_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_yyy = Derivative_yyy_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xxxx = Derivative_xxxx_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_yyyy = Derivative_yyyy_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_xxyy = Derivative_xxyy_T<VerletList_type, Dcpse>;
template <typename VerletList_type> using Derivative_G = Derivative_G_T<VerletList_type, Dcpse>;


#if defined(__NVCC__)
template <typename VerletList_type> using Derivative_x_gpu = Derivative_x_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_y_gpu = Derivative_y_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_z_gpu = Derivative_z_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Gradient_gpu = Gradient_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Curl2D_gpu = Curl2D_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Laplacian_gpu = Laplacian_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Divergence_gpu = Divergence_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Advection_gpu = Advection_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_xy_gpu = Derivative_xy_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_yz_gpu = Derivative_yz_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_xz_gpu = Derivative_xz_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_xx_gpu = Derivative_xx_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_yy_gpu = Derivative_yy_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_zz_gpu = Derivative_zz_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_xxx_gpu = Derivative_xxx_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_xxy_gpu = Derivative_xxy_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_yyx_gpu = Derivative_yyx_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_yyy_gpu = Derivative_yyy_T<VerletList_type, Dcpse_gpu>;
template <typename VerletList_type> using Derivative_G_gpu = Derivative_G_T<VerletList_type, Dcpse_gpu>;
#endif

#endif /*EIGEN */
#endif /* DCPSE_OP_HPP_ */
