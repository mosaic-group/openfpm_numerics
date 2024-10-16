/*
 * Vector_petsc.hpp
 *
 *  Created on: Apr 29, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_PETSC_HPP_
#define OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_PETSC_HPP_

#include "Vector/map_vector.hpp"
#include "Vector/vector_def.hpp"
#include <boost/mpl/int.hpp>
#include <petscvec.h>
#include "util/petsc_util.hpp"
#include <unordered_map>
#include "VTKWriter/VTKWriter.hpp"
#include "CSVWriter/CSVWriter.hpp"

#define PETSC_RVAL 2

/*! \brief It store one row value of a vector
 *
 * Given a row, store a value
 *
 *
 */
template<typename T>
class rval<T,PETSC_RVAL>
{
public:

	//! boost fusion that store the point
	typedef boost::fusion::vector<PetscInt,T> type;

	//! structure that store the data of the point
	type data;

	//! Property id for the row
	static const unsigned int row = 0;

	//! Property id for the value
	static const unsigned int value = 1;

	//! This object has 2 properties
	static const unsigned int max_prop = 2;

	/*! \brief Get the row
	 *
	 * \return the row
	 *
	 */
	long int & rw()
	{
		return boost::fusion::at_c<row>(data);
	}

	/*! \brief Get the value
	 *
	 * \return the value
	 *
	 */
	T & val()
	{
		return boost::fusion::at_c<value>(data);
	}

	/*! \brief Default constructor
	 *
	 */
	rval()	{}

	/*! \brief Constructor from row, column and value
	 *
	 * \param i row
	 * \param val value
	 *
	 */
	rval(long int i, T val)
	{
		rw() = i;
		val() = val;
	}

	/*! \brief Indicate that the structure has no pointer
	 *
	 * \return true
	 *
	 */
	static inline bool noPointers()
	{
		return true;
	}
};

constexpr unsigned int row_id = 0;
constexpr unsigned int val_id = 1;


/*! \brief PETSC vector for linear algebra
 *
 * This vector wrap the PETSC vector for solving linear systems
 *
 */
template<typename T>
class Vector<T,PETSC_BASE>
{
	//! Number of row the petsc vector has
	size_t n_row;

	//! Number of local rows
	size_t n_row_local;

	//! Indicate if v has been allocated
	mutable bool v_created = false;

	//! Mutable vector
	mutable Vec v;

	//! Mutable row value vector
	mutable openfpm::vector<rval<PetscScalar,PETSC_RVAL>,HeapMemory, memory_traits_inte > row_val;

	//! Global to local map
	mutable std::unordered_map<size_t,size_t> map;

	//! invalid
	T invalid;

	/*! \brief Set the Eigen internal vector
	 *
	 *
	 */
	void setPetsc() const
	{
		if (v_created == false)
		{PETSC_SAFE_CALL(VecSetType(v,VECMPI));}

		// set the vector

		if (row_val.size() != 0)
		{PETSC_SAFE_CALL(VecSetValues(v,row_val.size(),&row_val.template get<row_id>(0),&row_val.template get<val_id>(0),INSERT_VALUES))}

		PETSC_SAFE_CALL(VecAssemblyBegin(v));
		PETSC_SAFE_CALL(VecAssemblyEnd(v));

		v_created = true;
	}

public:

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(const Vector<T,PETSC_BASE> & v)
	:Vector()
	{
		this->operator=(v);
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector(Vector<T,PETSC_BASE> && v)
	:Vector()
	{
		this->operator=(v);
	}

	/*! \brief Destroy the vector
	 *
	 *
	 */
	~Vector()
	{
		if (is_openfpm_init() == true)
		{PETSC_SAFE_CALL(VecDestroy(&v));}
	}

	/*! \brief Create a vector with n elements
	 *
	 * \param n global number of elements in the vector
	 * \param n_row_local number
	 *
	 */
	Vector(size_t n, size_t n_row_local)
	:n_row_local(n_row_local),v(NULL),invalid(0)
	{
		// Create the vector
		PETSC_SAFE_CALL(VecCreate(PETSC_COMM_WORLD,&v));

		resize(n,n_row_local);
	}

	/*! \brief Create a vector with 0 elements
	 *
	 */
	Vector()
	:n_row(0),n_row_local(0),invalid(0)
	{
		// Create the vector
		PETSC_SAFE_CALL(VecCreate(PETSC_COMM_WORLD,&v));
	}

	/*! \brief Resize the Vector
	 *
	 * \param row numbers of row
	 * \param l_row number of local row
	 *
	 */
	void resize(size_t row, size_t l_row)
	{
		n_row = row;
		n_row_local = l_row;

		PETSC_SAFE_CALL(VecSetSizes(v,n_row_local,n_row));
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 * \param val value
	 *
	 */
	void insert(size_t i, T val)
	{
		row_val.add();

		// Map
		map[i] = row_val.size()-1;

		row_val.last().template get<row_id>() = i;
		row_val.last().template get<val_id>() = val;
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	inline PetscScalar & insert(size_t i)
	{
		row_val.add();

		// Map
		map[i] = row_val.size()-1;

		row_val.last().template get<row_id>() = i;
		return row_val.last().template get<val_id>();
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	inline const PetscScalar & insert(size_t i) const
	{
		row_val.add();

		// Map
		map[i] = row_val.size()-1;

		row_val.last().template get<row_id>() = i;
		return row_val.last().template get<val_id>();
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \warning The element must exist
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	const PetscScalar & operator()(size_t i) const
	{
		// Search if exist

		std::unordered_map<size_t,size_t>::iterator it = map.find(i);

		if ( it != map.end() )
			return row_val.template get<val_id>(it->second);

		return insert(i);
	}

	/*! \brief Return a reference to the vector element
	 *
	 * \warning The element must exist
	 *
	 * \param i element
	 *
	 * \return reference to the element vector
	 *
	 */
	PetscScalar & operator()(size_t i)
	{
		// Search if exist

		std::unordered_map<size_t,size_t>::iterator it = map.find(i);

		if ( it != map.end() )
			return row_val.template get<val_id>(it->second);

		return insert(i);
	}

	/*! \brief Get the PETSC Vector object
	 *
	 * \return the PETSC Vector
	 *
	 */
	const Vec & getVec() const
	{
		setPetsc();

		return v;
	}

	/*! \brief Get the PETSC Vector object
	 *
	 * \return the PETSC Vector
	 *
	 */
	Vec & getVec()
	{
		setPetsc();

		return v;
	}

	/*! \brief Update the Vector with the PETSC object
	 *
	 */
	void update()
	{
		PetscInt n_row;
		PetscInt n_row_local;

		// Get the size of the vector from PETSC
		VecGetSize(v,&n_row);
		VecGetLocalSize(v,&n_row_local);

		this->n_row = n_row;
		this->n_row_local = n_row_local;

		row_val.resize(n_row_local);

		//

		PetscInt low;
		PetscInt high;

		VecGetOwnershipRange(v,&low,&high);

		// Fill the index and construct the map

		size_t k = 0;
		for (size_t i = low ; i < (size_t)high ; i++)
		{
			row_val.template get<row_id>(k) = i;
			map[i] = k;
			k++;
		}

		PETSC_SAFE_CALL(VecGetValues(v,row_val.size(),&row_val.template get<row_id>(0),&row_val.template get<val_id>(0)))
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector<T,PETSC_BASE> & operator=(const Vector<T,PETSC_BASE> & v)
	{
		map = v.map;
		row_val = v.row_val;

		return *this;
	}

	/*! \brief Copy the vector
	 *
	 * \param v vector to copy
	 *
	 */
	Vector<T,PETSC_BASE> & operator=(Vector<T,PETSC_BASE> && v)
	{
		map.swap(v.map);
		row_val.swap(v.row_val);

		return *this;
	}

	/*! \brief Set to zero all the entries
	 *
	 *
	 */
	void setZero()
	{
		if (v_created == false)
		{PETSC_SAFE_CALL(VecSetType(v,VECMPI));}

		v_created = true;
	}

	/* Write vector on vtk
	 *
	 * \param out file to write into
	 *
	 */
	bool write(std::string out, size_t opt = VTK_WRITER)
	{
		Vcluster<> & v_cl = create_vcluster();

		openfpm::vector<Point<2, double>> row_col;
		openfpm::vector<aggregate<double>> values;

		row_col.resize(n_row_local);
		values.resize(n_row_local);

		int i = 0;
		for (auto it = map.begin() ; it != map.end() ; it++, i++)
		{
			row_col.template get<0>(i)[1] = it->first;
			row_col.template get<0>(i)[0] = 0.0;

			values.template get<0>(i) = row_val.template get<1>(it->second);
		}

		if (opt == VTK_WRITER)
		{
			auto ft = file_type::ASCII;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<2, double>>,
									   openfpm::vector<aggregate<double>>>,
									   VECTOR_POINTS> vtk_writer;

			vtk_writer.add(row_col,values,row_col.size());

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtk"));

			openfpm::vector<std::string> prp_names;
			prp_names.add("value");

			// Write the VTK file
			return vtk_writer.write(output,prp_names,"vector","",ft);
		}
		else
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<2,double>>,
			          openfpm::vector<aggregate<double>>> csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output,row_col,values);
		}
	}
};


#endif /* OPENFPM_NUMERICS_SRC_VECTOR_VECTOR_EIGEN_HPP_ */

