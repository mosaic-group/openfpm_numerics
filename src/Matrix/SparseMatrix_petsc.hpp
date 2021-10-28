/*
 * SparseMatrix_petsc.hpp
 *
 *  Created on: Apr 26, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_PETSC_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_PETSC_HPP_

#include "util/petsc_util.hpp"
#include "Vector/map_vector.hpp"
#include <boost/mpl/int.hpp>
#include <petscmat.h>
#include "VTKWriter/VTKWriter.hpp"
#include "CSVWriter/CSVWriter.hpp"

#define PETSC_BASE 2

/*! \brief It store one non-zero element in the sparse matrix
 *
 * Given a row, and a column, store a value
 *
 *
 */
template<typename T>
class triplet<T,PETSC_BASE>
{
	//! Row of the sparse matrix
	PetscInt row_;

	//! Colum of the sparse matrix
	PetscInt col_;

	//! Value of the Matrix
	PetscScalar val_;

public:

	/*! \brief Return the row of the triplet
	 *
	 * \return the row index
	 *
	 */
	PetscInt & row()
	{
		return row_;
	}

	/*! \brief Return the colum of the triplet
	 *
	 * \return the colum index
	 *
	 */
	PetscInt & col()
	{
		return col_;
	}

	/*! \brief Return the value of the triplet
	 *
	 * \return the value
	 *
	 */
	PetscScalar & value()
	{
		return val_;
	}

	/*! \brief Constructor from row, colum and value
	 *
	 * \param i row
	 * \param j colum
	 * \param val value
	 *
	 */
	triplet(long int i, long int j, T val)
	{
		row_ = i;
		col_ = j;
		val_ = val;
	}

	// Default constructor
	triplet()
	:row_(0),col_(0),val_(0)
	{};
};

/*! \brief Sparse Matrix implementation, that map over Eigen
 *
 * \tparam T Type of the sparse Matrix store on each row,colums
 * \tparam id_t type of id
 * \tparam impl implementation
 *
 */
template<typename T, typename id_t>
class SparseMatrix<T,id_t, PETSC_BASE>
{
public:

	//! Triplet implementation id
	typedef boost::mpl::int_<PETSC_BASE> triplet_impl;

	//! Triplet type
	typedef triplet<T,PETSC_BASE> triplet_type;

private:

	//! Number of matrix row (global)
	size_t g_row;
	//! Number of matrix colums (global)
	size_t g_col;

	//! Number of matrix row (local)
	size_t l_row;
	//! Number of matrix colums (local)
	size_t l_col;

	//! starting row for this processor
	size_t start_row;

	//! indicate if the matrix has been created
	bool m_created = false;

	//! PETSC Matrix
	Mat mat;

	//! Triplets of the matrix
	openfpm::vector<triplet_type> trpl;


	//! temporary list of values
	mutable openfpm::vector<PetscScalar> vals;
	//! temporary list of colums
	mutable openfpm::vector<PetscInt> cols;
	//! PETSC d_nnz
	mutable openfpm::vector<PetscInt> d_nnz;
	//! PETSC o_nnz
	mutable openfpm::vector<PetscInt> o_nnz;

	/*! \brief Fill the petsc Matrix
	 *
	 *
	 */
	void fill_petsc()
	{
		d_nnz.resize(l_row);
		o_nnz.resize(l_row);

		d_nnz.fill(0);
		o_nnz.fill(0);

		// Here we explore every row to count how many non zero we have in the diagonal matrix part,
		// and the non diagonal matrix part, needed by MatMPIAIJSetPreallocation

		size_t i = 0;

		// Set the Matrix from triplet
		while (i < trpl.size())
		{
			PetscInt row = trpl.get(i).row();

			while(i < trpl.size() && row == trpl.get(i).row())
			{
				if ((size_t)trpl.get(i).col() >= start_row && (size_t)trpl.get(i).col() < start_row + l_row)
					d_nnz.get(row - start_row)++;
				else
					o_nnz.get(row - start_row)++;
				i++;
			}
		}

		PETSC_SAFE_CALL(MatMPIAIJSetPreallocation(mat,0,static_cast<const PetscInt*>(d_nnz.getPointer()),0,
														static_cast<const PetscInt*>(o_nnz.getPointer())));

		// Counter i is zero
		i = 0;

		// Set the Matrix from triplet
		while (i < trpl.size())
		{
			vals.clear();
			cols.clear();

			PetscInt row = trpl.get(i).row();

			while(i < trpl.size() && row == trpl.get(i).row())
			{
				vals.add(trpl.get(i).value());
				cols.add(trpl.get(i).col());
				i++;
			}
			PETSC_SAFE_CALL(MatSetValues(mat,1,&row,cols.size(),static_cast<const PetscInt*>(cols.getPointer()),
					                                            static_cast<const PetscScalar *>(vals.getPointer()),
					                                            INSERT_VALUES));
		}

        PETSC_SAFE_CALL(MatAssemblyBegin(mat,MAT_FINAL_ASSEMBLY));
		PETSC_SAFE_CALL(MatAssemblyEnd(mat,MAT_FINAL_ASSEMBLY));

		m_created = true;
	}

	/*! \brief Disable copy constructor
	 *
	 *
	 */
	SparseMatrix(const SparseMatrix<T,id_t, PETSC_BASE> & spm)	{};


public:

	/*! \brief Create an empty Matrix
	 *
	 * \param N1 number of row
	 * \param N2 number of colums
	 * \param N1_loc number of local row
	 *
	 */
	SparseMatrix(size_t N1, size_t N2, size_t n_row_local)
	:g_row(N1),g_col(N2),l_row(n_row_local),l_col(n_row_local)
	{
		PETSC_SAFE_CALL(MatCreate(PETSC_COMM_WORLD,&mat));
		//PETSC_SAFE_CALL(MatSetType(mat,MATMPIAIJ));
        PETSC_SAFE_CALL(MatSetFromOptions(mat));
		PETSC_SAFE_CALL(MatSetSizes(mat,n_row_local,n_row_local,N1,N2));

		Vcluster<> & v_cl = create_vcluster();

		openfpm::vector<size_t> vn_row_local;
		v_cl.allGather(l_row,vn_row_local);
		v_cl.execute();

		// Calculate the starting row for this processor

		start_row = 0;
		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
			start_row += vn_row_local.get(i);
	}

	/*! \brief Create an empty Matrix
	 *
	 */
	SparseMatrix()
	:g_row(0),g_col(0),l_row(0l),l_col(0),start_row(0)
	{
		PETSC_SAFE_CALL(MatCreate(PETSC_COMM_WORLD,&mat));
        //PETSC_SAFE_CALL(MatSetType(mat,MATMPIAIJ));
        PETSC_SAFE_CALL(MatSetFromOptions(mat));

	}

	~SparseMatrix()
	{
		// Destroy the matrix
		if (is_openfpm_init() == true)
		{PETSC_SAFE_CALL(MatDestroy(&mat));}
	}

	/*! \brief Get the Matrix triplets buffer
	 *
	 * It return a buffer that can be filled with triplets
	 *
	 * \return Petsc Matrix
	 *
	 */
	openfpm::vector<triplet_type> & getMatrixTriplets()
	{
		m_created = false;

		return this->trpl;
	}

	/*! \brief Get the Patsc Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	const Mat & getMat() const
	{
		if (m_created == false)
		{fill_petsc();}

		return mat;
	}

	/*! \brief Get the Petsc Matrix object
	 *
	 * \return the Petsc Matrix
	 *
	 */
	Mat & getMat()
	{
		if (m_created == false)
		{fill_petsc();}

		return mat;
	}

	/*! \brief Resize the Sparse Matrix
	 *
	 * \param row number for row
	 * \param col number of colums
	 * \param local number of row
	 * \param local number of colums
	 *
	 */
	void resize(size_t row, size_t col, size_t l_row, size_t l_col)
	{
		if ((g_row != 0 && g_row != row) || (g_col != 0 && g_col != col) ||
			(this->l_row != 0 && this->l_row != l_row) || (this->l_col != 0 && this->l_col != l_col))
		{
			std::cout << __FILE__ << ":" << __LINE__ << " error you are resizing a PETSC matrix " << std::endl;
		}

		if (g_row != 0 && g_col != 0)	{return;}

		this->g_row = row;
		this->g_col = col;

		this->l_row = l_row;
		this->l_col = l_col;

		PETSC_SAFE_CALL(MatSetSizes(mat,l_row,l_col,g_row,g_col));

		Vcluster<> & v_cl = create_vcluster();

		openfpm::vector<size_t> vn_row_local;
		v_cl.allGather(l_row,vn_row_local);
		v_cl.execute();

		// Calculate the starting row for this processor

		start_row = 0;
		for (size_t i = 0 ; i < v_cl.getProcessUnitID() ; i++)
			start_row += vn_row_local.get(i);
	}

	/*! \brief Get the row i and the colum j of the Matrix
	 *
	 * \warning it is slow, consider to get blocks of the matrix
	 *
	 * \return the value of the matrix at row i colum j
	 *
	 */
	T operator()(id_t i, id_t j)
	{
		T v;

		MatGetValues(mat,1,&i,1,&j,&v);

		return v;
	}

	/*! \brief Get the value from triplet
	 *
	 * \warning It is extremly slow because it do a full search across the triplets elements
	 *
	 * \param r row
	 * \param c colum
	 *
	 */
	T getValue(size_t r, size_t c)
	{
		for (size_t i = 0 ; i < trpl.size() ; i++)
		{
			if (r == (size_t)trpl.get(i).row() && c == (size_t)trpl.get(i).col())
				return trpl.get(i).value();
		}

		return 0;
	}

	/* Write matrix on vtk
	 *
	 * \param out file to write into
	 *
	 */
	bool write(std::string out,size_t opt = VTK_WRITER)
	{
		Vcluster<> & v_cl = create_vcluster();

		openfpm::vector<Point<2, double>> row_col;
		openfpm::vector<aggregate<double>> values;

		row_col.resize(trpl.size());
		values.resize(trpl.size());

		for (int i = 0 ; i < trpl.size() ; i++)
		{
			row_col.template get<0>(i)[1] = trpl.get(i).row();
			row_col.template get<0>(i)[0] = trpl.get(i).col();

			values.template get<0>(i) = trpl.get(i).value();
		}

		if (opt == VTK_WRITER)
		{
			auto ft = file_type::BINARY;

			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<2, double>>,
									   openfpm::vector<aggregate<double>>>,
									   VECTOR_POINTS> vtk_writer;

			vtk_writer.add(row_col,values,row_col.size());

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtk"));

			openfpm::vector<std::string> prp_names;
			prp_names.add("value");

			// Write the VTK file
			return vtk_writer.write(output,prp_names,"matrix","",ft);
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


#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_PETSC_HPP_ */
