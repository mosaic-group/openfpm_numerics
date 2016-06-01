/*
 * SparseMatrix_Eigen.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: i-bird
 */

#ifndef OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_
#define OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_

#include "Vector/map_vector.hpp"
#include <boost/mpl/int.hpp>
#include "VCluster.hpp"

#define EIGEN_TRIPLET 1

/*! \brief It store one non-zero element in the sparse matrix
 *
 * Given a row, and a column, store a value
 *
 *
 */
template<typename T>
class triplet<T,EIGEN_TRIPLET>
{
	Eigen::Triplet<T,long int> triplet_ei;

public:

	long int & row()
	{
		size_t ptr = (size_t)&triplet_ei.row();

		return (long int &)*(long int *)ptr;
	}

	long int & col()
	{
		size_t ptr = (size_t)&triplet_ei.col();

		return (long int &)*(long int *)ptr;
	}

	T & value()
	{
		size_t ptr = (size_t)&triplet_ei.value();

		return (T &)*(T *)ptr;
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
		row() = i;
		col() = j;
		value() = val;
	}

	// Default constructor
	triplet()	{};
};

/* ! \brief Sparse Matrix implementation, that map over Eigen
 *
 * \tparam T Type of the sparse Matrix store on each row,colums
 * \tparam id_t type of id
 * \tparam impl implementation
 *
 */
template<typename T, typename id_t>
class SparseMatrix<T,id_t,EIGEN_BASE>
{
public:

	//! Triplet implementation id
	typedef boost::mpl::int_<EIGEN_BASE> triplet_impl;

	//! Triplet type
	typedef triplet<T,EIGEN_BASE> triplet_type;

private:

	Eigen::SparseMatrix<T,0,id_t> mat;
	openfpm::vector<triplet_type> trpl;
	openfpm::vector<triplet_type> trpl_recv;

	/*! \brief Assemble the matrix
	 *
	 *
	 */
	void assemble()
	{
		Vcluster & vcl = create_vcluster();

		////// On Master and only here
		// we assemble the Matrix from the collected data
		if (vcl.getProcessingUnits() != 1)
		{
			collect();
			// only master assemble the Matrix
			if (vcl.getProcessUnitID() == 0)
				mat.setFromTriplets(trpl_recv.begin(),trpl_recv.end());
		}
		else
			mat.setFromTriplets(trpl.begin(),trpl.end());
	}

	/*! \brief Here we collect the full matrix on master
	 *
	 */
	void collect()
	{
		Vcluster & vcl = create_vcluster();

		trpl_recv.clear();

		// here we collect all the triplet in one array on the root node
		vcl.SGather(trpl,trpl_recv,0);
	}

public:



	/*! \brief Create an empty Matrix
	 *
	 * \param N1 number of row
	 * \param N2 number of colums
	 *
	 */
	SparseMatrix(size_t N1, size_t N2)
	:mat(N1,N2)
	{
	}


	/*! \brief Create an empty Matrix
	 *
	 */
	SparseMatrix()	{}

	/*! \brief Get the Matrix triplets bugger
	 *
	 * It return a buffer that can be filled with triplets
	 *
	 * \return triplet buffer
	 *
	 */
	openfpm::vector<triplet_type> & getMatrixTriplets()
	{
		return this->trpl;
	}

	/*! \brief Get the Eigen Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	const Eigen::SparseMatrix<T,0,id_t> & getMat() const
	{
		// Here we collect the information on master
		assemble();

		return mat;
	}

	/*! \brief Get the Eigen Matrix object
	 *
	 * \return the Eigen Matrix
	 *
	 */
	Eigen::SparseMatrix<T,0,id_t> & getMat()
	{
		assemble();

		return mat;
	}

	/*! \brief Resize the Sparse Matrix
	 *
	 * \param row number for row
	 * \param col number of colums
	 *
	 */
	void resize(size_t row, size_t col,size_t l_row, size_t l_col)
	{
		mat.resize(row,col);
	}

	/*! \brief Get the row i and the colum j of the Matrix
	 *
	 * \return the value of the matrix at row i colum j
	 *
	 */
	T operator()(id_t i, id_t j)
	{
		return mat.coeffRef(i,j);
	}

	/*! \brief Save this object into file
	 *
	 * \param file filename
	 *
	 * \return true if succed
	 *
	 */
	bool save(const std::string & file) const
	{
		std::vector<size_t> pap_prp;

		Packer<decltype(trpl),HeapMemory>::packRequest(trpl,pap_prp);

		// Calculate how much preallocated memory we need to pack all the objects
		size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

		// allocate the memory
		HeapMemory pmem;
		pmem.allocate(req);
		ExtPreAlloc<HeapMemory> mem(pap_prp,pmem);

		//Packing

		Pack_stat sts;
		Packer<decltype(trpl),HeapMemory>::pack(mem,trpl,sts);

		// Save into a binary file
	    std::ofstream dump (file, std::ios::out | std::ios::binary);
	    if (dump.is_open() == false)
	    	return false;
	    dump.write ((const char *)pmem.getPointer(), pmem.size());

	    return true;
	}

	/*! \brief Load this object from file
	 *
	 * \param file filename
	 *
	 * \return true if succed
	 *
	 */
	bool load(const std::string & file)
	{
	    std::ifstream fs (file, std::ios::in | std::ios::binary | std::ios::ate );
	    if (fs.is_open() == false)
	    	return false;

	    // take the size of the file
	    size_t sz = fs.tellg();

	    fs.close();

	    // reopen the file without ios::ate to read
	    std::ifstream input (file, std::ios::in | std::ios::binary );
	    if (input.is_open() == false)
	    	return false;

	    // Create the HeapMemory and the ExtPreAlloc memory
	    std::vector<size_t> pap_prp;
	    pap_prp.push_back(sz);
	    HeapMemory pmem;
		ExtPreAlloc<HeapMemory> mem(pap_prp,pmem);

		// read
	    input.read((char *)pmem.getPointer(), sz);

	    //close the file
	    input.close();

		//Unpacking
		Unpack_stat ps;

	 	Unpacker<decltype(trpl),HeapMemory>::unpack(mem,trpl,ps);

	 	return true;
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
};



#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_ */
