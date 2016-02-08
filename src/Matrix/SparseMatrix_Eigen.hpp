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
struct triplet<T,EIGEN_TRIPLET>
{
	Eigen::Triplet<T,long int> triplet;

	long int & row()
	{
		size_t ptr = (size_t)&triplet.row();

		return (long int &)*(long int *)ptr;
	}

	long int & col()
	{
		size_t ptr = (size_t)&triplet.col();

		return (long int &)*(long int *)ptr;
	}

	T & value()
	{
		size_t ptr = (size_t)&triplet.value();

		return (T &)*(T *)ptr;
	}
};

/* ! \brief Sparse Matrix implementation, that map over Eigen
 *
 * \tparam T Type of the sparse Matrix store on each row,colums
 * \tparam id_t type of id
 * \tparam impl implementation
 *
 */
template<typename T, typename id_t>
class SparseMatrix<T,id_t,Eigen::SparseMatrix<T,0,id_t>>
{
public:

	//! Triplet implementation id
	typedef boost::mpl::int_<EIGEN_TRIPLET> triplet_impl;

	//! Triplet type
	typedef triplet<T,EIGEN_TRIPLET> triplet_type;

private:

	Eigen::SparseMatrix<T,0,id_t> mat;
	openfpm::vector<triplet_type> trpl;

	/*! \brief Call-back to allocate buffer to receive incoming elements (particles)
	 *
	 * \param msg_i size required to receive the message from i
	 * \param total_msg total size to receive from all the processors
	 * \param total_p the total number of processor that want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr a pointer to the vector_dist structure
	 *
	 * \return the pointer where to store the message for the processor i
	 *
	 */
	template<typename triplet> static void * msg_alloc(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		openfpm::vector<openfpm::vector<triplet> *> * trpl = (openfpm::vector<openfpm::vector<triplet> *> *)ptr;

		if (trpl == NULL)
			std::cerr << __FILE__ << ":" << __LINE__ << " Internal error this processor is not suppose to receive\n";

		// We need one more slot
		trpl->add();

		trpl->last()->resize(msg_i/sizeof(triplet));

		// return the pointer
		return trpl->last()->getPointer();
	}

	/*! \brief Assemble the matrix
	 *
	 * \param trpl Matrix in form of triplets
	 * \param mat Matrix to assemble
	 *
	 */
	template<typename triplet, typename mat_impl> void assembleMatrix(openfpm::vector<openfpm::vector<triplet> *> & trpl, SparseMatrix<double,int,mat_impl> & mat)
	{
		Vcluster & vcl = *global_v_cluster;

		////// On Master and only here
		// we assemble the Matrix from the collected data
		if (vcl.getProcessingUnits() != 1)
		{
			// count the total triplet we have
			size_t tot_trpl = 0;

			for (size_t i = 0 ; i < trpl.size() ; i++)
				tot_trpl += trpl.get(i).size();

			openfpm::vector<triplet> mat_t;
			mat_t.resize(tot_trpl);

			// id zero
			size_t id = 0;

			// Add all the received triplet in one array
			for (size_t i = 0 ; i < trpl.size() ; i++)
			{
				for (size_t j = 0 ; j < trpl.get(i).size() ; j++)
				{
					mat_t.get(id) = trpl.get(i).get(j);
					id++;
				}
			}

			mat.setFromTriplets(mat_t);
		}
		else
		{
			//
			if (trpl.size() != 1)
				std::cerr << "Internal error: " << __FILE__ << ":" << __LINE__ << " in case of single processor we must have a single triplet set\n";

			mat.setFromTriplets(*trpl.get(0));
		}
	}

	/*! \brief Here we collect the full matrix on master
	 *
	 */
	void collect()
	{
		Vcluster & vcl = global_v_cluster;

		// If we are on master collect the information
		if (vcl.getProcessUnitID() == 0)
		{
			// send buffer (master does not send anything) so send req and send_buf
			// remain buffer with size 0
			openfpm::vector<size_t> send_req;
			openfpm::vector<openfpm::vector<triplet_type>> send_buf;

			// for each processor we are going to receive a set of triplet
			openfpm::vector<openfpm::vector<triplet_type>> trpl;

			// Send and recv multiple messages
			vcl.sendrecvMultipleMessagesNBX(send_req, send_buf,msg_alloc<triplet>,&trpl);

			assembleMatrix<triplet,Eigen::SparseMatrix<T,0,id_t>>(trpl);
		}
		else
		{
			// send buffer (master does not send anything) so send req and send_buf
			// remain buffer with size 0
			openfpm::vector<size_t> send_req;
			send_req.add(0);
			openfpm::vector<openfpm::vector<triplet_type> *> send_buf;
			send_buf.add(&A);

			// for each processor we are going to receive a set of triplet
			openfpm::vector<openfpm::vector<triplet_type>> trpl;

			// Send and recv multiple messages
			vcl.sendrecvMultipleMessagesNBX(send_req, send_buf,msg_alloc<triplet_type>,NULL);
		}
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

	/*! \brief Create a Matrix from a set of triplet
	 *
	 * \param N1 number of row
	 * \param N2 number of colums
	 * \param trpl triplet set
	 *
	 */
	SparseMatrix(size_t N1, size_t N2, const openfpm::vector<Eigen::Triplet<T,id_t>> & trpl)
	:mat(N1,N2)
	{
		this->trpl = trpl;
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
		collect();
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
		collect();
		assemble();

		return mat;
	}

	/*! \brief Resize the Sparse Matrix
	 *
	 * \param row number for row
	 * \param col number of colums
	 *
	 */
	void resize(size_t row, size_t col)
	{
		mat.resize(row,col);
	}
};



#endif /* OPENFPM_NUMERICS_SRC_MATRIX_SPARSEMATRIX_EIGEN_HPP_ */
