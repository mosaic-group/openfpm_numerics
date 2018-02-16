/*
 * EMatrix.hpp
 *
 *  Created on: Feb 12, 2018
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPACE_EMATRIX_HPP_
#define OPENFPM_DATA_SRC_SPACE_EMATRIX_HPP_

#include "config.h"

#ifdef HAVE_EIGEN
#include <Eigen/Dense>
#include "memory/ExtPreAlloc.hpp"
#include "memory/HeapMemory.hpp"
#include "Packer_Unpacker/Packer_util.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"

/*! \brief it is just a wrapper for Eigen matrix compatible for openfpm serialization
 *
 *
 */
template<typename _Scalar, int _Rows, int _Cols>
class EMatrix : public Eigen::Matrix<_Scalar,_Rows,_Cols>
{
public:

	/*! \brief It calculate the number of byte required to serialize the object
	 *
	 * \tparam prp list of properties
	 *
	 * \param req reference to the total counter required to pack the information
	 *
	 */
	template<int ... prp> inline void packRequest(size_t & req) const
	{
		// Memory required to serialize the Matrix
		req += sizeof(_Scalar)*this->rows()*this->cols() + 2*sizeof(_Scalar);
	}


	/*! \brief pack a vector selecting the properties to pack
	 *
	 * \param mem preallocated memory where to pack the vector
	 * \param sts pack-stat info
	 *
	 */
	template<int ... prp> inline void pack(ExtPreAlloc<HeapMemory> & mem, Pack_stat & sts) const
	{
		//Pack the number of rows and colums
		Packer<size_t, HeapMemory>::pack(mem,this->rows(),sts);
		Packer<size_t, HeapMemory>::pack(mem,this->cols(),sts);

		mem.allocate(sizeof(_Scalar)*this->rows()*this->cols());
		void * dst = mem.getPointer();
		void * src = (void *)this->data();

		memcpy(dst,src,sizeof(_Scalar)*this->rows()*this->cols());
		sts.incReq();
	}

	/*! \brief unpack a vector
	 *
	 * \param mem preallocated memory from where to unpack the vector
	 * \param ps unpack-stat info
	 */
	template<int ... prp> inline void unpack(ExtPreAlloc<HeapMemory> & mem, Unpack_stat & ps)
	{
		size_t row;
		size_t col;

		//Pack the number of rows and colums
		Unpacker<size_t, HeapMemory>::unpack(mem,row,ps);
		Unpacker<size_t, HeapMemory>::unpack(mem,col,ps);

		this->resize(row,col);

		void * dst = this->data();
		void * src = mem.getPointerOffset(ps.getOffset());

		memcpy(dst,src,sizeof(_Scalar)*row*col);

		ps.addOffset(sizeof(_Scalar)*row*col);
	}

};

// Standard typedef from eigen
typedef EMatrix< std::complex< double >, 2, 2 > 	EMatrix2cd;
typedef EMatrix< std::complex< float >, 2, 2 > 	EMatrix2cf;
typedef EMatrix< double, 2, 2 > 	EMatrix2d;
typedef EMatrix< float, 2, 2 > 	EMatrix2f;
typedef EMatrix< int, 2, 2 > 	EMatrix2i;
typedef EMatrix< std::complex< double >, 2, Eigen::Dynamic > 	EMatrix2Xcd;
typedef EMatrix< std::complex< float >, 2, Eigen::Dynamic > 	EMatrix2Xcf;
typedef EMatrix< double, 2, Eigen::Dynamic > 	EMatrix2Xd;
typedef EMatrix< float, 2, Eigen::Dynamic > 	EMatrix2Xf;
typedef EMatrix< int, 2, Eigen::Dynamic > 	EMatrix2Xi;
typedef EMatrix< std::complex< double >, 3, 3 > 	EMatrix3cd;
typedef EMatrix< std::complex< float >, 3, 3 > 	EMatrix3cf;
typedef EMatrix< double, 3, 3 > 	EMatrix3d;
typedef EMatrix< float, 3, 3 > 	EMatrix3f;
typedef EMatrix< int, 3, 3 > 	EMatrix3i;
typedef EMatrix< std::complex< double >, 3, Eigen::Dynamic > 	EMatrix3Xcd;
typedef EMatrix< std::complex< float >, 3, Eigen::Dynamic > 	EMatrix3Xcf;
typedef EMatrix< double, 3, Eigen::Dynamic > 	EMatrix3Xd;
typedef EMatrix< float, 3, Eigen::Dynamic > 	EMatrix3Xf;
typedef EMatrix< int, 3, Eigen::Dynamic > 	EMatrix3Xi;
typedef EMatrix< std::complex< double >, 4, 4 > 	EMatrix4cd;
typedef EMatrix< std::complex< float >, 4, 4 > 	EMatrix4cf;
typedef EMatrix< double, 4, 4 > 	EMatrix4d;
typedef EMatrix< float, 4, 4 > 	EMatrix4f;
typedef EMatrix< int, 4, 4 > 	EMatrix4i;
typedef EMatrix< std::complex< double >, 4, Eigen::Dynamic > 	EMatrix4Xcd;
typedef EMatrix< std::complex< float >, 4, Eigen::Dynamic > 	EMatrix4Xcf;
typedef EMatrix< double, 4, Eigen::Dynamic > 	EMatrix4Xd;
typedef EMatrix< float, 4, Eigen::Dynamic > 	EMatrix4Xf;
typedef EMatrix< int, 4, Eigen::Dynamic > 	EMatrix4Xi;
typedef EMatrix< std::complex< double >, Eigen::Dynamic, 2 > 	EMatrixX2cd;
typedef EMatrix< std::complex< float >, Eigen::Dynamic, 2 > 	EMatrixX2cf;
typedef EMatrix< double, Eigen::Dynamic, 2 > 	EMatrixX2d;
typedef EMatrix< float, Eigen::Dynamic, 2 > 	EMatrixX2f;
typedef EMatrix< int, Eigen::Dynamic, 2 > 	EMatrixX2i;
typedef EMatrix< std::complex< double >, Eigen::Dynamic, 3 > 	EMatrixX3cd;
typedef EMatrix< std::complex< float >, Eigen::Dynamic, 3 > 	EMatrixX3cf;
typedef EMatrix< double, Eigen::Dynamic, 3 > 	EMatrixX3d;
typedef EMatrix< float, Eigen::Dynamic, 3 > 	EMatrixX3f;
typedef EMatrix< int, Eigen::Dynamic, 3 > 	EMatrixX3i;
typedef EMatrix< std::complex< double >, Eigen::Dynamic, 4 > 	EMatrixX4cd;
typedef EMatrix< std::complex< float >, Eigen::Dynamic, 4 > 	EMatrixX4cf;
typedef EMatrix< double, Eigen::Dynamic, 4 > 	EMatrixX4d;
typedef EMatrix< float, Eigen::Dynamic, 4 > 	EMatrixX4f;
typedef EMatrix< int, Eigen::Dynamic, 4 > 	EMatrixX4i;
typedef EMatrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic > 	EMatrixXcd;
typedef EMatrix< std::complex< float >, Eigen::Dynamic, Eigen::Dynamic > 	EMatrixXcf;
typedef EMatrix< double, Eigen::Dynamic, Eigen::Dynamic > 	EMatrixXd;
typedef EMatrix< float, Eigen::Dynamic, Eigen::Dynamic > 	EMatrixXf;
typedef EMatrix< int, Eigen::Dynamic, Eigen::Dynamic > 	EMatrixXi;
typedef EMatrix< std::complex< double >, 1, 2 > 	ERowVector2cd;
typedef EMatrix< std::complex< float >, 1, 2 > 	ERowVector2cf;
typedef EMatrix< double, 1, 2 > 	ERowVector2d;
typedef EMatrix< float, 1, 2 > 	ERowVector2f;
typedef EMatrix< int, 1, 2 > 	ERowVector2i;
typedef EMatrix< std::complex< double >, 1, 3 > 	ERowVector3cd;
typedef EMatrix< std::complex< float >, 1, 3 > 	ERowVector3cf;
typedef EMatrix< double, 1, 3 > 	ERowVector3d;
typedef EMatrix< float, 1, 3 > 	ERowVector3f;
typedef EMatrix< int, 1, 3 > 	ERowVector3i;
typedef EMatrix< std::complex< double >, 1, 4 > 	ERowVector4cd;
typedef EMatrix< std::complex< float >, 1, 4 > 	ERowVector4cf;
typedef EMatrix< double, 1, 4 > 	ERowVector4d;
typedef EMatrix< float, 1, 4 > 	ERowVector4f;
typedef EMatrix< int, 1, 4 > 	ERowVector4i;
typedef EMatrix< std::complex< double >, 1, Eigen::Dynamic > 	ERowVectorXcd;
typedef EMatrix< std::complex< float >, 1, Eigen::Dynamic > 	ERowVectorXcf;
typedef EMatrix< double, 1, Eigen::Dynamic > 	ERowVectorXd;
typedef EMatrix< float, 1, Eigen::Dynamic > 	ERowVectorXf;
typedef EMatrix< int, 1, Eigen::Dynamic > 	ERowVectorXi;
typedef EMatrix< std::complex< double >, 2, 1 > 	EVector2cd;
typedef EMatrix< std::complex< float >, 2, 1 > 	EVector2cf;
typedef EMatrix< double, 2, 1 > 	EVector2d;
typedef EMatrix< float, 2, 1 > 	EVector2f;
typedef EMatrix< int, 2, 1 > 	EVector2i;
typedef EMatrix< std::complex< double >, 3, 1 > 	EVector3cd;
typedef EMatrix< std::complex< float >, 3, 1 > 	EVector3cf;
typedef EMatrix< double, 3, 1 > 	EVector3d;
typedef EMatrix< float, 3, 1 > 	EVector3f;
typedef EMatrix< int, 3, 1 > 	EVector3i;
typedef EMatrix< std::complex< double >, 4, 1 > 	EVector4cd;
typedef EMatrix< std::complex< float >, 4, 1 > 	EVector4cf;
typedef EMatrix< double, 4, 1 > 	EVector4d;
typedef EMatrix< float, 4, 1 > 	EVector4f;
typedef EMatrix< int, 4, 1 > 	EVector4i;
typedef EMatrix< std::complex< double >, Eigen::Dynamic, 1 > 	EVectorXcd;
typedef EMatrix< std::complex< float >, Eigen::Dynamic, 1 > 	EVectorXcf;
typedef EMatrix< double, Eigen::Dynamic, 1 > 	EVectorXd;
typedef EMatrix< float, Eigen::Dynamic, 1 > 	EVectorXf;
typedef EMatrix< int, Eigen::Dynamic, 1 > 	EVectorXi;

#else

// There is not EMatrix if we do not have Eigen

#endif

#endif /* OPENFPM_DATA_SRC_SPACE_EMATRIX_HPP_ */
