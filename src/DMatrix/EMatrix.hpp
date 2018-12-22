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

template<typename _Scalar, int _Rows, int _Cols,
         int _Options = Eigen::AutoAlign |
#if EIGEN_GNUC_AT(3,4)
    // workaround a bug in at least gcc 3.4.6
    // the innermost ?: ternary operator is misparsed. We write it slightly
    // differently and this makes gcc 3.4.6 happy, but it's ugly.
    // The error would only show up with EIGEN_DEFAULT_TO_ROW_MAJOR is defined
    // (when EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION is RowMajor)
                          ( (_Rows==1 && _Cols!=1) ? Eigen::RowMajor
                          : !(_Cols==1 && _Rows!=1) ?  EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION
                          : Eigen::ColMajor ),
#else
                          ( (_Rows==1 && _Cols!=1) ? Eigen::RowMajor
                          : (_Cols==1 && _Rows!=1) ? Eigen::ColMajor
                          : EIGEN_DEFAULT_MATRIX_STORAGE_ORDER_OPTION ),
#endif
         int _MaxRows = _Rows,
         int _MaxCols = _Cols
> class EMatrix;

namespace Eigen {

namespace internal {
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct traits<EMatrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> >
{
private:
  enum { size = internal::size_at_compile_time<_Rows,_Cols>::ret };
  typedef typename find_best_packet<_Scalar,size>::type PacketScalar;
  enum {
      row_major_bit = _Options&RowMajor ? RowMajorBit : 0,
      is_dynamic_size_storage = _MaxRows==Dynamic || _MaxCols==Dynamic,
      max_size = is_dynamic_size_storage ? Dynamic : _MaxRows*_MaxCols,
      default_alignment = compute_default_alignment<_Scalar,max_size>::value,
      actual_alignment = ((_Options&DontAlign)==0) ? default_alignment : 0,
      required_alignment = unpacket_traits<PacketScalar>::alignment,
      packet_access_bit = (packet_traits<_Scalar>::Vectorizable && (EIGEN_UNALIGNED_VECTORIZE || (actual_alignment>=required_alignment))) ? PacketAccessBit : 0
    };

public:
  typedef _Scalar Scalar;
  typedef Dense StorageKind;
  typedef Eigen::Index StorageIndex;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = _Rows,
    ColsAtCompileTime = _Cols,
    MaxRowsAtCompileTime = _MaxRows,
    MaxColsAtCompileTime = _MaxCols,
    Flags = compute_matrix_flags<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::ret,
    Options = _Options,
    InnerStrideAtCompileTime = 1,
    OuterStrideAtCompileTime = (Options&RowMajor) ? ColsAtCompileTime : RowsAtCompileTime,

    // FIXME, the following flag in only used to define NeedsToAlign in PlainObjectBase
    EvaluatorFlags = LinearAccessBit | DirectAccessBit | packet_access_bit | row_major_bit,
    Alignment = actual_alignment
  };
};
}
}

/*! \brief it is just a wrapper for Eigen matrix compatible for openfpm serialization
 *
 *
 */
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
class EMatrix : public Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>
{

	typedef typename Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>::Base Base;

public:

	//typedef typename Eigen::internal::traits<Eigen::Matrix<_Scalar,_Rows,_Cols>>::Scalar Scalar;

	EIGEN_DENSE_PUBLIC_INTERFACE(EMatrix)

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

	////// wrap all eigen operator

    /** \brief Default constructor.
      *
      * For fixed-size matrices, does nothing.
      *
      * For dynamic-size matrices, creates an empty matrix of size 0. Does not allocate any array. Such a matrix
      * is called a null matrix. This constructor is the unique way to create null matrices: resizing
      * a matrix to 0 is not supported.
      *
      * \sa resize(Index,Index)
      */
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix()
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>()
    {}

    /**
      * \brief Assigns matrices to each other.
      *
      * \note This is a special case of the templated operator=. Its purpose is
      * to prevent a default operator= from hiding the templated operator=.
      *
      * \callgraph
      */
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix& operator=(const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>& other)
    {
    	((Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> *)this)->operator=(other);
    	return *this;
    }

    /** \internal
      * \brief Copies the value of the expression \a other into \c *this with automatic resizing.
      *
      * *this might be resized to match the dimensions of \a other. If *this was a null matrix (not already initialized),
      * it will be initialized.
      *
      * Note that copying a row-vector into a vector (and conversely) is allowed.
      * The resizing, if any, is then done in the appropriate way so that row-vectors
      * remain row-vectors and vectors remain vectors.
      */
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix& operator=(const Eigen::DenseBase<OtherDerived>& other)
    {
    	((Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> *)this)->operator=(other);
    	return *this;
    }

    /**
      * \brief Copies the generic expression \a other into *this.
      * \copydetails DenseBase::operator=(const EigenBase<OtherDerived> &other)
      */
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix& operator=(const Eigen::EigenBase<OtherDerived> &other)
    {
    	((Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> *)this)->operator=(other);
    	return *this;
    }

    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix& operator=(const Eigen::ReturnByValue<OtherDerived>& func)
    {
    	((Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> *)this)->operator=(func);
    	return *this;
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    EIGEN_DEVICE_FUNC
    EMatrix(Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> && other) EIGEN_NOEXCEPT_IF(std::is_nothrow_move_constructible<Scalar>::value)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(other)
    {}

    EIGEN_DEVICE_FUNC
    EMatrix& operator=(Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> && other) EIGEN_NOEXCEPT_IF(std::is_nothrow_move_assignable<Scalar>::value)
    {return *this;}

#endif

    #ifndef EIGEN_PARSED_BY_DOXYGEN

    // This constructor is for both 1x1 matrices and dynamic vectors
    template<typename T>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE explicit EMatrix(const T& x)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(x)
    {}

    template<typename T0, typename T1>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix(const T0& x, const T1& y)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(x,y)
    {}

    #else
    /** \brief Constructs a fixed-sized matrix initialized with coefficients starting at \a data */
    EIGEN_DEVICE_FUNC
    explicit EMatrix(const Scalar *data)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(data)
	 {}

    /** \brief Constructs a vector or row-vector with given dimension. \only_for_vectors
      *
      * This is useful for dynamic-size vectors. For fixed-size vectors,
      * it is redundant to pass these parameters, so one should use the default constructor
      * Matrix() instead.
      *
      * \warning This constructor is disabled for fixed-size \c 1x1 matrices. For instance,
      * calling Matrix<double,1,1>(1) will call the initialization constructor: Matrix(const Scalar&).
      * For fixed-size \c 1x1 matrices it is therefore recommended to use the default
      * constructor Matrix() instead, especially when using one of the non standard
      * \c EIGEN_INITIALIZE_MATRICES_BY_{ZERO,\c NAN} macros (see \ref TopicPreprocessorDirectives).
      */
    EIGEN_STRONG_INLINE explicit EMatrix(Index dim)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(dim)
	 {}

    /** \brief Constructs an initialized 1x1 matrix with the given coefficient */
    EMatrix(const Scalar& x)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(x)
    {}
    /** \brief Constructs an uninitialized matrix with \a rows rows and \a cols columns.
      *
      * This is useful for dynamic-size matrices. For fixed-size matrices,
      * it is redundant to pass these parameters, so one should use the default constructor
      * Matrix() instead.
      *
      * \warning This constructor is disabled for fixed-size \c 1x2 and \c 2x1 vectors. For instance,
      * calling Matrix2f(2,1) will call the initialization constructor: Matrix(const Scalar& x, const Scalar& y).
      * For fixed-size \c 1x2 or \c 2x1 vectors it is therefore recommended to use the default
      * constructor Matrix() instead, especially when using one of the non standard
      * \c EIGEN_INITIALIZE_MATRICES_BY_{ZERO,\c NAN} macros (see \ref TopicPreprocessorDirectives).
      */
    EIGEN_DEVICE_FUNC
    EMatrix(Index rows, Index cols)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(rows,cols)
	 {}

    /** \brief Constructs an initialized 2D vector with given coefficients */
    EMatrix(const Scalar& x, const Scalar& y)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(x,y)
	 {};
    #endif

    /** \brief Constructs an initialized 3D vector with given coefficients */
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix(const Scalar& x, const Scalar& y, const Scalar& z)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(x,y,z)
    {}

    /** \brief Constructs an initialized 4D vector with given coefficients */
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix(const Scalar& x, const Scalar& y, const Scalar& z, const Scalar& w)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(x,y,z,w)
    {}


    /** \brief Copy constructor */
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix(const Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols> & other)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(other)
    {}

    /** \brief Copy constructor for generic expressions.
      * \sa MatrixBase::operator=(const EigenBase<OtherDerived>&)
      */
    template<typename OtherDerived>
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE EMatrix(const Eigen::EigenBase<OtherDerived> &other)
    :Eigen::Matrix<_Scalar,_Rows,_Cols,_Options,_MaxRows,_MaxCols>(other)
    {}
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
