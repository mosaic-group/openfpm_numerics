cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

if (POLICY CMP0074)
        cmake_policy(SET CMP0074 OLD)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake_modules/)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(PETSC_ROOT CACHE PATH "If compiling with linear algebra indicate the PETSC root directory")
set(PARMETIS_ROOT CACHE PATH "Parmetis root directory")
set(METIS_ROOT CACHE PATH "Metis root directory")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(HDF5_ROOT CACHE PATH "HDF5 root path")
set(EIGEN3_ROOT CACHE PATH "Eigen3 include path")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(SUITESPARSE_ROOT CACHE PATH "The suitesparse root path")
set(TINYOBJLOADER_ROOT CACHE PATH "TinyObjLoader library path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(BLITZ_ROOT CACHE PATH "Blitz root directory")
set(ALGOIM_ROOT CACHE PATH "Algoim root directory")
set(PROFILE_WITH_SCOREP CACHE BOOL "Enable profiling with scorep")
set(ENV{PETSC_DIR} ${PETSC_ROOT})
set(ENV{HDF5_ROOT} ${HDF5_ROOT})
set(ENV{EIGEN3_ROOT} ${EIGEN3_ROOT})
set(METIS_DIR ${METIS_ROOT})
set(PARMETIS_DIR ${PARMETIS_ROOT})
set(OPENBLAS_ROOT CACHE PATH "Root path for blas library")
set(SuiteSparse_ROOT ${SUITESPARSE_ROOT})
set(ALPAKA_ROOT CACHE PATH "Alpaka root path")

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CUDA_STANDARD 14)

if (OPENBLAS_ROOT)
	set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${OPENBLAS_ROOT}/lib")
	set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:${OPENBLAS_ROOT}/lib")
endif()

set(ENV{PATH} "$ENV{PATH}:${HDF5_ROOT}/bin")
set(HDF5_PREFER_PARALLEL TRUE)

set (ENV{BOOST_ROOT} ${BOOST_ROOT})
set (Boost_NO_BOOST_CMAKE OFF)
find_package(Boost 1.72.0 REQUIRED COMPONENTS unit_test_framework iostreams program_options OPTIONAL_COMPONENTS fiber context)
find_package(MPI REQUIRED)
find_package(PETSc)
find_package(HDF5)
find_package(Eigen3)
find_package(LibHilbert)
find_package(TinyObjLoader )
find_package(BLAS)
find_package(LAPACK)
find_package(Eigen3)
find_package(SuiteSparse OPTIONAL_COMPONENTS UMFPACK)
find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system thread)

set(CMAKE_SKIP_BUILD_RPATH TRUE)

if(PROFILE_WITH_SCOREP)
	set(CMAKE_CXX_COMPILER_LAUNCHER "scorep")
	set(CMAKE_CC_COMPILER_LAUNCHER "scorep")
	set(CMAKE_CUDA_COMPILER_LAUNCHER "scorep")
endif()


if(CUDA_FOUND)
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cuda.cu")
else()
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cpu.cpp")
endif()
###### CONFIG.h FILE ######

if(SE_CLASS1)
	set(DEFINE_SE_CLASS1 "#define SE_CLASS1")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(PETSC_FOUND)
	set(DEFINE_HAVE_PETSC "#define HAVE_PETSC")
endif()

if(SCAN_COVERTY)
	set(DEFINE_SCAN_COVERTY "#define COVERTY_SCAN")
endif()

if(CUDA_FOUND)
	set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if (METIS_FOUND)
	set(DEFINE_HAVE_METIS "#define HAVE_METIS 1")
else()
	file(WRITE error_code "201")
	message( FATAL_ERROR "Metis is required in order to install OpenFPM" )
endif()

if (PARMETIS_FOUND)
	set(DEFINE_HAVE_PARMETIS "#define HAVE_PARMETIS 1")
else()
	file(WRITE error_code "208")
	message( FATAL_ERROR "ParMetis is required in order to install OpenFPM")
endif()

if(MPI_FOUND)
	set(DEFINE_HAVE_MPI "#define HAVE_MPI")
else()
	file(WRITE error_code "200")
	message( FATAL_ERROR "MPI is required in order to install OpenFPM" )
endif()

if(PETSC_FOUND)
        set(DEFINE_HAVE_PETSC "#define HAVE_PETSC")
endif()

if(HDF5_FOUND)
	if (HDF5_IS_PARALLEL)
		set(DEFINE_HAVE_HDF5 "#define HAVE_HDF5")
	else()
		file(WRITE error_code "207")
		message( FATAL_ERROR "HDF5 found ${HDF5_INCLUDE_DIRS} does not have parallel support, OpenFPM require it" )
	endif()
else()
	file(WRITE error_code "207")
	message( FATAL_ERROR "HDF5 with parallel support is required in order to install OpenFPM" )
endif()

if(EIGEN3_FOUND)
	set(DEFINE_HAVE_EIGEN "#define HAVE_EIGEN")
endif()


if(EIGEN3_FOUND)
        set(DEFINE_HAVE_EIGEN "#define HAVE_EIGEN")
endif()

if(LIBHILBERT_FOUND)
        set(DEFINE_HAVE_LIBHILBERT "#define HAVE_LIBHILBERT 1")
else()
        file(WRITE error_code "210")
        message( FATAL_ERROR "LibHilbert is required in order to install OpenFPM")
endif()

if(SUITESPARSE_FOUND AND SuiteSparse_UMFPACK_FOUND)
        set(DEFINE_HAVE_SUITESPARSE "#define HAVE_SUITESPARSE")
endif()

if(SUITESPARSE_FOUND AND SuiteSparse_UMFPACK_FOUND)
	set(DEFINE_HAVE_SUITESPARSE "#define HAVE_SUITESPARSE")
endif()

if(APPLE)
	set(DEFINE_HAVE_OSX "#define HAVE_OSX")
endif()

if(TINYOBJLOADER_FOUND)
	set(DEFINE_HAVE_TINYOBJLOADER "#define HAVE_TINYOBJLOADER 1")
endif()

if (TEST_COVERAGE)
        set(DEFINE_TEST_COVERAGE_MODE "#define TEST_COVERAGE_MODE")
endif()

file(WRITE error_code "0")
file(WRITE cuda_lib "${CUDA_cudart_static_LIBRARY} ${CUDA_cudadevrt_LIBRARY}")
file(WRITE cuda_include "-I${CUDA_INCLUDE_DIRS}")
file(WRITE mpi_include "-I${MPI_C_INCLUDE_DIRS}")
file(WRITE mpi_libs "${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.h)

if (CUDA_ON_BACKEND STREQUAL "CUDA")
        enable_language(CUDA)
endif()

add_subdirectory (src)


