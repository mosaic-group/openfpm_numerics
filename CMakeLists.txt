cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_numerics LANGUAGES C CXX)


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake_modules/)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(PETSC_ROOT CACHE PATH "If compiling with linear algebra indicate the PETSC root directory")
set(PARMETIS_ROOT CACHE PATH "Parmetis root directory")
set(METIS_ROOT CACHE PATH "Metis root directory")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(HDF5_ROOT CACHE PATH "HDF5 root path")
set(EIGEN3_ROOT CACHE PATH "Eigen3 include path")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS2 CACHE BOOL "Activate compilation with SE_CLASS2")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(PROFILE_WITH_SCOREP CACHE BOOL "Enable profiling with scorep")
set(ENV{PETSC_DIR} ${PETSC_ROOT})
set(ENV{HDF5_ROOT} ${HDF5_ROOT})
set(ENV{EIGEN3_ROOT} ${EIGEN3_ROOT})
set(METIS_DIR ${METIS_ROOT})
set(PARMETIS_DIR ${PARMETIS_ROOT})

set(ENV{PATH} "$ENV{PATH}:${HDF5_ROOT}/bin")
set(HDF5_PREFER_PARALLEL TRUE)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CUDA_STANDARD 11)

if(ENABLE_GPU)
        enable_language(CUDA)
        find_package(CUDA)

        if (CUDA_VERSION_MAJOR EQUAL 9 AND CUDA_VERSION_MINOR EQUAL 2)
                message("CUDA is compatible")
                set(WARNING_SUPPRESSION_AND_OPTION_NVCC  -Xcudafe "--display_error_number --diag_suppress=611 --diag_suppress=2885 --diag_suppress=2886  --diag_suppress=2887  --diag_suppress=2888 --diag_suppress=186 --diag_suppress=111" --expt-extended-lambda)
                FILE(WRITE cuda_options " -Xcudafe \"--display_error_number --diag_suppress=611 --diag_suppress=2885 --diag_suppress=2886  --diag_suppress=2887  --diag_suppress=2888 --diag_suppress=186 --diag_suppress=111\" --expt-extended-lambda ")
        elseif ( CUDA_VERSION_MAJOR EQUAL 10 AND CUDA_VERSION_MINOR EQUAL 1 )
                message("CUDA is compatible")
                set(WARNING_SUPPRESSION_AND_OPTION_NVCC  -Xcudafe "--display_error_number --diag_suppress=2915 --diag_suppress=2914 --diag_suppress=2912 --diag_suppress=2913 --diag_suppress=111 --diag_suppress=186 --diag_suppress=611 " --expt-extended-lambda )
                FILE(WRITE cuda_options "-Xcudafe \"--display_error_number --diag_suppress=2915 --diag_suppress=2914  --diag_suppress=2912 --diag_suppress=2913 --diag_suppress=111 --diag_suppress=186 --diag_suppress=611 \" --expt-extended-lambda")
        else()
                message(FATAL_ERROR "CUDA is incompatible, version 9.2 and 10.1 is only supported")
        endif()
endif()


find_package(Boost 1.66.0 REQUIRED COMPONENTS unit_test_framework iostreams program_options)
find_package(MPI REQUIRED)
find_package(PETSc)
find_package(HDF5 REQUIRED)
find_package(Eigen3)
find_package(LibHilbert REQUIRED)
find_package(Metis REQUIRED)
find_package(ParMetis REQUIRED)

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

if(SE_CLASS2)
	set(DEFINE_SE_CLASS2 "#define SE_CLASS2")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(CUDA_FOUND)
	set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if (METIS_FOUND)
	set(DEFINE_HAVE_METIS "#define HAVE_METIS 1")
else()
	message( FATAL_ERROR "Metis is required in order to install OpenFPM" )
endif()

if (PARMETIS_FOUND)
	set(DEFINE_HAVE_PARMETIS "#define HAVE_PARMETIS 1")
else()
	message( FATAL_ERROR "ParMetis is required in order to install OpenFPM")
endif()

if(MPI_FOUND)
	set(DEFINE_HAVE_MPI "#define HAVE_MPI")
else()
	message( FATAL_ERROR "MPI is required in order to install OpenFPM" )
endif()

if(PETSC_FOUND)
        set(DEFINE_HAVE_PETSC "#define HAVE_PETSC")
endif()

if (Boost_FOUND)
	set(DEFINE_HAVE_BOOST "#define HAVE_BOOST")
	set(DEFINE_HAVE_BOOST_IOSTREAMS "#define HAVE_BOOST_IOSTREAMS")
	set(DEFINE_HAVE_BOOST_PROGRAM_OPTIONS "#define HAVE_BOOST_PROGRAM_OPTIONS")
	set(DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK "#define HAVE_BOOST_UNIT_TEST_FRAMEWORK")
else()
	message( FATAL_ERROR "BOOST is required in order to install OpenFPM" )
endif()

if(HDF5_FOUND)
	if (HDF5_IS_PARALLEL)
		set(DEFINE_HAVE_HDF5 "#define HAVE_HDF5")
	else()
		message( STATUS "HDF5 found ${HDF5_INCLUDE_DIRS} does not have parallel support, OpenFPM require it" )
	endif()
else()
	message( FATAL_ERROR "HDF5 with parallel support is required in order to install OpenFPM" )
endif()

if(EIGEN_FOUND)
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.h)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory (src)

