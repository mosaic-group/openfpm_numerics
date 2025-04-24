cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 OLD)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake_modules/)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CUDA_STANDARD 14)

if (OPENBLAS_ROOT)
	set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${OPENBLAS_ROOT}/lib")
	set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:${OPENBLAS_ROOT}/lib")
endif()

if(CUDA_FOUND)
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cuda.cu")
else()
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cpu.cpp")
endif()

if(PETSC_FOUND)
	set(DEFINE_HAVE_PETSC "#define HAVE_PETSC")
endif()

if(EIGEN3_FOUND)
	set(DEFINE_HAVE_EIGEN "#define HAVE_EIGEN")
endif()

if(SUITESPARSE_FOUND AND SuiteSparse_UMFPACK_FOUND)
	set(DEFINE_HAVE_SUITESPARSE "#define HAVE_SUITESPARSE")
endif()

add_subdirectory (src)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(DEFINE_HAVE_SUITESPARSE ${DEFINE_HAVE_SUITESPARSE} CACHE INTERNAL "")
	set(DEFINE_HAVE_EIGEN ${DEFINE_HAVE_EIGEN} CACHE INTERNAL "")
	set(DEFINE_HAVE_PETSC ${DEFINE_HAVE_PETSC} CACHE INTERNAL "")
endif()

