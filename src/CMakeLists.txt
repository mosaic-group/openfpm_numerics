cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables


if (CUDA_ON_CPU)
        add_definitions(-DCUDA_ON_CPU)
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDA_GPU -DCUDART_VERSION=11000")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
                add_definitions("-x c++")
        endif()
endif()

if (CUDA_FOUND)
	set(CUDA_SOURCES Operators/Vector/vector_dist_operators_unit_tests.cu
					 Operators/Vector/vector_dist_operators_apply_kernel_unit_tests.cu)
endif()

add_executable(numerics ${OPENFPM_INIT_FILE} ${CUDA_SOURCES}
						main.cpp 
						Matrix/SparseMatrix_unit_tests.cpp
						interpolation/interpolation_unit_tests.cpp
						Vector/Vector_unit_tests.cpp  
						Solvers/petsc_solver_unit_tests.cpp 
						FiniteDifference/FDScheme_unit_tests.cpp 
						FiniteDifference/eq_unit_test_3d.cpp 
						FiniteDifference/eq_unit_test.cpp  
						Operators/Vector/vector_dist_operators_unit_tests.cpp 
						Operators/Vector/vector_dist_operators_apply_kernel_unit_tests.cpp
						../../src/lib/pdata.cpp)

add_dependencies(numerics ofpmmemory)
add_dependencies(numerics vcluster)

###########################

if ( CMAKE_COMPILER_IS_GNUCC )
    target_compile_options(numerics PRIVATE "-Wno-deprecated-declarations")
    if (TEST_COVERAGE)
        target_compile_options(numerics PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fprofile-arcs -ftest-coverage>)
    endif()
endif()

if (CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
        add_definitions(-D__STRICT_ANSI__)
endif()

if(CUDA_FOUND)
        if (TEST_COVERAGE)
                target_compile_options(numerics PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fprofile-arcs -ftest-coverage" >)
        endif()
endif()

target_include_directories (numerics PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_devices/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_vcluster/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_data/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (numerics PUBLIC ${PARMETIS_ROOT}/include)
target_include_directories (numerics PUBLIC ${METIS_ROOT}/include)
target_include_directories (numerics PUBLIC ${HDF5_ROOT}/include)
target_include_directories (numerics PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
target_include_directories (numerics PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories (numerics PUBLIC ${Vc_INCLUDE_DIR})
target_include_directories (numerics PUBLIC ${ALPAKA_ROOT}/include)
if(EIGEN3_FOUND)
	target_include_directories (numerics PUBLIC ${EIGEN3_INCLUDE_DIR})
endif()

link_directories(${PARMETIS_ROOT} ${METIS_ROOT})
target_link_libraries(numerics ${Boost_LIBRARIES})
target_link_libraries(numerics ${HDF5_LIBRARIES})
target_link_libraries(numerics ${PARMETIS_LIBRARIES})
target_link_libraries(numerics -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(numerics ${Vc_LIBRARIES})
if(PETSC_FOUND)
	target_include_directories (numerics PUBLIC ${PETSC_INCLUDES})
	target_link_libraries(numerics ${PETSC_LIBRARIES})
endif()
target_link_libraries(numerics vcluster)
target_link_libraries(numerics ofpmmemory)

if(SuiteSparse_FOUND)
	target_include_directories (numerics PUBLIC ${SuiteSparse_INCLUDE_DIRS})
        target_link_libraries(numerics ${SuiteSparse_LIBRARIES})
	target_include_directories (numerics PUBLIC ${BLAS_INCLUDE_DIRS})
	target_link_libraries(numerics ${BLAS_LIBRARIES})
endif()

target_link_libraries(numerics -L${METIS_ROOT}/lib metis)

if (TEST_COVERAGE)
    target_link_libraries(numerics -lgcov)
endif()

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(numerics PUBLIC cxx_std_11)
target_link_libraries(numerics ${MPI_C_LIBRARIES})
if (NOT APPLE)
    target_link_libraries(numerics rt)
endif ()

install(FILES Matrix/SparseMatrix.hpp 
	      Matrix/SparseMatrix_Eigen.hpp 
	      Matrix/SparseMatrix_petsc.hpp
	      DESTINATION openfpm_numerics/include/Matrix
	      COMPONENT OpenFPM)

install(FILES Vector/Vector_eigen.hpp
	      Vector/Vector_petsc.hpp 
	      Vector/Vector_util.hpp 
	      Vector/Vector.hpp
	      DESTINATION openfpm_numerics/include/Vector 
	      COMPONENT OpenFPM)

install(FILES Solvers/umfpack_solver.hpp 
	      Solvers/petsc_solver.hpp
	      Solvers/petsc_solver_AMG_report.hpp
	      DESTINATION openfpm_numerics/include/Solvers 
	      COMPONENT OpenFPM)

install(FILES util/petsc_util.hpp 
	      util/linalgebra_lib.hpp 
	      util/util_num.hpp 
	      util/grid_dist_testing.hpp
	      DESTINATION openfpm_numerics/include/util 
	      COMPONENT OpenFPM)

install(FILES FiniteDifference/Average.hpp 
      	      FiniteDifference/Derivative.hpp 
      	      FiniteDifference/FD_util_include.hpp  
      	      FiniteDifference/eq.hpp 
      	      FiniteDifference/FDScheme.hpp 
      	      FiniteDifference/Laplacian.hpp 
      	      FiniteDifference/mul.hpp 
      	      FiniteDifference/sum.hpp 
      	      DESTINATION openfpm_numerics/include/FiniteDifference 
	      COMPONENT OpenFPM)
	      
install(FILES FiniteDifference/util/common.hpp
      	      DESTINATION openfpm_numerics/include/FiniteDifference/util 
	      COMPONENT OpenFPM)

install(FILES PSE/Kernels.hpp PSE/Kernels_test_util.hpp 
        DESTINATION openfpm_numerics/include/PSE 
	COMPONENT OpenFPM)

install(FILES Operators/Vector/vector_dist_operators_extensions.hpp 
      	       Operators/Vector/vector_dist_operators.hpp 
      	       Operators/Vector/vector_dist_operators_apply_kernel.hpp 
      	       Operators/Vector/vector_dist_operators_functions.hpp 
      	       Operators/Vector/vector_dist_operator_assign.hpp
      	       DESTINATION openfpm_numerics/include/Operators/Vector 
	       COMPONENT OpenFPM)

       install(FILES Operators/Vector/cuda/vector_dist_operators_cuda.cuh
	      DESTINATION openfpm_numerics/include/Operators/Vector/cuda 
	      COMPONENT OpenFPM)

install(FILES Draw/DrawParticles.hpp 
      	      Draw/PointIterator.hpp 
      	      Draw/PointIteratorSkin.hpp
	      Draw/DrawCircle.hpp
	      Draw/DrawSphere.hpp
      	      DESTINATION openfpm_numerics/include/Draw 
	      COMPONENT OpenFPM)

install(FILES interpolation/interpolation.hpp 
      	      interpolation/mp4_kernel.hpp 
      	      interpolation/z_spline.hpp
      	      DESTINATION openfpm_numerics/include/interpolation 
	      COMPONENT OpenFPM)

install(FILES level_set/redistancing_Sussman/ComputeGradient.hpp
	      level_set/redistancing_Sussman/HelpFunctions.hpp
	      level_set/redistancing_Sussman/HelpFunctionsForGrid.hpp
	      level_set/redistancing_Sussman/NarrowBand.hpp
	      level_set/redistancing_Sussman/RedistancingSussman.hpp
	      DESTINATION openfpm_numerics/include/level_set/redistancing_Sussman
	      COMPONENT OpenFPM)

install(FILES DMatrix/EMatrix.hpp
	DESTINATION openfpm_numerics/include/DMatrix 
	COMPONENT OpenFPM)

#if(BUILD_TESTING)

#  add_executable(particle_test test.cu)

#  set_target_properties(particle_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#  target_link_libraries(particle_test PRIVATE particles)

#  add_test(NAME particles_10k COMMAND particle_test 10000 )
#  add_test(NAME particles_256k COMMAND particle_test 256000 )

#  if(APPLE)
    # We need to add the default path to the driver (libcuda.dylib) as an rpath,
    # so that the static cuda runtime can find it at runtime.
    #    set_property(TARGET particle_test PROPERTY BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    #  endif()
    #endif()

