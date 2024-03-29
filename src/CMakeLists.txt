cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables

if (NOT CUDA_ON_BACKEND STREQUAL "None")
	set(CUDA_SOURCES Operators/Vector/vector_dist_operators_unit_tests.cu
		Operators/Vector/vector_dist_operators_apply_kernel_unit_tests.cu)
endif()

if (CUDA_ON_BACKEND STREQUAL "CUDA")
	set(CUDA_SOURCES ${CUDA_SOURCES} DCPSE/DCPSE_op/tests/DCPSE_op_Solver_test.cu
		DCPSE/DCPSE_op/tests/DCPSE_op_test_base_tests.cu
		#DCPSE/DCPSE_op/tests/DCPSE_op_subset_test.cu
		#OdeIntegrators/tests/Odeintegrators_test_gpu.cu
		DCPSE/DCPSE_op/tests/DCPSE_op_test_temporal.cu)
endif()

if (CUDA_ON_BACKEND STREQUAL "CUDA")
	set(CUDA_SOURCES ${CUDA_SOURCES} DCPSE/DCPSE_op/tests/DCPSE_op_Solver_test.cu
		DCPSE/DCPSE_op/tests/DCPSE_op_test_base_tests.cu
		#DCPSE/DCPSE_op/tests/DCPSE_op_subset_test.cu
		#OdeIntegrators/tests/Odeintegrators_test_gpu.cu
		DCPSE/DCPSE_op/tests/DCPSE_op_test_temporal.cu)
endif()

if (CUDA_ON_BACKEND STREQUAL "OpenMP" OR CUDA_ON_BACKEND STREQUAL "SEQUENTIAL")
	set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)
	set_source_files_properties(${CUDA_SOURCES} PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDART_VERSION=11000")

	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
		add_definitions("-x c++")
	endif()
endif()

if ( CUDA_ON_BACKEND STREQUAL "HIP" AND HIP_FOUND )

	list(APPEND HIP_HIPCC_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		list(APPEND HIP_HIPCC_FLAGS -O0)
	endif()

	list(APPEND HIP_HIPCC_FLAGS -D__NVCC__ -D__HIP__  -DCUDART_VERSION=11000 -D__CUDACC__ -D__CUDACC_VER_MAJOR__=11 -D__CUDACC_VER_MINOR__=0 -D__CUDACC_VER_BUILD__=0 --std=c++14)
	set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)

	hip_add_executable(numerics ${OPENFPM_INIT_FILE} ${CUDA_SOURCES}
		regression/regression_test.cpp
		regression/poly_levelset_test.cpp
		OdeIntegrators/tests/OdeIntegratores_base_tests.cpp
		OdeIntegrators/tests/OdeIntegrator_grid_tests.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_subset_test.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_test_base_tests.cpp
		FiniteDifference/FD_Solver_test.cpp
		FiniteDifference/FD_op_Tests.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_test3d.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_Solver_test.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_test_temporal.cpp
		DCPSE/tests/Dcpse_unit_tests.cpp
		DCPSE/tests/DcpseRhs_unit_tests.cpp
		DCPSE/tests/MonomialBasis_unit_tests.cpp
		DCPSE/tests/Support_unit_tests.cpp
		DCPSE/tests/Vandermonde_unit_tests.cpp
		main.cpp
		Matrix/SparseMatrix_unit_tests.cpp
		#interpolation/interpolation_unit_tests.cpp
		Vector/Vector_unit_tests.cpp
		Solvers/petsc_solver_unit_tests.cpp
		FiniteDifference/FDScheme_unit_tests.cpp
		FiniteDifference/eq_unit_test_3d.cpp
		FiniteDifference/eq_unit_test.cpp
		FiniteDifference/tests/Eno_Weno_unit_test.cpp
		FiniteDifference/tests/Upwind_gradient_unit_test.cpp
		FiniteDifference/tests/FD_simple_unit_test.cpp
		Operators/Vector/vector_dist_operators_unit_tests.cpp
		Operators/Vector/vector_dist_operators_apply_kernel_unit_tests.cpp
		../../openfpm_pdata/src/lib/pdata.cpp
		#BoundaryConditions/tests/method_of_images_cylinder_unit_test.cpp
		#level_set/closest_point/closest_point_unit_tests.cpp
		level_set/redistancing_Sussman/tests/redistancingSussman_fast_unit_test.cpp
		#level_set/redistancing_Sussman/tests/help_functions_unit_test.cpp
		level_set/redistancing_Sussman/tests/narrowBand_unit_test.cpp
		#level_set/redistancing_Sussman/tests/redistancingSussman_unit_test.cpp
		#level_set/redistancing_Sussman/tests/convergence_test.cpp
		)

else()

	add_executable(numerics ${OPENFPM_INIT_FILE} ${CUDA_SOURCES}
		regression/regression_test.cpp
		regression/poly_levelset_test.cpp
		OdeIntegrators/tests/OdeIntegratores_base_tests.cpp
		OdeIntegrators/tests/OdeIntegrator_grid_tests.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_subset_test.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_test_base_tests.cpp
		FiniteDifference/FD_Solver_test.cpp
		FiniteDifference/FD_op_Tests.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_test3d.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_Solver_test.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_test_temporal.cpp
		DCPSE/tests/Dcpse_unit_tests.cpp
		DCPSE/tests/DcpseRhs_unit_tests.cpp
		DCPSE/tests/MonomialBasis_unit_tests.cpp
		DCPSE/tests/Support_unit_tests.cpp
		DCPSE/tests/Vandermonde_unit_tests.cpp
		main.cpp
		Matrix/SparseMatrix_unit_tests.cpp
		#interpolation/interpolation_unit_tests.cpp
		Vector/Vector_unit_tests.cpp
		Solvers/petsc_solver_unit_tests.cpp
		FiniteDifference/FDScheme_unit_tests.cpp
		FiniteDifference/eq_unit_test_3d.cpp
		FiniteDifference/eq_unit_test.cpp
		FiniteDifference/tests/Eno_Weno_unit_test.cpp
		FiniteDifference/tests/Upwind_gradient_unit_test.cpp
		FiniteDifference/tests/FD_simple_unit_test.cpp
		Operators/Vector/vector_dist_operators_unit_tests.cpp
		Operators/Vector/vector_dist_operators_apply_kernel_unit_tests.cpp
		../../openfpm_pdata/src/lib/pdata.cpp
		BoundaryConditions/tests/method_of_images_cylinder_unit_test.cpp
		#level_set/closest_point/closest_point_unit_tests.cpp
		#level_set/redistancing_Sussman/tests/redistancingSussman_unit_test.cpp
		#level_set/redistancing_Sussman/tests/convergence_test.cpp
		level_set/particle_cp/pcp_unit_tests.cpp
		DCPSE/DCPSE_op/tests/DCPSE_op_Surface_tests.cpp
		#BoundaryConditions/tests/method_of_images_cylinder_unit_test.cpp
		level_set/redistancing_Sussman/tests/redistancingSussman_fast_unit_test.cpp
		#level_set/redistancing_Sussman/tests/help_functions_unit_test.cpp
		level_set/redistancing_Sussman/tests/narrowBand_unit_test.cpp)

	set_property(TARGET numerics PROPERTY CUDA_ARCHITECTURES OFF)

endif()

if (HIP_FOUND)
	add_dependencies(numerics ofpmmemory_dl)
	add_dependencies(numerics vcluster_dl)
else()
	add_dependencies(numerics ofpmmemory)
	add_dependencies(numerics vcluster)
endif()

###########################
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(numerics PRIVATE "-Wno-undefined-var-template")
	target_compile_options(numerics PRIVATE "-Wno-macro-redefined")
endif()

if ( CMAKE_COMPILER_IS_GNUCC )
	#target_compile_options(numerics PRIVATE "-Wno-deprecated-declarations")
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
	target_link_libraries(numerics -lcublas)
endif()

target_include_directories (numerics PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_devices/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_vcluster/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_pdata/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_data/src/)
target_include_directories (numerics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../openfpm_io/src/)
target_include_directories (numerics PUBLIC ${CMAKE_BINARY_DIR}/config)
target_include_directories (numerics PUBLIC ${PARMETIS_ROOT}/include)
target_include_directories (numerics PUBLIC ${METIS_ROOT}/include)
target_include_directories (numerics PUBLIC ${HDF5_ROOT}/include)
target_include_directories (numerics PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
target_include_directories (numerics PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories (numerics PUBLIC ${Vc_INCLUDE_DIR})
target_include_directories (numerics PUBLIC ${BLITZ_ROOT}/include)
target_include_directories (numerics PUBLIC ${ALGOIM_ROOT}/include)
target_include_directories (numerics PUBLIC ${ALPAKA_ROOT}/include)
target_include_directories (numerics PUBLIC ${MPI_C_INCLUDE_DIRS})

if(EIGEN3_FOUND)
	target_include_directories (numerics PUBLIC /usr/local/include) 
	target_include_directories (numerics PUBLIC ${EIGEN3_INCLUDE_DIR})
endif()

link_directories(${PARMETIS_ROOT} ${METIS_ROOT})
target_link_libraries(numerics ${Boost_LIBRARIES})
target_link_libraries(numerics ${HDF5_LIBRARIES})
target_link_libraries(numerics ${PARMETIS_LIBRARIES})
target_link_libraries(numerics -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(numerics ${Vc_LIBRARIES})
target_link_libraries(numerics ${MPI_C_LIBRARIES})
target_link_libraries(numerics ${MPI_CXX_LIBRARIES})

if(PETSC_FOUND)
	target_include_directories (numerics PUBLIC ${PETSC_INCLUDES})
	target_link_libraries(numerics ${PETSC_LIBRARIES})
endif()

if (HIP_FOUND)
	target_link_libraries(numerics vcluster_dl)
	target_link_libraries(numerics ofpmmemory_dl)
else()
	target_link_libraries(numerics vcluster)
	target_link_libraries(numerics ofpmmemory)
endif()

if(SuiteSparse_FOUND)
	target_include_directories (numerics PUBLIC ${SuiteSparse_INCLUDE_DIRS})
	target_link_libraries(numerics ${SuiteSparse_LIBRARIES})
	target_include_directories (numerics PUBLIC ${BLAS_INCLUDE_DIRS})
	target_link_libraries(numerics ${BLAS_LIBRARIES})
endif()

if(OpenMP_CXX_FOUND)
	target_link_libraries(numerics OpenMP::OpenMP_CXX)
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
	util/SphericalHarmonics.hpp
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
	FiniteDifference/Upwind_gradient.hpp
	FiniteDifference/Eno_Weno.hpp
	FiniteDifference/FD_op.hpp
	FiniteDifference/FD_Solver.hpp
	FiniteDifference/FD_simple.hpp
	FiniteDifference/FD_expressions.hpp
	DESTINATION openfpm_numerics/include/FiniteDifference
	COMPONENT OpenFPM)

install(FILES FiniteDifference/util/common.hpp
	FiniteDifference/util/EqnsStructFD.hpp
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

install(FILES DCPSE/Dcpse.hpp
	DCPSE/Dcpse.cuh
	DCPSE/DcpseDiagonalScalingMatrix.hpp
	DCPSE/DcpseRhs.hpp
	DCPSE/Monomial.hpp
	DCPSE/Monomial.cuh
	DCPSE/MonomialBasis.hpp
	DCPSE/Support.hpp
	DCPSE/SupportBuilder.cuh
	DCPSE/SupportBuilder.hpp
	DCPSE/Vandermonde.hpp
	DCPSE/VandermondeRowBuilder.hpp
	DCPSE/DcpseInterpolation.hpp
	DESTINATION openfpm_numerics/include/DCPSE
	COMPONENT OpenFPM)

install(FILES DMatrix/EMatrix.hpp
	DESTINATION openfpm_numerics/include/DMatrix
	COMPONENT OpenFPM)

install(FILES util/eq_solve_common.hpp
	DESTINATION openfpm_numerics/include/util
	COMPONENT OpenFPM)

install(FILES DCPSE/DCPSE_op/DCPSE_op.hpp
	DCPSE/DCPSE_op/DCPSE_Solver.hpp
	DCPSE/DCPSE_op/EqnsStruct.hpp
	DCPSE/DCPSE_op/DCPSE_surface_op.hpp
	DESTINATION openfpm_numerics/include/DCPSE/DCPSE_op
	COMPONENT OpenFPM)

install(FILES OdeIntegrators/OdeIntegrators.hpp
	OdeIntegrators/vector_algebra_ofp.hpp
	OdeIntegrators/vector_algebra_ofp_gpu.hpp
	DESTINATION openfpm_numerics/include/OdeIntegrators
	COMPONENT OpenFPM)

install(FILES Draw/DrawParticles.hpp
	Draw/PointIterator.hpp
	Draw/PointIteratorSkin.hpp
	Draw/DrawDisk.hpp
	Draw/DrawSphere.hpp
	DESTINATION openfpm_numerics/include/Draw
	COMPONENT OpenFPM)

install(FILES interpolation/interpolation.hpp 
	interpolation/mp4_kernel.hpp
	interpolation/lambda_kernel.hpp
	interpolation/z_spline.hpp
	DESTINATION openfpm_numerics/include/interpolation
	COMPONENT OpenFPM)

install(FILES level_set/redistancing_Sussman/HelpFunctions.hpp
	level_set/redistancing_Sussman/HelpFunctionsForGrid.hpp
	level_set/redistancing_Sussman/NarrowBand.hpp
	level_set/redistancing_Sussman/RedistancingSussman.hpp
	level_set/redistancing_Sussman/tests/l_norms/LNorms.hpp
	level_set/redistancing_Sussman/tests/analytical_SDF/AnalyticalSDF.hpp
	DESTINATION openfpm_numerics/include/level_set/redistancing_Sussman
	COMPONENT OpenFPM)

install(FILES BoundaryConditions/MethodOfImages.hpp
	BoundaryConditions/SurfaceNormal.hpp
	DESTINATION openfpm_numerics/include/BoundaryConditions
	COMPONENT OpenFPM)

install(FILES DMatrix/EMatrix.hpp
	DESTINATION openfpm_numerics/include/DMatrix
	COMPONENT OpenFPM)

install(FILES level_set/particle_cp/particle_cp.hpp
	DESTINATION openfpm_numerics/include/level_set/particle_cp
	COMPONENT OpenFPM)

install(FILES level_set/closest_point/closest_point.hpp
	DESTINATION openfpm_numerics/include/level_set/closest_point
	COMPONENT OpenFPM)

install(FILES regression/regression.hpp
	regression/poly_levelset.hpp
	DESTINATION openfpm_numerics/include/regression
	COMPONENT OpenFPM)

install(DIRECTORY regression/minter
	DESTINATION openfpm_numerics/include/regression
	COMPONENT OpenFPM)

