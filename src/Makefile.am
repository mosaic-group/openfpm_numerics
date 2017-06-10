
LINKLIBS = $(HDF5_LDFLAGS)  $(HDF5_LIBS) $(OPENMP_LDFLAGS) $(LIBHILBERT_LIB) $(PETSC_LIB) $(SUITESPARSE_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)  $(METIS_LIB) $(PARMETIS_LIB) $(DEFAULT_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_IOSTREAMS_LIB) $(LIBQUADMATH) $(OPENMP_LDFLAGS) $(LIBIFCORE)

noinst_PROGRAMS = numerics
numerics_SOURCES = main.cpp ../../openfpm_vcluster/src/VCluster/VCluster.cpp ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp ../../openfpm_devices/src/Memleak_check.cpp
numerics_CXXFLAGS = $(HDF5_CPPFLAGS) $(OPENMP_CFLAGS) $(LIBHILBERT_INCLUDE) $(AM_CXXFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) $(SUITESPARSE_INCLUDE) $(METIS_INCLUDE) $(PARMETIS_INCLUDE) $(EIGEN_INCLUDE) $(PETSC_INCLUDE) -Wno-deprecated-declarations -Wno-unused-local-typedefs
numerics_CFLAGS = $(CUDA_CFLAGS)
numerics_LDADD = $(LINKLIBS) -lparmetis -lmetis
nobase_include_HEADERS = Matrix/SparseMatrix.hpp Matrix/SparseMatrix_Eigen.hpp Matrix/SparseMatrix_petsc.hpp \
Vector/Vector_eigen.hpp Vector/Vector_petsc.hpp Vector/Vector_util.hpp Vector/Vector.hpp \
Solvers/umfpack_solver.hpp Solvers/petsc_solver.hpp \
util/petsc_util.hpp util/linalgebra_lib.hpp util/util_num.hpp util/grid_dist_testing.hpp  \
FiniteDifference/Average.hpp FiniteDifference/Derivative.hpp FiniteDifference/eq.hpp FiniteDifference/FDScheme.hpp FiniteDifference/Laplacian.hpp FiniteDifference/mul.hpp FiniteDifference/sum.hpp FiniteDifference/util/common.hpp \
PSE/Kernels.hpp PSE/Kernels_test_util.hpp Operators/Vector/vector_dist_operators_extensions.hpp Operators/Vector/vector_dist_operators.hpp Operators/Vector/vector_dist_operators_apply_kernel.hpp Operators/Vector/vector_dist_operators_functions.hpp Operators/Vector/vector_dist_operator_assign.hpp \
Draw/DrawParticles.hpp Draw/PointIterator.hpp Draw/PointIteratorSkin.hpp \
interpolation/interpolation.hpp interpolation/mp4_kernel.hpp interpolation/z_spline.hpp

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

test: numerics
	source $(HOME)/openfpm_vars && cd .. && mpirun -np 2 ./src/numerics && mpirun -np 3 ./src/numerics && mpirun -np 4 ./src/numerics

