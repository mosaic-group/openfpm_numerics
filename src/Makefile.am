
LINKLIBS =  $(SUITESPARSE_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)  $(METIS_LIB) $(PARMETIS_LIB) $(DEFAULT_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_IOSTREAMS_LIB) $(HDF5_LDFLAGS)  $(HDF5_LIBS) $(LIBQUADMATH)

noinst_PROGRAMS = numerics
numerics_SOURCES = main.cpp ../../openfpm_vcluster/src/VCluster.cpp ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp ../../openfpm_devices/src/Memleak_check.cpp
numerics_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CPPFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) $(SUITESPARSE_INCLUDE) $(METIS_INCLUDE) $(PARMETIS_INCLUDE)  $(EIGEN_INCLUDE) -I/usr/local/libhilbert/include -Wno-deprecated-declarations -Wno-unused-local-typedefs
numerics_CFLAGS = $(CUDA_CFLAGS)
numerics_LDADD = $(LINKLIBS) -lparmetis -lmetis
nobase_include_HEADERS = PSE/Kernels.hpp PSE/Kernels_test_util.hpp

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
