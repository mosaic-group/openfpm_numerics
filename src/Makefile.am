
LINKLIBS =  $(SUITESPARSE_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)  $(METIS_LIB) $(DEFAULT_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_IOSTREAMS_LIB)

noinst_PROGRAMS = numerics
numerics_SOURCES = main.cpp ../../openfpm_vcluster/src/VCluster.cpp ../../openfpm_devices/src/memory/HeapMemory.cpp ../../openfpm_devices/src/memory/PtrMemory.cpp ../../openfpm_devices/src/Memleak_check.cpp
numerics_CXXFLAGS = $(INCLUDES_PATH) $(BOOST_CPPFLAGS) $(SUITESPARSE_INCLUDE) $(METIS_INCLUDE) $(EIGEN_INCLUDE) -Wno-deprecated-declarations -Wunused-local-typedefs
numerics_CFLAGS = $(CUDA_CFLAGS)
numerics_LDADD = $(LINKLIBS) -lmetis

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
