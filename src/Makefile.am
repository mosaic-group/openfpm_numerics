
LINKLIBS = $(DEFAULT_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS)

noinst_PROGRAMS = numerics
numerics_SOURCES = main.cpp
numerics_CXXFLAGS = $(INCLUDES_PATH) $(BOOST_CPPFLAGS)
numerics_CFLAGS = $(CUDA_CFLAGS)
numerics_LDADD = $(LINKLIBS)

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
