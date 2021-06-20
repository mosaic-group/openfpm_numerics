#include "initialize_wrapper.hpp"
#include "VCluster/VCluster.hpp"

#ifndef NO_INIT_AND_MAIN

void openfpm_init_wrapper(int * argc, char *** argv)
{
	openfpm_init(argc,argv);
}

void openfpm_finalize_wrapper()
{
	openfpm_finalize();
}

#endif