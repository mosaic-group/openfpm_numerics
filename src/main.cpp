#include "config.h"
#undef VERSION
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "unit_test_init_cleanup.hpp"

#include "config.h"
#include "FiniteDifference/util/common_test.hpp"
#include "util/util_num_unit_tests.hpp"
#include "PSE/Kernels_unit_tests.hpp"
#include "Draw/DrawParticles_unit_tests.hpp"

// initialization function:
bool init_unit_test()
{
  return true;
}

// entry point:
int main(int argc, char* argv[])
{
  return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}
