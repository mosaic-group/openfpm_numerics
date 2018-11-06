#include "config.h"
#undef VERSION
#define BOOST_TEST_MODULE "C++ test module for OpenFPM_numerics project"
#include <boost/test/included/unit_test.hpp>

#include "unit_test_init_cleanup.hpp"

#include "config.h"
#include "FiniteDifference/util/common_test.hpp"
#include "util/util_num_unit_tests.hpp"
#include "PSE/Kernels_unit_tests.hpp"
#include "Draw/DrawParticles_unit_tests.hpp"
