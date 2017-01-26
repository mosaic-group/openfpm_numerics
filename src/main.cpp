#define BOOST_TEST_MODULE "C++ test module for OpenFPM_numerics project"
#include <boost/test/included/unit_test.hpp>

#include "unit_test_init_cleanup.hpp"

#include "config.h"
#include "Vector/Vector_unit_tests.hpp"
#include "Matrix/SparseMatrix_unit_tests.hpp"
#include "FiniteDifference/FDScheme_unit_tests.hpp"
#include "FiniteDifference/util/common_test.hpp"
#include "FiniteDifference/eq_unit_test.hpp"
#include "FiniteDifference/eq_unit_test_3d.hpp"
#include "util/util_num_unit_tests.hpp"
#include "PSE/Kernels_unit_tests.hpp"
#include "Operators/Vector/vector_dist_operators_unit_tests.hpp"
#include "Draw/DrawParticles_unit_tests.hpp"
