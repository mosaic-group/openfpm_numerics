#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "random/SprngEngine.hpp"

BOOST_AUTO_TEST_SUITE (SprngEngine_test)

int f() {
  return 0;
}

BOOST_AUTO_TEST_CASE(SprngEngine_normal_distribution_test) {
  BOOST_REQUIRE_SMALL(f(), 0.001f);
}

BOOST_AUTO_TEST_CASE(SprngEngine_pi_test) {
  BOOST_REQUIRE_SMALL(f(), 0.001f);
}

BOOST_AUTO_TEST_SUITE_END()
