#define BOOST_TEST_MODULE SprngEngine_test
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <boost/math/distributions/chi_squared.hpp>
#include <random>

#include "random/SprngEngine.hpp"

using BigN = unsigned long long int;
using int2int = std::map<BigN, BigN>;

BOOST_AUTO_TEST_SUITE(SprngEngine_test)

template <typename RandomNumberEngine>
auto probPi(RandomNumberEngine& gen, const BigN n) -> pointnumber {
  std::uniform_real_distribution<double> dist(-1, 1);

  BigN nInside = 0;
  for (auto i = 0; i < n; ++i) {
    pointnumber x = dist(gen);
    pointnumber y = dist(gen);

    if ((pow(x, 2) + pow(y, 2)) < 1) {
      nInside += 1;
    }
  }

  pointnumber p = (pointnumber)nInside / n;
  return p;
}

auto errStdPi(const pointnumber p, const BigN n) -> pointnumber {
  return 2.0 * sqrt(p * (1 - p) / n);
}

auto errPi(const pointnumber p) -> pointnumber { return abs(4.0 * p - M_PI); }

template <typename RandomNumberEngine>
auto rollDice(RandomNumberEngine& gen, const BigN nSides) -> BigN {
  std::uniform_int_distribution<BigN> dist(1, nSides);
  return dist(gen);
}

template <typename RandomNumberEngine>
auto rollDice(RandomNumberEngine& gen, const BigN nSides, const BigN n)
    -> std::map<BigN, BigN> {
  std::map<BigN, BigN> bins{};  // faces: # rolls

  for (auto i = 0; i < n; ++i) {
    ++bins[rollDice(gen, nSides)];
  }

  return bins;
}

BOOST_AUTO_TEST_CASE(SprngEngine_die_test) {
  std::random_device rd;
  SprngEngine<std::uint_fast32_t, 1> gen(rd());

  constexpr BigN nSides = 6;
  constexpr BigN n = 1e6;
  auto rolls = rollDice(gen, nSides, n);

  const double expectedRollPerSide = (double)n / nSides;
  const double chiSquaredResult = std::accumulate(
      std::begin(rolls),
      std::end(rolls),
      0.0,
      [&expectedRollPerSide](double value, const int2int::value_type& p) {
        return value + pow((double)p.second - expectedRollPerSide, 2) /
                           expectedRollPerSide;
      });

  const BigN degFreedom = nSides - 1;
  boost::math::chi_squared dist(degFreedom);

  const double alpha = 0.05;
  const double chiSquareTest = boost::math::quantile(dist, alpha);

  BOOST_CHECK(chiSquaredResult < chiSquareTest);
}

BOOST_AUTO_TEST_CASE(SprngEngine_pi_test) {
  std::random_device rd;
  SprngEngine<std::uint_fast32_t, 2> gen(rd());

  constexpr BigN n = 1e7;
  pointnumber p = probPi(gen, n);
  pointnumber stdError = errStdPi(p, n);
  pointnumber error = errPi(p);

  BOOST_CHECK(error < stdError);  // error is within std error
}

BOOST_AUTO_TEST_CASE(SprngEngine_VS_STD) {
  std::random_device rd, rd1;
  SprngEngine<std::uint_fast32_t, 3> genSPRNG(rd());
  std::default_random_engine genSTD(rd1());

  BigN betterThanSTD = 0, matches = 0;

  for (auto n = 1e1; n < 1e7; n *= 5) {
    pointnumber pSPRNG = probPi(genSPRNG, n);
    pointnumber errorSPRNG = errPi(pSPRNG);

    pointnumber pSTD = probPi(genSTD, n);
    pointnumber errorSTD = errPi(pSTD);

    if (errorSPRNG < errorSTD) {
      betterThanSTD += 1;
    }
    matches += 1;
  }

  BOOST_CHECK(betterThanSTD > 0.5 * matches);
}

BOOST_AUTO_TEST_SUITE_END()
