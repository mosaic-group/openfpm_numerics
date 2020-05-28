#ifndef OPENFPM_NUMERICS_SRC_RANDOM_SPRNG_ENGINE_HPP_
#define OPENFPM_NUMERICS_SRC_RANDOM_SPRNG_ENGINE_HPP_

#include <fstream>
#include <random>

#ifdef FLOAT_GEN
using pointnumber = float;
#else
using pointnumber = double;
#endif

#include "sprng_cpp.h"

#ifdef USE_MPI
#include "util/mpi_utils.hpp"
#endif

/**
 * Adheres to https://en.cppreference.com/w/cpp/named_req/RandomNumberEngine
 */
template <typename UIntType, UIntType gen>
class SprngEngine {
public:
  // NOLINTNEXTLINE
  typedef UIntType result_type;  // the type of the generated random value

  static constexpr UIntType GARBAGE_SEED = 0;
  static constexpr UIntType SPRNG_GEN_TYPE = gen;
  static constexpr UIntType DEFAULT_STREAM_ID = 0;
  static constexpr UIntType DEFAULT_NUM_STREAMS = 1;
  static constexpr int DEFAULT_PARAM = SPRNG_DEFAULT;

#ifdef USE_MPI
  SprngEngine() : SprngEngine(makeSeed(), DEFAULT_PARAM) {}

  template <typename Sseq>
  explicit SprngEngine(Sseq& q) : SprngEngine((UIntType)q, DEFAULT_PARAM) {}

  explicit SprngEngine(UIntType seed) : SprngEngine(seed, DEFAULT_PARAM) {}

  SprngEngine(const UIntType seed, const int param)
    : SprngEngine(seed, getRank(), getSize(), param) {}
#else
  SprngEngine() : SprngEngine(makeSeed()) {}

  template <typename Sseq>
  explicit SprngEngine(Sseq& q) : SprngEngine((UIntType)q) {}

  explicit SprngEngine(UIntType seed)
    : SprngEngine(seed, DEFAULT_STREAM_ID, DEFAULT_NUM_STREAMS, DEFAULT_PARAM) {
  }
#endif

  SprngEngine(const UIntType seed,
              const unsigned int streamId,
              const unsigned int nStreams,
              const int param) {
    _s = SelectType(gen);
    init(streamId, nStreams, seed, param);
  }

  SprngEngine(char* buffer) {
    _s = SelectType(gen);
    unPack(buffer);
  }

  ~SprngEngine() {
    if (getCount() > 0 && _s) {
      getCounter() = _s->free_sprng();
    }
  }

  static constexpr auto min() -> UIntType { return 0; }

  static constexpr auto max() -> UIntType { return (1ULL << 31) - 1; }

  void discard(const unsigned long long z) {
    for (auto i = 0; i < z; ++i) {
      (*this)();
    }
  }

  auto operator()() -> UIntType { return randInt(); }

  /**
   * @brief Compares two linear congruential random number generator
   * objects of the same type for equality.
   */
  auto operator==(const SprngEngine& other) const -> bool {
    // will not check streams left, it is inferred by state

    char* myBuffer;
    char* otherBuffer;
    const int mySize = pack(&myBuffer);
    const int otherSize = other.pack(&otherBuffer);
    bool out = false;

    if (mySize == otherSize) {
      for (auto i = 0; i < mySize; ++i) {
        if (myBuffer[i] != otherBuffer[i]) {
          free(myBuffer);
          free(otherBuffer);
          return false;
        }
      }

      out = true;
    }

    free(myBuffer);
    free(otherBuffer);
    return out;
  }

  /**
   * @brief Compares two linear congruential random number generator
   * objects of the same type for inequality.
   */
  auto operator!=(const SprngEngine& other) const -> bool {
    return !(*this == other);
  }

  /**
   * @brief Writes the textual representation
   */
  template <typename _CharT, typename _Traits>
  friend auto operator<<(std::basic_ostream<_CharT, _Traits>& os,
                         const SprngEngine& x)
      -> std::basic_ostream<_CharT, _Traits>& {
    x.printInfo();  // todo

    os << "minimum: " << min() << ", maxium: " << max();
    return os;
  }

  /**
   * @brief Sets the state of the engine by reading its textual representation
   */
  template <typename _CharT, typename _Traits>
  friend auto operator>>(std::basic_istream<_CharT, _Traits>& is,
                         const SprngEngine& x)
      -> std::basic_istream<_CharT, _Traits>& {
    // out.unPack(buffer); // todo
    return is;
  }

  static auto makeSeed() -> UIntType { return make_sprng_seed(); }

  static auto getCount() -> std::size_t { return getCounter(); }

  auto toFile(std::string fileName) const -> int {
    FILE* fp = fopen(fileName.c_str(), "w");
    if (fp) {
      char* bytes;
      const int size = pack(&bytes);
      fwrite(&size, 1, sizeof(int), fp);
      fwrite(bytes, 1, size, fp);

      fclose(fp);
      free(bytes);
      return size;
    }

    return 0;
  }

  static auto fromString(char* buffer) -> SprngEngine {
    return SprngEngine(buffer);
  }

  auto fromFile(const std::string fileName) -> SprngEngine {
    char buffer[MAX_PACKED_LENGTH];
    FILE* fp = fopen(fileName.c_str(), "r");
    unsigned int numBytes;

    fseek(fp, 0L, SEEK_END);  // get # bytes
    numBytes = ftell(fp);

    fseek(fp, 0L, SEEK_SET);  // reset file position indicator

    fread(buffer, sizeof(char), numBytes, fp);  // read file
    fclose(fp);

    return SprngEngine(buffer);
  }

  void printInfo() const { _s->print_sprng(); }

  auto pack(char** buffer) const -> size_t { return _s->pack_sprng(buffer); }

  auto unPack(char* buffer) const -> size_t { return _s->unpack_sprng(buffer); }

private:
  void init(const int streamId,
            const int nStreams,
            const int seed,
            const int param) {
    _s->init_sprng(streamId, nStreams, seed, param);
    getCounter() += 1;
  }

  auto randInt() const -> UIntType { return _s->isprng(); }

  static auto getCounter() -> std::size_t& {
    static std::size_t counter = 0;
    return counter;
  }

  Sprng* _s = nullptr;
};

#endif /* OPENFPM_NUMERICS_SRC_RANDOM_SPRNG_ENGINE_HPP_ */
