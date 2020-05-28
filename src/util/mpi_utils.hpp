#ifndef OPENFPM_NUMERICS_SRC_UTIL_MPI_UTILS_HPP_
#define OPENFPM_NUMERICS_SRC_UTIL_MPI_UTILS_HPP_

#include <mpi.h>

auto getRank() -> int {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,
                &rank);  // store rank of MPI process inside communicator
  return rank;
}

auto getSize() -> int {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,
                &size);  // store size MPI processes inside communicator
  return size;
}

auto amIMaster() -> bool { return getRank() == 0; }

#endif /* OPENFPM_NUMERICS_SRC_UTIL_MPI_UTILS_HPP_ */
