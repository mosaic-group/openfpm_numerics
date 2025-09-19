# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter"
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-build"
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter"
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/tmp"
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp"
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src"
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp${cfgdir}") # cfgdir has leading slash
endif()
