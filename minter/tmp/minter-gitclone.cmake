# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitclone-lastrun.txt" AND EXISTS "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitinfo.txt" AND
  "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitclone-lastrun.txt" IS_NEWER_THAN "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitinfo.txt")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitclone-lastrun.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/Users/foggia/Homebrew/bin/git"
            clone --no-checkout --config "advice.detachedHead=false" "https://git.mpi-cbg.de/mosaic/software/math/minter.git" "minter"
    WORKING_DIRECTORY "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://git.mpi-cbg.de/mosaic/software/math/minter.git'")
endif()

execute_process(
  COMMAND "/Users/foggia/Homebrew/bin/git"
          checkout "origin/header_only" --
  WORKING_DIRECTORY "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'origin/header_only'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/Users/foggia/Homebrew/bin/git" 
            submodule update --recursive --init 
    WORKING_DIRECTORY "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitinfo.txt" "/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/Users/foggia/opt/extra/openfpm/openfpm_numerics/minter/src/minter-stamp/minter-gitclone-lastrun.txt'")
endif()
