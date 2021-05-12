cmake_minimum_required(VERSION 3.13)
project(Boost_Tests_run)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER_LAUNCHER="/usr/bin/ccache")
set(DCMAKE_C_COMPILER_LAUNCHER="/usr/bin/ccache")


set (BOOST_ROOT "/Users/jstark/OPENFPM_DEPENDENCIES_dev/BOOST")
set (Boost_NO_BOOST_CMAKE ON CACHE BOOL "" FORCE)

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)


# For unit testing
set (Boost_USE_STATIC_LIBS OFF)


# Find MPI
find_package(MPI REQUIRED)
# Need unit_test_framework for unit testing
find_package(Boost REQUIRED COMPONENTS unit_test_framework filesystem system log log_setup thread)

link_directories(/Users/jstark/OpenFPM_install_dev/openfpm_devices/lib
		/Users/jstark/OpenFPM_install_dev/openfpm_pdata/lib  /Users/jstark/OpenFPM_install_dev/openfpm_vcluster/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/VCDEVEL/lib  /Users/jstark/OPENFPM_DEPENDENCIES_dev/METIS/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/PARMETIS/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/BOOST/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/HDF5/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/LIBHILBERT/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/PETSC/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/OPENBLAS/lib /Users/jstark/OPENFPM_DEPENDENCIES_dev/SUITESPARSE/lib )


add_executable (Boost_Tests_run redistancingSussman_unit_test.cpp computeGradient_unit_test.cpp)

# Supress warning
target_compile_options(Boost_Tests_run PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wno-macro-redefined >)

# Link against MPI
target_link_libraries(Boost_Tests_run ${MPI_LIBRARIES} ${Boost_LIBRARIES} vcluster ofpm_pdata ofpmmemory parmetis metis boost_iostreams boost_program_options hdf5 libhilbert Vc petsc openblas umfpack amd btf camd ccolamd cholmod colamd cxsparse klu dl rbio spqr suitesparseconfig dl boost_filesystem boost_system boost_filesystem boost_system)

target_include_directories(Boost_Tests_run PUBLIC /. /Users/jstark/OpenFPM_install_dev/openfpm_numerics/include /Users/jstark/OpenFPM_install_dev/openfpm_pdata/include/config /Users/jstark/OpenFPM_install_dev/openfpm_pdata/include /Users/jstark/OpenFPM_install_dev/openfpm_data/include /Users/jstark/OpenFPM_install_dev/openfpm_vcluster/include /Users/jstark/OpenFPM_install_dev/openfpm_io/include /Users/jstark/OpenFPM_install_dev/openfpm_devices/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/VCDEVEL/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/METIS/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/PARMETIS/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/BOOST/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/HDF5/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/LIBHILBERT/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/PETSC/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/OPENBLAS/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/SUITESPARSE/include /Users/jstark/OPENFPM_DEPENDENCIES_dev/EIGEN)


# Include MPI and BOOST includes
include_directories(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})
if(MPI_COMPILE_FLAGS)
	set_target_properties(Boost_Tests_run PROPERTIES
			COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
	set_target_properties(Boost_Tests_run PROPERTIES
			LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
