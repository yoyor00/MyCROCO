###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
cmake_minimum_required(VERSION 3.1)

###########################################################
project(CROCO
        VERSION 1.3.0
        DESCRIPTION "Oceanic simulation : Coastal and Regional Ocean COmmunity (CROCO)"
        LANGUAGES Fortran C)

###########################################################
# to support old versions
cmake_policy(SET CMP0057 NEW)
# enable testing
enable_testing()

###########################################################
#Add internal search path for scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/)

###########################################################
#Delegate some cmake scripts
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros_psyclone.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros_agrif.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/croco_compiler_options.cmake)

###########################################################
# Handle configuration
set(WITH_PARALLEL  OFF     CACHE STRING "Select parallelisation mode : OFF, openmp, openacc-native, openacc-psyclone")
set(WITH_SPLITTING "1x4"   CACHE STRING "Define the splitting for parallel=openmp : {NP_X}x{NP_Y}")
set(WITH_THREADS   4       CACHE STRING "Number of threads to use with parallel=openmp")
set(WITH_CASE      "BASIN" CACHE STRING "Select the case to run.")
set(ENABLE_AGRIF   OFF     CACHE BOOL   "Enable AGRIF for mesh refinement.")

###########################################################
# Enable case so it lands in config.h.in
add_definitions(-DHAVE_CMAKE_CONFIG)

# Trigger the case definition (macro variable having the name of the case)
set(${WITH_CASE} ON)
set(AGRIF ${ENABLE_AGRIF})

# Enable some vars for config_post.h.in
if (WITH_PARALLEL STREQUAL "openmp")
	set(OPENMP ON)
endif()
if (WITH_PARALLEL STREQUAL "openacc-native" OR WITH_PARALLEL STREQUAL "openacc-psyclone")
	set(OPENACC ON)
endif()
if (WITH_PARALLEL STREQUAL "mpi")	
	set(MPI ON)
endif ()

# enable some static config
set(NB_THREADS ${WITH_THREADS})
# extract NP_x & NP_eta
string(REPLACE "x" ";" WITH_SPLITTING ${WITH_SPLITTING})
list(GET WITH_SPLITTING 0 SPLITTING_X)
list(GET WITH_SPLITTING 1 SPLITTING_ETA)

###########################################################
# Search dependences
find_package(NetCDFF REQUIRED)
if (MPI)
	find_package(MPI REQUIRED)
endif ()
if (WITH_PARALLEL STREQUAL "openacc-psyclone")
	find_package(PSyClone REQUIRED)
endif ()
if (OPENMP)
	find_package(OpenMP REQUIRED)
endif ()

###########################################################
# Compiler build flags
croco_tune_compile_flags()

###########################################################
# Check
croco_last_checkings()

###########################################################
# some croco trick to keep same way to work than with jobcomp
croco_trick_create_cpp_def_override()
include_directories(${CMAKE_BINARY_DIR})

###########################################################
# add subdirs
add_subdirectory(AGRIF)
add_subdirectory(OCEAN)

###########################################################
# Print status to help ensure everything is correct
croco_print_status()

###########################################################
if (WITH_PARALLEL STREQUAL "openacc-psyclone")
	add_test(NAME psyclone-pytest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/PSYCLONE COMMAND pytest)
endif()
