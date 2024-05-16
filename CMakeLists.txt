###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
cmake_minimum_required(VERSION 3.5)

###########################################################
project(CROCO
        VERSION 1.3.0
        DESCRIPTION "Oceanic simulation : Coastal and Regional Ocean COmmunity (CROCO)"
        LANGUAGES Fortran C CXX)

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
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros_twin_checker.cmake)

###########################################################
# Handle configuration
set(WITH_OPTIM          OFF     CACHE STRING "Select parallelisation mode : OFF, openmp, openacc-native, openacc-psyclone")
set(WITH_SPLITTING      "1x4"   CACHE STRING "Define the splitting for parallel=openmp : {NP_X}x{NP_Y}")
set(WITH_THREADS        4       CACHE STRING "Number of threads to use with parallel=openmp")
set(WITH_CASE           "BASIN" CACHE STRING "Select the case to run.")
set(WITH_KEYS           ""      CACHE STRING "List of CPPKEYS, comma separated to enable or disable (+NBQ,-IBM,FLOATS). Lands in cppdefs_edit.h")
set(WITH_TWIN_CHECKER   OFF     CACHE STRING "Enable usage of the twin checker and possibly say in which prefix to find it.")
set(ENABLE_AGRIF        OFF     CACHE BOOL   "Enable AGRIF for mesh refinement.")
set(ENABLE_TWIN_CHECKER OFF     CACHE BOOL   "Enable the twin checker to compare on the fly two croco runs made at same time.")

###########################################################
# Enable case so it lands in config.h.in
add_definitions(-DHAVE_CMAKE_CONFIG)
add_definitions(-DHAVE_CMAKE_CONFIG_CPPDEF_EDIT)
add_definitions(-DHAVE_CMAKE_CONFIG_OVERRIDE)

###########################################################
# Set default vars
set(MPI OFF)
set(OPENMP OFF)
set(OPENACC OFF)

###########################################################
# Trigger the case definition (macro variable having the name of the case)
set(${WITH_CASE} ON)
set(AGRIF ${ENABLE_AGRIF})

# Enable some vars for config_post.h.in
if (WITH_OPTIM STREQUAL "openmp")
	set(OPENMP ON)
	set(PARALLELISM_SUMMARY "openmp / threads=${WITH_THREADS}")
endif()
if (WITH_OPTIM STREQUAL "openacc-native" OR WITH_OPTIM STREQUAL "openacc-psyclone")
	set(OPENACC ON)
	set(PARALLELISM_SUMMARY ${WITH_OPTIM})
endif()
if (WITH_OPTIM STREQUAL "mpi")	
	set(MPI ON)
	set(PARALLELISM_SUMMARY "mpi / splitting=${WITH_SPLITTING}")
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
if (ENABLE_TWIN_CHECKER)
	find_package(TwinChecker REQUIRED)
endif (ENABLE_TWIN_CHECKER)
if (MPI)
	find_package(MPI REQUIRED)
endif ()
if (WITH_OPTIM STREQUAL "openacc-psyclone")
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
croco_copy_case_files()
croco_trick_create_cpp_def_override()

###########################################################
# add subdirs
add_subdirectory(AGRIF)
add_subdirectory(OCEAN)

###########################################################
# Print status to help ensure everything is correct
croco_print_status()

###########################################################
if (WITH_OPTIM STREQUAL "openacc-psyclone")
	add_test(NAME psyclone-pytest WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/PSYCLONE COMMAND pytest)
endif()
