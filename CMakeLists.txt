######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
cmake_minimum_required(VERSION 3.1)

######################################################
project(CROCO
        VERSION 1.3.0
        DESCRIPTION "Oceanic simulation : Coastal and Regional Ocean COmmunity (CROCO)"
        LANGUAGES Fortran C)

######################################################
#Add internal search path for scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/)

######################################################
#Delegate some cmake scripts
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros_psyclone.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/croco_macros_agrif.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/croco_compiler_options.cmake)

######################################################
# Handle configuration
set(WITH_OPENACC   OFF    CACHE STRING "Enable OenACC version, either by using 'native' or 'psyclone'.")
set(CROCO_CASE    "BASIN" CACHE STRING "Select the case to run.")
set(ENABLE_OPENMP OFF     CACHE BOOL   "Enable OpenMP for CPU parallelism.")
set(ENABLE_AGRIF  OFF     CACHE BOOL   "Enable AGRIF for mesh refinement.")

######################################################
# Enable case so it lands in config.h.in
add_definitions(-DHAVE_CMAKE_CONFIG)
# Trigger the case definition (macro variable having the name of the case)
set(${CROCO_CASE} ON)
# Enable some vars for config.parallelism.h.in
set(OPENMP ${ENABLE_OPENMP})
set(AGRIF ${ENABLE_AGRIF})
if (WITH_OPENACC)
	set(OPENACC ON)
endif()

######################################################
# Search dependences
find_package(NetCDFF REQUIRED)
if (WITH_OPENACC STREQUAL "psyclone")
	find_package(PSyClone REQUIRED)
endif ()
if (ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
endif ()

######################################################
# Compiler build flags
croco_tune_compile_flags()

######################################################
# Check
croco_last_checkings()

######################################################
# add subdirs
add_subdirectory(AGRIF)
add_subdirectory(OCEAN)

######################################################
# Print status to help ensure everything is correct
croco_print_status()
