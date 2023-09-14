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

######################################################
# Handle configuration
set(ENABLE_PSYCLONE_OPENACC OFF CACHE BOOL "Enable usage of PSyClone to genrate the GPU/OpenACC version of the code.")
set(WITH_CASE "BASIN" CACHE STRING "Select the case to run.")

######################################################
# Enable case so it lands in config.h.in
add_definitions(-DHAS_CMAKE_CONFIG)
# Trigger the case definition (macro variable having the name of the case)
set(${WITH_CASE} ON)

######################################################
# Search dependences
find_package(NetCDFF REQUIRED)
find_package(PSyClone QUIET)

######################################################
# Check
if (ENABLE_PSYCLONE_OPENACC AND NOT PSYCLONE_FOUND)
	message(FATAL_ERROR "Fail to find PSyClone, required if enabling PSyClone OpenACC via -DWITH_PSYCLONE_VENV !")
endif (ENABLE_PSYCLONE_OPENACC AND NOT PSYCLONE_FOUND)

######################################################
# add subdirs
add_subdirectory(OCEAN)
