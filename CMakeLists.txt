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
include(${CMAKE_SOURCE_DIR}/cmake/croco_compiler_options.cmake)

######################################################
# Handle configuration
set(ENABLE_PSYCLONE_OPENACC OFF CACHE BOOL "Enable usage of PSyClone to genrate the GPU/OpenACC version of the code.")
set(CROCO_CASE "BASIN" CACHE STRING "Select the case to run.")
set(ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP for CPU parallelism.")

######################################################
# Enable case so it lands in config.h.in
add_definitions(-DHAS_CMAKE_CONFIG)
# Trigger the case definition (macro variable having the name of the case)
set(${CROCO_CASE} ON)
set(OPENACC ${ENABLE_PSYCLONE_OPENACC})
set(OPENMP ${ENABLE_OPENMP})

######################################################
# Search dependences
find_package(NetCDFF REQUIRED)
if (ENABLE_PSYCLONE_OPENACC)
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
if (ENABLE_PSYCLONE_OPENACC AND NOT PSYCLONE_FOUND)
	message(FATAL_ERROR "Fail to find PSyClone, required if enabling PSyClone OpenACC via -DWITH_PSYCLONE_VENV !")
endif (ENABLE_PSYCLONE_OPENACC AND NOT PSYCLONE_FOUND)

######################################################
# add subdirs
add_subdirectory(AGRIF)
add_subdirectory(OCEAN)

######################################################
# Prepare some vars
list(JOIN CROCO_FORTRAN_FLAGS " " CROCO_FORTRAN_FLAGS_STR)
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
if (CMAKE_BUILD_TYPE STREQUAL "")
	list(JOIN CMAKE_Fortran_FLAGS " " CMAKE_BUILD_TYPE_FLAGS)
else()
	list(JOIN CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPPER} " " CMAKE_BUILD_TYPE_FLAGS)
endif()

# Print summary
message(STATUS "==============================================================")
message(STATUS "|  OS               : ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "|  Compiler familly : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "|  Compiler         : ${CMAKE_Fortran_COMPILER}")
message(STATUS "--------------------------------------------------------------")
message(STATUS "|  NetCDF           : ${NETCDFF_FOUND}")
message(STATUS "|  PsyClone         : ${PSYCLONE_FOUND}")
message(STATUS "|  OpenACC          : ${OPENACC}")
message(STATUS "|  MPI              : ${MPI_FOUND}")
message(STATUS "|  OpenMP           : ${OPENMP_FOUND}")
message(STATUS "--------------------------------------------------------------")
message(STATUS "|  CMake build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "|  User fflags      : ${CMAKE_Fortran_FLAGS}")
message(STATUS "|  CMake fflags     : ${CMAKE_BUILD_TYPE_FLAGS}")
message(STATUS "|  Croco fflags     : ${CROCO_FORTRAN_FLAGS_STR}")
message(STATUS "--------------------------------------------------------------")
message(STATUS "|  Case             : ${CROCO_CASE}")
message(STATUS "==============================================================")
