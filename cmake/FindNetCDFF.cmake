###########################################################
# - Try to find NETCDFF (https://www.unidata.ucar.edu/software/netcdf/)
# Once done this will define
#  NETCDFF_FOUND - System has NETCDFF
#  NETCDFF_INCLUDE_DIRS - The NETCDFF include directories
#  NETCDFF_LIBRARIES - The libraries needed to use NETCDFF
#  NETCDFF_DEFINITIONS - Compiler switches required for using NETCDFF

###########################################################
set(NETCDF_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Help cmake to find NetCDF-Fortran library (https://www.unidata.ucar.edu/software/netcdf/) into your system.")

###########################################################
find_path(NETCDFF_INCLUDE_DIR
	NAMES netcdf.inc
	HINTS ${NETCDF_PREFIX}/include)

###########################################################
find_library(NETCDFF_LIBRARY NAMES netcdff
	HINTS ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)

###########################################################
set(NETCDFF_LIBRARIES ${NETCDFF_LIBRARY} )
set(NETCDFF_INCLUDE_DIRS ${NETCDFF_INCLUDE_DIR} )

###########################################################
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDFF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NetCDFF DEFAULT_MSG
	NETCDFF_LIBRARY NETCDFF_INCLUDE_DIR)

###########################################################
mark_as_advanced(NETCDFF_INCLUDE_DIR NETCDFF_LIBRARY)
