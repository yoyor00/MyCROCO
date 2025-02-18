###########################################################
# - Try to find TWIN_CHECKER (https://gitlab.inria.fr/svalat/twin-checker/)
# Once done this will define
#  TWIN_CHECKER_FOUND - System has TWIN_CHECKER
#  TWIN_CHECKER_INCLUDE_DIRS - The TWIN_CHECKER include directories
#  TWIN_CHECKER_LIBRARIES - The libraries needed to use TWIN_CHECKER
#  TWIN_CHECKER_DEFINITIONS - Compiler switches required for using TWIN_CHECKER

###########################################################
set(TWIN_CHECKER_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Help cmake to find twin-checker tool (https://gitlab.inria.fr/svalat/twin-checker/) into your system.")

###########################################################
find_program(TWIN_CHECKER_RUN
	NAMES twin-checker-run
	HINTS ${TWIN_CHECKER_PREFIX}/bin
)

###########################################################
find_program(TWIN_CHECKER_PSYCLONE_INSTRU
	NAMES twin-checker-psyclone-instru
	HINTS ${TWIN_CHECKER_PREFIX}/bin
)

###########################################################
find_program(TWIN_CHECKER_FILE_LIST
	NAMES twin-checker-file-list
	HINTS ${TWIN_CHECKER_PREFIX}/bin
)

###########################################################
find_program(TWIN_CHECKER_COMPILER_WRAPPER
	NAMES twin-checker-compiler-wrap
	HINTS ${TWIN_CHECKER_PREFIX}/bin
)

###########################################################
get_filename_component(bin_path ${TWIN_CHECKER_PSYCLONE_INSTRU} DIRECTORY)
get_filename_component(TWIN_CHECKER_BIN_PREFIX ${bin_path} DIRECTORY)

###########################################################
message(STATUS "Found twin-checker - ${TWIN_CHECKER_BIN_PREFIX}")
message(STATUS "Found twin-checker-psyclone-instru - ${TWIN_CHECKER_PSYCLONE_INSTRU}")

###########################################################
find_path(TWIN_CHECKER_INCLUDE_DIR
	NAMES twin-checker/CheckerApi.h
	HINTS ${TWIN_CHECKER_BIN_PREFIX}/include
          ${TWIN_CHECKER_PREFIX}/include)

###########################################################
message(STATUS "Found twin-checker-include - ${TWIN_CHECKER_INCLUDE_DIR}")

###########################################################
find_library(TWIN_CHECKER_LIBRARY NAMES twin-checker
	HINTS ${TWIN_CHECKER_BIN_PREFIX}/lib
	      ${TWIN_CHECKER_BIN_PREFIX}/lib64
	      ${TWIN_CHECKER_PREFIX}/lib
	      ${TWIN_CHECKER_PREFIX}/lib64)

###########################################################
set(TWIN_CHECKER_LIBRARIES ${TWIN_CHECKER_LIBRARY} )
set(TWIN_CHECKER_INCLUDE_DIRS ${TWIN_CHECKER_INCLUDE_DIR} )

###########################################################
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TWIN_CHECKER_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TwinChecker DEFAULT_MSG
    TWIN_CHECKER_LIBRARY
    TWIN_CHECKER_INCLUDE_DIR
    TWIN_CHECKER_RUN
    TWIN_CHECKER_PSYCLONE_INSTRU
    TWIN_CHECKER_FILE_LIST
    TWIN_CHECKER_COMPILER_WRAPPER
    TWIN_CHECKER_BIN_PREFIX
)

###########################################################
mark_as_advanced(TWIN_CHECKER_INCLUDE_DIR
                 TWIN_CHECKER_LIBRARY
                 TWIN_CHECKER_RUN
                 TWIN_CHECKER_PSYCLONE_INSTRU
                 TWIN_CHECKER_FILE_LIST
                 TWIN_CHECKER_COMPILER_WRAPPER
)
