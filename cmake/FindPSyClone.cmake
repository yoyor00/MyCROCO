######################################################
# - Try to find PSyClone (https://psyclone.readthedocs.io/)
# Once done this will define
#  PSYCLONE_FOUND   - System has PSYCLONE
#  PSYCLONE_COMMAND - Path to the psyclone command
#  PSYCLONE_VENV    - Path to the psyclone venv

######################################################
set(PSYCLONE_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Help cmake to find PSyClone (https://psyclone.readthedocs.io/) into your system.")

######################################################
find_program(PSYCLONE_COMMAND 
	NAMES psyclone
	HINTS ${PSYCLONE_PREFIX}/bin)

######################################################
find_path(PSYCLONE_VENV 
	NAMES activate
	HINTS ${PSYCLONE_PREFIX}/bin)

######################################################
get_filename_component(PSYCLONE_VENV ${PSYCLONE_VENV} DIRECTORY)

######################################################
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PSYCLONE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PSyClone DEFAULT_MSG
	PSYCLONE_COMMAND PSYCLONE_VENV)

######################################################
mark_as_advanced(PSYCLONE_COMMAND PSYCLONE_VENV)