###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Loop on all values of the list and make them absolute path
#
# Parameters
# ----------
# list_to_update:
#     Name of the list to loop in and to update.
# path:
#     Path to preprend to each sub paths.
function(croco_make_absolute_paths list_to_update path)
	# reset the list
	set(_newfiles)

	# loop on all files
	foreach(file IN LISTS ${list_to_update})
		list(APPEND _newfiles ${path}/${file})
	endforeach()

	# export
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()

###########################################################
# Apply the CPP pre-processor and MPC code reshaper
# to prepare the sources before build.
#
# Parameters
# ----------
# list_to_update : list [IN,OUT]
#     Name of the list to loop in and in which to replace all the file
#     names after building the make rules.
#
# Inspiration
# -----------
# https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
function(croco_cpp_and_mpc_preprocess list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# build abs path
		get_filename_component(oldfile_name ${oldfile} NAME)
		set(oldfile_abs ${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${oldfile_name})

		# build new name
		string(REGEX REPLACE "\\.F" ".cpp.F"   newfile_cpp           ${oldfile_abs})
		string(REGEX REPLACE "\\.F" ".mpc.F"   newfile_cpp_mpc       ${newfile_cpp})

		# build command
		add_custom_command(
			OUTPUT ${newfile_cpp_mpc}
			COMMAND ${CROCO_FORTRAN_CPP} ${CROCO_FORTRAN_CPP_FLAGS} ${oldfile} -o ${newfile_cpp}
			COMMAND ${CROCO_MPC} < ${newfile_cpp} > ${newfile_cpp_mpc}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc ${CROCO_CPP_H}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile_cpp_mpc})
	endforeach()

	# export
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()

###########################################################
# Print a summary status to help ensuring everything
# is correct
function(croco_print_status)
	# Prepare some vars for printing
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
	message(STATUS "|  NetCDF           : ${NETCDFF_LIBRARY}")
	message(STATUS "|  Parallelism      : ${WITH_PARALLEL}")
	message(STATUS "|  PSyClone         : ${PSYCLONE_COMMAND}")
	message(STATUS "|  OpenACC          : ${OPENACC}")
	message(STATUS "|  MPI              : ${MPI_FOUND}")
	message(STATUS "|  OpenMP           : ${ENABLE_OPENMP}")
	message(STATUS "|  AGRIF            : ${ENABLE_AGRIF}")
	message(STATUS "--------------------------------------------------------------")
	message(STATUS "|  CMake build type : ${CMAKE_BUILD_TYPE}")
	message(STATUS "|  User fflags      : ${CMAKE_Fortran_FLAGS}")
	message(STATUS "|  CMake fflags     : ${CMAKE_BUILD_TYPE_FLAGS}")
	message(STATUS "|  Croco fflags     : ${CROCO_FORTRAN_FLAGS_STR}")
	message(STATUS "--------------------------------------------------------------")
	message(STATUS "|  Case             : ${WITH_CASE}")
	message(STATUS "==============================================================")
endfunction()

###########################################################
# Perform some extra checks on variables to see if everything is correct
function(croco_last_checkings)
	# allowed
	list(APPEND para_allowed OFF openmp openacc-native openacc-psyclone)
	if (NOT ${WITH_PARALLEL} IN_LIST para_allowed)
		message(FATAL_ERROR "Select an invalid parallelism mode : -DWITH_PARALLEL=${WITH_PARALLEL}, should be in (${para_allowed})")
	endif()

	# Require psyclone
	if (WITH_PARALLEL STREQUAL "openacc-psyclone" AND NOT PSYCLONE_FOUND)
		message(FATAL_ERROR "Fail to find PSyClone, required if enabling PSyClone OpenACC via -DWITH_PSYCLONE_VENV !")
	endif ()
endfunction()
