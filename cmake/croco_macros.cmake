###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
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
# Create subdirectories for prepared sources and extract
# destination basename
function(croco_calc_prepared_sources_path oldfile newfile_var_name)
	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# build abs path (extract dirname & fname)
	get_filename_component(oldfile_fname ${oldfile} NAME)
	get_filename_component(oldfile_dir ${oldfile} DIRECTORY)
	get_filename_component(oldfile_dirname ${oldfile_dir} NAME)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${oldfile_dirname})

	# export final path
	set(newfile_name ${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${oldfile_dirname}/${oldfile_fname})

	# export
	set(${newfile_var_name} ${newfile_name} PARENT_SCOPE)
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
		# get target path to build name from
		croco_calc_prepared_sources_path(${oldfile} oldfile_abs)

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
	message(STATUS "|  Optimization     : ${WITH_OPTIM}")
	message(STATUS "|  PSyClone         : ${PSYCLONE_COMMAND}")
	message(STATUS "|  OpenACC          : ${OPENACC}")
	message(STATUS "|  MPI              : ${MPI} - ${MPI_INCLUDE_DIR}")
	message(STATUS "|  OpenMP           : ${OPENMP}")
	message(STATUS "|  AGRIF            : ${AGRIF}")
	message(STATUS "|  Twin-checker     : ${ENABLE_TWIN_CHECKER} - ${TWIN_CHECKER_LIBRARY}")
	message(STATUS "--------------------------------------------------------------")
	message(STATUS "|  Parallelism      : ${PARALLELISM_SUMMARY}")
	message(STATUS "|  Splitting        : ${SPLITTING_X}x${SPLITTING_ETA}")
	message(STATUS "--------------------------------------------------------------")
	message(STATUS "|  CMake build type : ${CMAKE_BUILD_TYPE}")
	message(STATUS "|  CMake fflags     : ${CMAKE_BUILD_TYPE_FLAGS}")
	message(STATUS "|  Croco fflags     : ${CROCO_FORTRAN_FLAGS_STR}")
	message(STATUS "|  User fflags      : ${CMAKE_Fortran_FLAGS}")
	message(STATUS "|  All fflags       : ${CMAKE_BUILD_TYPE_FLAGS} ${CROCO_FORTRAN_FLAGS_STR} ${CMAKE_Fortran_FLAGS}")
	message(STATUS "--------------------------------------------------------------")
	message(STATUS "|  Case             : ${WITH_CASE}")
	message(STATUS "|  Keys             : ${WITH_KEYS}")
	message(STATUS "==============================================================")
endfunction()

###########################################################
# Perform some extra checks on variables to see if everything is correct
function(croco_last_checkings)
	# Require psyclone
	if (WITH_OPTIM STREQUAL "openacc-psyclone" AND NOT PSYCLONE_FOUND)
		message(FATAL_ERROR "Fail to find PSyClone, required if enabling PSyClone OpenACC via -DWITH_PSYCLONE_VENV !")
	endif ()
endfunction()

###########################################################
# create an empty cppdef.h & cppdef_dev.h in build directory
# to override (as with jobcomp) the content of the default one
#
# Todo
# ----
# If we fully move to cmake we can keep the current user way by using cppdef.h
# instead of cppdef_override.h but it requires to rename the orginal file in
# OCEAN source dir (which break jobcomp). Can also patch jobcomp to get both.
function(croco_trick_create_cpp_def_override)
	# set file names
	# ideally should be renamed (but break compat with jobcomp if we do that)
	#  - cppdefs_override.h
	#  - cppdefs_dev_override.h
	#  - param_override.h
	set(CPPDEF_OVERRIDE ${CMAKE_BINARY_DIR}/cppdefs_override.h)
	set(CPPDEF_EDIT ${CMAKE_BINARY_DIR}/cppdefs_edit.h)
	set(CPPDEF_DEV_OVERRIDE ${CMAKE_BINARY_DIR}/cppdefs_dev_override.h)
	set(PARAM_H_OVERRIDE ${CMAKE_BINARY_DIR}/param_override.h)
	set(CPPDEF_EDITOR ${CMAKE_SOURCE_DIR}/cmake/cppdef_handler_simple.py)

	# give access to them via -I
	include_directories(${CMAKE_BINARY_DIR})

	# create if needs
	execute_process(
		COMMAND ${CPPDEF_EDITOR} --format extract --case ${WITH_CASE} --keys "${WITH_KEYS}" ${CMAKE_SOURCE_DIR}/OCEAN/cppdefs.h
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		OUTPUT_FILE ${CPPDEF_EDIT}
	)

	# create if needs
	if (NOT EXISTS ${CPPDEF_OVERRIDE})
		write_file(${CPPDEF_OVERRIDE} "/* please just override the keys you needs here */")
		#file(COPY_FILE ${CMAKE_SOURCE_DIR}/OCEAN/cppdefs.h ${CPPDEF_OVERRIDE})
	endif ()

	# create if needs
	if (NOT EXISTS ${CPPDEF_DEV_OVERRIDE})
		file(COPY_FILE ${CMAKE_SOURCE_DIR}/OCEAN/cppdefs_dev.h ${CPPDEF_DEV_OVERRIDE})
	endif ()

	# create if needs
	if (NOT EXISTS ${PARAM_H_OVERRIDE})
		file(COPY_FILE ${CMAKE_SOURCE_DIR}/OCEAN/param.h ${PARAM_H_OVERRIDE})
	endif ()
endfunction(croco_trick_create_cpp_def_override)

###########################################################
# Copy the required files in the build dir so we are ready to run croco
#
# Todo
# ----
# There is currently more stuff done by the script /create_config.bash for jobcomp
# which should also bringed here. Or we should re-arranged this bash script to
# be able to all its sub-part from here.
function(croco_copy_case_files)
	if (NOT EXISTS ${CMAKE_BINARY_DIR}/TEST_CASES)
		file(COPY ${CMAKE_SOURCE_DIR}/TEST_CASES DESTINATION ${CMAKE_BINARY_DIR})
	endif()
endfunction(croco_copy_case_files)

###########################################################
# copy sources files
function(croco_copy_sources patterns_list)
  make_directory(${CROCO_SOURCE_DIR})
  foreach(pattern IN LISTS ${patterns_list})
    file(GLOB files ${CMAKE_CURRENT_SOURCE_DIR}/${pattern})
    if(files)
      # message(STATUS "Copying ${files} to ${CROCO_SOURCE_DIR}")
      file(COPY ${files} DESTINATION ${CROCO_SOURCE_DIR})
    else()
      message(WARNING "Warning : cannot find files matching: ${pattern}")
    endif()
  endforeach()
endfunction()


# User part. A no copy_list is needed as some .h are generated in CMAKE_BINARY_DIR
function(croco_copy_user_sources patterns_list no_copy_list)
  make_directory(${CROCO_SOURCE_DIR})

  foreach(pattern IN LISTS ${patterns_list})
    file(GLOB files ${CMAKE_BINARY_DIR}/${pattern})

    if(files)
      set(filtered_files)
      foreach(file IN LISTS files)
        get_filename_component(bfile ${file} NAME)
        list(FIND ${no_copy_list} ${bfile} index)
        if(index EQUAL -1)
          list(APPEND filtered_files ${file})
        endif()
      endforeach()

      if(filtered_files)
        message(STATUS "Copying ${filtered_files} to ${CROCO_SOURCE_DIR}")
        file(COPY ${filtered_files} DESTINATION ${CROCO_SOURCE_DIR})
      endif()
    endif()
  endforeach()
endfunction()

