######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
# Macro to replace some flage, like -O3 => -O2
# Params:
#   - variable : Name of the variable to affect
#   - flag     : Flag to search and replace.
#   - new_flag : New value (can be empty to erase)
macro(replace_compile_flag variable flag new_flag)
	set(tmp_list)
	foreach(lvalue IN LISTS variable)
		if (lvalue STREQUAL ${flag})
			list(APPEND tmp_list ${new_flag})
		else (lvalue STREQUAL ${flag})
			list(APPEND tmp_list ${value})
		endif (lvalue STREQUAL ${flag})
	endforeach(lvalue IN LISTS variable)
	#string(REPLACE "${flag}" "${new_flag}" ${variable} "${CMAKE_Fortran_FLAGS_RELEASE}")
endmacro()

######################################################
# Replace the include to use the one without OPENACC
# Inspired from https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
#
# TODO: REMOVED WHEN FINISHED
function(croco_psyclone_pre_filter_acc list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# build new name
		string(REGEX REPLACE "\\.F" ".no-acc.F" newfile ${oldfile})

		# build new path
		get_filename_component(newfile ${newfile} NAME)
		set(newfile ${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${newfile})

		# build command
		add_custom_command(
			OUTPUT ${newfile}
			COMMAND ${CMAKE_SOURCE_DIR}/psyclone/psyclone.preprocessor.skip.acc.sh ${oldfile} ${newfile}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS psyclone/psyclone.preprocessor.skip.acc.sh
					psyclone/skip.openacc.rules.lst
					
			VERBATIM
		)
		list(APPEND _newfiles ${newfile})
	endforeach()

	# move to user variable
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()

######################################################
# Apply cpp + mpc on the sources to prepare them before build
# Inspired from https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
function(croco_cpp_and_mpc_preprocess list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# build new name
		string(REGEX REPLACE "\\.F" ".mpc.F" newfile ${oldfile})

		# build new full path
		get_filename_component(newfile ${newfile} NAME)
		set(newfile ${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${newfile})

		# build command
		add_custom_command(
			OUTPUT ${newfile}
			COMMAND ${CROCO_FORTAN_CPP} ${CROCO_FORTRAN_CPP_FLAGS} ${oldfile} | ${CROCO_MPC} > ${newfile}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc ${OCEAN_CPP_H}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile})
	endforeach()
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()

######################################################
# Apply cpp + mpc on the sources to prepare them before build
# Inspired from https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
function(croco_change_loop_preprocess list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# build new name
		string(REGEX REPLACE "\\.F" ".loops.F" newfile ${oldfile})

		# build new full path
		get_filename_component(newfile ${newfile} NAME)
		set(newfile ${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${newfile})

		# build final name
		string(REGEX REPLACE "\\.F" ".mpc.F" newfile_mpc ${newfile})

		# build command
		add_custom_command(
			OUTPUT ${newfile_mpc}
			COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/change_loops.py ${oldfile} ${newfile}
			COMMAND ${CROCO_MPC} < ${newfile}  > ${newfile_mpc}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc ${OCEAN_CPP_H}
			BYPRODUCTS ${newfile}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile_mpc})
	endforeach()
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()
