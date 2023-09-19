######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

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
			COMMAND ${CMAKE_SOURCE_DIR}/PSYCLONE/psyclone.preprocessor.skip.acc.sh ${oldfile} ${newfile}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS ${CMAKE_SOURCE_DIR}/PSYCLONE/psyclone.preprocessor.skip.acc.sh
			        ${CMAKE_SOURCE_DIR}/PSYCLONE/skip.openacc.rules.lst
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
			COMMAND ${CROCO_FORTRAN_CPP} ${CROCO_FORTRAN_CPP_FLAGS} ${oldfile} -o ${newfile}.tmp.F 
			COMMAND ${CROCO_MPC} < ${newfile}.tmp.F > ${newfile}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc ${CROCO_CPP_H}
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
			COMMAND ${CROCO_MPC} < ${newfile}  > ${newfile_mpc}.to_fix.F
			COMMAND egrep -v "^ +& *$" ${newfile_mpc}.to_fix.F > ${newfile_mpc} || cp ${newfile_mpc}.to_fix.F ${newfile_mpc} 
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc ${OCEAN_CPP_H}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile_mpc})
	endforeach()

	# export
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()

######################################################
# Loop on all values of the list and make them absolute path
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