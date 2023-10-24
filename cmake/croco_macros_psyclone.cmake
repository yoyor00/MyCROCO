###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Replace the includes to make the files listing in 
# PSYCLONE/skip.openacc.rules.lst not applying the OPENACC
# macro to keep them sequential.
#
# It is a temporary solution to keep working with psyclone
# from the relady GPU code.
#
# TODO
# ----
# Remove this when going back to the master branch
# if not keeping the manual GPU code.
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
function(croco_psyclone_pre_filter_acc list_to_update)
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
		string(REGEX REPLACE "\\.F" ".no-acc.F" newfile_no_acc ${oldfile_abs})

		# get deps
		file(GLOB_RECURSE PSYCLONE_SCRIPTS "${CMAKE_SOURCE_DIR}/PSYCLONE/*.py")

		# build command
		add_custom_command(
			OUTPUT ${newfile_no_acc}
			COMMAND ${CMAKE_SOURCE_DIR}/PSYCLONE/psyclone.preprocessor.skip.acc.py ${oldfile} ${newfile_no_acc}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS ${CMAKE_SOURCE_DIR}/PSYCLONE/psyclone.preprocessor.skip.acc.py
			        ${CMAKE_SOURCE_DIR}/PSYCLONE/psyclone.rules.json
					${PSYCLONE_SCRIPTS}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile_no_acc})
	endforeach()

	# move to user variable
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()

###########################################################
# Apply the change_loop.py script on the manual GPU
# code to reshape some loops.
#
# TODO
# ----
# Remove this when going back to the master branch
# if not keeping the manual GPU code.
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
function(croco_change_loop_preprocess list_to_update)
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
		string(REGEX REPLACE "\\.F" ".loops.F" newfile_loops         ${oldfile_abs})
		string(REGEX REPLACE "\\.F" ".mpc.F"   newfile_loops_mpc     ${newfile_loops})
		string(REGEX REPLACE "\\.F" ".fix.F"   newfile_loops_mpc_fix ${newfile_loops_mpc})

		# get deps
		file(GLOB_RECURSE PSYCLONE_SCRIPTS "${CMAKE_SOURCE_DIR}/PSYCLONE/*.py")

		# build command
		add_custom_command(
			OUTPUT ${newfile_loops_mpc_fix}
			COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/change_loops.py ${oldfile} ${newfile_loops}
			COMMAND ${CROCO_MPC} < ${newfile_loops}  > ${newfile_loops_mpc}
			COMMAND egrep -v "^ +& *$" ${newfile_loops_mpc} > ${newfile_loops_mpc_fix} || cp ${newfile_loops_mpc} ${newfile_loops_mpc_fix} 
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc
			        ${OCEAN_CPP_H}
			        ${CMAKE_CURRENT_SOURCE_DIR}/change_loops.py
			        ${PSYCLONE_SCRIPTS}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile_loops_mpc_fix})
	endforeach()

	# export
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()
