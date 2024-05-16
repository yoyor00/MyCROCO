###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Apply the AGRIF transformation steps to prepare the sources before building them.
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
function(croco_agrif_preprocess list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/AGRIF_MODELFILES)
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/AGRIF_INCLUDES)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# get target path to build name from
		croco_calc_prepared_sources_path(${oldfile} oldfile_abs)

		# build new name
		string(REGEX REPLACE "\\.F" ".cpp.F"   newfile_cpp           ${oldfile_abs})
		string(REGEX REPLACE "\\.F" ".agrif.F" newfile_cpp_agrif     ${newfile_cpp})
		string(REGEX REPLACE "\\.F" ".mpc.F"   newfile_cpp_agrif_mpc ${newfile_cpp_agrif})

		# build command
		add_custom_command(
			OUTPUT ${newfile_cpp_agrif_mpc}
			COMMAND ${CROCO_FORTRAN_CPP} ${CROCO_FORTRAN_CPP_FLAGS} ${oldfile} -o ${newfile_cpp}
			COMMAND ${CROCO_AGRIF_CONV} amr.scrum -rm -comdirin ${CMAKE_CURRENT_BINARY_DIR} -incdir ${CMAKE_CURRENT_BINARY_DIR}/AGRIF_INCLUDES/ -comdirout ${CMAKE_CURRENT_BINARY_DIR}/AGRIF_MODELFILES/ -convfile ${newfile_cpp} > ${newfile_cpp_agrif}
			COMMAND ${CROCO_MPC} < ${newfile_cpp_agrif}  > ${newfile_cpp_agrif_mpc}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS conv ${CROCO_CPP_H}
			VERBATIM
		)
		list(APPEND _newfiles ${newfile_cpp_agrif_mpc})
	endforeach()

	# export
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()
