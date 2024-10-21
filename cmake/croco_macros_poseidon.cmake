###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Use de `poseidon` command to regenerate the math sources
# and to tune them in another way.
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
function(croco_poseidon_pre_process list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/poseidon)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# extract dir
		get_filename_component(oldfile_dir ${oldfile} DIRECTORY)
		get_filename_component(oldfile_dirname ${oldfile_dir} NAME)

		# Here we remove the ts_2d/* files from the list to replace all of them
		# by the poseidon regenerated single file (step2d.poseidon.F90)
		# [done after the loop].
		#
		# LIMITATION : Currently working only on BASIN with SOLVE3D enabled (default)
		if (oldfile_dirname STREQUAL "ts_2d")
			list(APPEND step2d_files ${oldfile})
			list(REMOVE_ITEM ${list_to_update} ${oldfile})
		else()
			list(APPEND _newfiles ${oldfile})
		endif()
	endforeach()

	# get deps
	# TODO: see if we keep this dep, this is more for POSEIDON dev than really usefull
	# ony work if posidon was installed with `pip install -e .`
	file(GLOB_RECURSE POSEIDON_SCRIPTS "${CMAKE_SOURCE_DIR}/POSEIDON/*.py")

	# build to generate with posidon for step2d -> step2d_poseidon.F90 & step2d_poseidon_decl_croco_vars.h
	croco_calc_prepared_sources_path("${oldfile_dir}/prepared_sources/poseidon/step2d.poseidon" step2d_new_single_file)
	add_custom_command(
		OUTPUT ${step2d_new_single_file}.F90
		COMMAND poseidon dag regen -o ${step2d_new_single_file} ${step2d_files}
		DEPENDS ${step2d_files}
				${POSEIDON_SCRIPTS}
				${CROCO_CPP_H}
		VERBATIM
	)
	list(APPEND _newfiles ${step2d_new_single_file}.F90)

	# move to user variable
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()
