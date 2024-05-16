###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Apply the twin checker routine injection via psyclone.
#
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
function(croco_psyclone_twin_checker list_to_update)
	# reset the list
	set(_newfiles)

	# set
	set(twin_config ${CMAKE_SOURCE_DIR}/PSYCLONE/twin-checker-config.jsonc)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# get list of files to treat
	exec_program(twin-checker-file-list ARGS --config ${twin_config} OUTPUT_VARIABLE allow_files_str)
	string(REPLACE " " ";" allowed_files ${allow_files_str})

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# get target path to build name from
		croco_calc_prepared_sources_path(${oldfile} oldfile_abs)

		# build new name
		string(REGEX REPLACE "\\.F" ".twin_check.F90" newfile_twin_check ${oldfile_abs})

		# simplify names
		# extract base name
		get_filename_component(oldfile_fname ${oldfile} NAME)
		string(REPLACE "." ";" simplified_name_parts ${oldfile_fname})
		list(GET simplified_name_parts 0 simplified_name)
		set(simplified_name ${simplified_name}.F)

		# build command
		if (${simplified_name} IN_LIST allowed_files)
			add_custom_command(
				OUTPUT ${newfile_twin_check}
				COMMAND twin-checker-psyclone-instru --config ${twin_config} --output ${newfile_twin_check} ${oldfile}
				MAIN_DEPENDENCY ${oldfile}
				DEPENDS ${twin_config}
				VERBATIM
			)
		else()
			set(newfile_twin_check ${oldfile})
		endif()
		list(APPEND _newfiles ${newfile_twin_check})
	endforeach()

	# move to user variable
	set(${list_to_update} ${_newfiles} PARENT_SCOPE)
endfunction()