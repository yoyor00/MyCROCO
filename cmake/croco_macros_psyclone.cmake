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
		# get target path to build name from
		croco_calc_prepared_sources_path(${oldfile} oldfile_abs)

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
		# get target path to build name from
		croco_calc_prepared_sources_path(${oldfile} oldfile_abs)

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


function(croco_psyclone_twin_checker list_to_update)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS ${list_to_update})
		# get target path to build name from
		croco_calc_prepared_sources_path(${oldfile} oldfile_abs)

		# build new name
		string(REGEX REPLACE "\\.F" ".twin_check.F90" newfile_twin_check ${oldfile_abs})

		# set
		set(twin_script ${CMAKE_SOURCE_DIR}/TWINCHECKER/twin_psyclone_script.py)

		# set filteres
		get_filename_component(oldfile_fname ${oldfile} NAME)
		set(allowed_files step2d.cpp.mpc.F zetabc.cpp.mpc.F u2dbc.cpp.mpc.F v2dbc.cpp.mpc.F step3d_t.cpp.mpc.F t3dbc.cpp.mpc.F step3d_uv2.cpp.mpc.F step3d_uv1.cpp.mpc.F rhs3d.cpp.mpc.F pre_step3d.cpp.mpc.F u3dbc.cpp.mpc.F v3dbc.cpp.mpc.F t3dmix.cpp.mpc.F wvlcty.cpp.mpc.F setup_grid1.cpp.mpc.F setup_grid2.cpp.mpc.F rho_eos.cpp.mpc.F prsgrd.cpp.mpc.F omega.cpp.mpc.F grid_stiffness.cpp.mpc.F get_stflux.cpp.mpc.F get_vbc.cpp.mpc.F diag.cpp.mpc.F analytical.cpp.mpc.F ana_initial.cpp.mpc.F get_srflux.cpp.mpc.F ana_grid.cpp.mpc.F)
		set(allowed_files ${allowed_files} set_avg.cpp.mpc.F set_bio_diags_avg.cpp.mpc.F set_cycle.cpp.mpc.F set_depth.cpp.mpc.F set_diags_avg.cpp.mpc.F set_diags_eddy_avg.cpp.mpc.F set_diags_ek_avg.cpp.mpc.F set_diags_ek.cpp.mpc.F set_diagsM_avg.cpp.mpc.F set_diags_pv_avg.cpp.mpc.F set_diags_pv.cpp.mpc.F set_diags_vrt_avg.cpp.mpc.F set_diags_vrt.cpp.mpc.F set_nudgcof.cpp.mpc.F set_nudgcof_fine.cpp.mpc.F set_scoord.cpp.mpc.F set_surf_avg.cpp.mpc.F set_weights.cpp.mpc.F)
		set(allowed_files ${allowed_files} init_arrays.cpp.mpc.F)
		#set(allowed_files omega.cpp.mpc.F)
		#set(allowed_files analytical.cpp.mpc.F)

		# build command
		if (${oldfile_fname} IN_LIST allowed_files)
			add_custom_command(
				OUTPUT ${newfile_twin_check}
				COMMAND ${twin_script} ${oldfile} --output ${newfile_twin_check}
				MAIN_DEPENDENCY ${oldfile}
				DEPENDS ${twin_script}
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