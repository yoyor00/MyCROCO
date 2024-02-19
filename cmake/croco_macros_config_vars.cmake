###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Handle the WITH_OPTIM value the correct way by checking its value
# and applying some side effects we requires
#
# Parameters
# ----------
# optim:
#     The optimization values received via WITH_OPTIM to check and apply.
function(croco_check_and_config_optim optim)
	# avail
	set(optims_avail OFF seq openmp mpi openacc-native openacc-psyclone)

	# check in
	if (NOT optim IN_LIST optims_avail)
		message(FATAL_ERROR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		                    "Invalid optimization via -DWITH_OPTIM : ${optim} !\n"
		                    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	endif()

	# apply some side effects
	# Enable some vars for config_post.h.in
	if (optim STREQUAL "openmp")
		set(OPENMP ON PARENT_SCOPE)
	endif()
	if (optim STREQUAL "openacc-native" OR WITH_OPTIM STREQUAL "openacc-psyclone")
		set(OPENACC ON PARENT_SCOPE)
	endif()
	if (optim STREQUAL "mpi")
		set(MPI ON PARENT_SCOPE)
	endif ()

	# more for printing summary as currently there is match between optim & parallel
	# but will change in future
	set(PARALLELISM_SUMMARY ${optim} PARENT_SCOPE)
endfunction()

###########################################################
# Set some vars depending on what has been asked by users
function(croco_transpose_config_to_internal_vars)
	# Trigger the case definition (macro variable having the name of the case)
	set(${WITH_CASE} ON PARENT_SCOPE)
	set(AGRIF ${ENABLE_AGRIF} PARENT_SCOPE)

	# extract NP_x & NP_eta
	string(REPLACE "x" ";" WITH_SPLITTING_ARRAY ${WITH_SPLITTING})
	list(GET WITH_SPLITTING_ARRAY 0 SPLITTING_X)
	list(GET WITH_SPLITTING_ARRAY 1 SPLITTING_ETA)

	# export
	set(SPLITTING_X ${SPLITTING_X} PARENT_SCOPE)
	set(SPLITTING_ETA ${SPLITTING_ETA} PARENT_SCOPE)
endfunction()
