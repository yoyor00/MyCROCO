######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
# Add the source files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})

######################################################
# -I directives
#include_directories(${CMAKE_CURRENT_BUILD_DIR})
include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(/home/svalat/usr/include/)
set(CPP_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_BINARY_DIR} -I/home/svalat/usr/include/)

######################################################
macro(replace_compile_flag variable flag new_flag)
	string(REPLACE "${flag}" "${new_flag}" ${variable} "${CMAKE_Fortran_FLAGS_RELEASE}")
endmacro()

######################################################
# configure build options
add_definitions(-DLinux)
if (CMAKE_Fortran_COMPILER MATCHES "nvfortran$")
	add_compile_options(-g -fast -r8 -i4 -mcmodel=medium -acc=gpu)
	add_link_options(-g -fast -r8 -i4 -mcmodel=medium -acc=gpu)
	# Note:
	# there is a buf with nvfortran which makes it crashing with O3 on somes files
	replace_compile_flag(CMAKE_Fortran_FLAGS_RELEASE "-O3" "-O2")
	replace_compile_flag(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3" "-O2")
else(CMAKE_Fortran_COMPILER MATCHES "nvfortran$")
	add_compile_options(-mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy -P)
endif(CMAKE_Fortran_COMPILER MATCHES "nvfortran$")

######################################################
# Preprocessor options
set(CROCO_FORTRAN_CPP_FLAGS -traditional -DLinux -P ${COMPILE_DEFINITIONS} ${CPP_INCLUDES})
set(CROCO_FORTAN_CPP cpp)

######################################################
# path to "mpc", the ocean internal source preprocessor
set(CROCO_MPC ${CMAKE_CURRENT_BINARY_DIR}/mpc)

######################################################
# Get list of files
file(GLOB OCEAN_CPP_H *.h)
file(GLOB OCEAN_SRC_F *.F)
file(GLOB OCEAN_SRC_F90 *.F90)

# we remove the unused files or used for other binaries
list(REMOVE_ITEM OCEAN_SRC_F 
	${CMAKE_CURRENT_SOURCE_DIR}/mpc.F
	${CMAKE_CURRENT_SOURCE_DIR}/lmd_skpp1994.F
	${CMAKE_CURRENT_SOURCE_DIR}/uv3dmix4_GP.F
	${CMAKE_CURRENT_SOURCE_DIR}/uv3dmix4_S.F
	${CMAKE_CURRENT_SOURCE_DIR}/bio_BioEBUS.F
	${CMAKE_CURRENT_SOURCE_DIR}/bio_NChlPZD.F
	${CMAKE_CURRENT_SOURCE_DIR}/bio_N2ChlPZD2.F
	${CMAKE_CURRENT_SOURCE_DIR}/cppcheck.F
	${CMAKE_CURRENT_SOURCE_DIR}/cross_matrix.F
	${CMAKE_CURRENT_SOURCE_DIR}/srcscheck.F
	${CMAKE_CURRENT_SOURCE_DIR}/checkkwds.F
	${CMAKE_CURRENT_SOURCE_DIR}/partit.F
	${CMAKE_CURRENT_SOURCE_DIR}/ncjoin.F
	${CMAKE_CURRENT_SOURCE_DIR}/ncrename.F
	${CMAKE_CURRENT_SOURCE_DIR}/lmd_bkpp1994.F
	${CMAKE_CURRENT_SOURCE_DIR}/lmd_bkpp2005.F
	${CMAKE_CURRENT_SOURCE_DIR}/lmd_skpp2005.F
	${CMAKE_CURRENT_SOURCE_DIR}/t3dmix_ISO.F
	${CMAKE_CURRENT_SOURCE_DIR}/t3dmix_S.F
	${CMAKE_CURRENT_SOURCE_DIR}/uv3dmix_S.F
	${CMAKE_CURRENT_SOURCE_DIR}/srcscheck.F
	${CMAKE_CURRENT_SOURCE_DIR}/setup_kwds.F
	${CMAKE_CURRENT_SOURCE_DIR}/check_srcs.F
	${CMAKE_CURRENT_SOURCE_DIR}/check_switches1.F
	${CMAKE_CURRENT_SOURCE_DIR}/testkeys.F
	${CMAKE_CURRENT_SOURCE_DIR}/uv3dmix_GP.F
	${CMAKE_CURRENT_SOURCE_DIR}/lenstr.F
	${CMAKE_CURRENT_SOURCE_DIR}/insert_node.F
)
list(REMOVE_ITEM OCEAN_SRC_F90 
	${CMAKE_CURRENT_SOURCE_DIR}/module_checkmpi.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/module_qsort.F90
	${CMAKE_CURRENT_SOURCE_DIR}/module_tracetxt_out.F90
)

######################################################
# Replace the include to use the one without OPENACC
# Inspired from https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
#
# TODO: REMOVED WHEN FINISHED
function(croco_psyclone_pre_filter_acc oldfiles newfiles)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS oldfiles)
		# build new name
		string(REGEX REPLACE "\\.F" ".no-acc.F" newfile ${oldfile})

		# build new path
		get_filename_component(newfile ${newfile} NAME)
		set(newfile ${CMAKE_CURRENT_BINARY_DIR}/prepared_sources/${newfile})

		# build command
		add_custom_command(
			OUTPUT ${newfile}
			COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/psyclone/psyclone.preprocessor.skip.acc.sh ${oldfile} ${newfile}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS psyclone/psyclone.preprocessor.skip.acc.sh
			        psyclone/skip.openacc.rules.lst
					
			VERBATIM
		)
		list(APPEND _newfiles ${newfile})
	endforeach()

	# move to user variable
	set(${newfiles} ${_newfiles} PARENT_SCOPE)
endfunction()

######################################################
# Apply cpp + mpc on the sources to prepare them before build
# Inspired from https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
function(croco_cpp_and_mpc_preprocess oldfiles newfiles)
	# reset the list
	set(_newfiles)

	# create dir not to trash everything in the root dir
	make_directory(${CMAKE_CURRENT_BINARY_DIR}/prepared_sources)

	# loop on all files
	foreach(oldfile IN LISTS oldfiles)
		# build new name
		string(REGEX REPLACE "\\.F" ".preprocessed.F" newfile ${oldfile})

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
	set(${newfiles} ${_newfiles} PARENT_SCOPE)
endfunction()

######################################################
# generate some sources files byb using some binary programs

# setup_kwds.F
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/setup_kwds.F
	COMMAND checkkwds && mv setup_kwds.F ${CMAKE_CURRENT_BINARY_DIR}/setup_kwds.F
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS checkkwds
	VERBATIM
)

# check_srcs.F
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/check_srcs.F
	COMMAND srcscheck
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	DEPENDS srcscheck
	VERBATIM
)

# check_switches1.F
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/check_switches1.F
	COMMAND cat cppdefs.h cppdefs_dev.h > mergcpp.txt && ${CMAKE_CURRENT_BINARY_DIR}/cppcheck && mv check_switches1.F ${CMAKE_CURRENT_BINARY_DIR}/check_switches1.F
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS cppcheck cppdefs.h cppdefs_dev.h
	VERBATIM
)

# Append the generated file to ocean croco
list(APPEND OCEAN_SRC_F ${CMAKE_CURRENT_BINARY_DIR}/setup_kwds.F)
list(APPEND OCEAN_SRC_F ${CMAKE_CURRENT_BINARY_DIR}/check_srcs.F)
list(APPEND OCEAN_SRC_F ${CMAKE_CURRENT_BINARY_DIR}/check_switches1.F)

######################################################
# ISSUE: remark, this is to solve an issue due to the trick
# of pre-applying MPC for "ocean" which is we let lenstr.F
# in the list, make add_executable(partit) skipping it.
add_library(croco_common STATIC lenstr.F insert_node.F)

######################################################
# build ncjoin utility
set(OCEAN_NCJOIN_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncjoin.F)
croco_cpp_and_mpc_preprocess("${OCEAN_NCJOIN_SRC_F}" OCEAN_NCJOIN_SRC_F_PRE_PROCESSED)
add_executable(ncjoin ${OCEAN_NCJOIN_SRC_F_PRE_PROCESSED})
target_link_libraries(ncjoin /home/svalat/usr/lib/libnetcdff.so)
target_link_libraries(ncjoin croco_common)
set_property(TARGET ncjoin PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build ncrename utility
set(OCEAN_NCRENAME_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncrename.F)
croco_cpp_and_mpc_preprocess("${OCEAN_NCRENAME_SRC_F}" OCEAN_NCRENAME_SRC_F_PRE_PROCESSED)
add_executable(ncrename ${OCEAN_NCRENAME_SRC_F_PRE_PROCESSED})
target_link_libraries(ncrename /home/svalat/usr/lib/libnetcdff.so)
target_link_libraries(ncrename croco_common)
set_property(TARGET ncrename PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build mpc utility
add_executable(mpc mpc.F)
set_property(TARGET mpc PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build checkkwds utility (used to produce setup_kwds.F)
add_executable(checkkwds checkkwds.F)
set_property(TARGET checkkwds PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build checkkwds utility (used to produce setup_kwds.F)
add_executable(cppcheck cppcheck.F)
set_property(TARGET cppcheck PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build checkkwds utility (used to produce check_srcs.F)
add_executable(srcscheck srcscheck.F)
set_property(TARGET srcscheck PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build cross_matrix utility
add_executable(cross_matrix cross_matrix.F)
set_property(TARGET cross_matrix PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build partit utility
add_executable(partit partit.F)
target_link_libraries(partit /home/svalat/usr/lib/libnetcdff.so)
target_link_libraries(partit croco_common)
set_property(TARGET partit PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build croco ocean simulation
croco_psyclone_pre_filter_acc("${OCEAN_SRC_F}" CROCO_SRC_F_PRE_PROCESSED_NO_ACC)
croco_cpp_and_mpc_preprocess("${CROCO_SRC_F_PRE_PROCESSED_NO_ACC}" CROCO_SRC_F_PRE_PROCESSED_MPC)
set(CMAKE_Fortran_PREPROCESS ON)
add_executable(croco ${CROCO_SRC_F90} ${CROCO_SRC_F_PRE_PROCESSED_MPC})
target_link_libraries(croco /home/svalat/usr/lib/libnetcdff.so)
target_link_libraries(croco croco_common)
target_compile_options(croco PRIVATE -acc=gpu)
set_target_properties(croco PROPERTIES Fortran_COMPILER_LAUNCHER ${CMAKE_CURRENT_SOURCE_DIR}/psyclone/psyclone.compiler.wrapper.py)
set_property(TARGET croco PROPERTY LINKER_LANGUAGE Fortran)
