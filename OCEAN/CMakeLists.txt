######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
# Add the source files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})

######################################################
# -I directo
#include_directories(${CMAKE_CURRENT_BUILD_DIR})
include_directories(.)
include_directories(/home/svalat/usr/include/)
set(CPP_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR} -I/home/svalat/usr/include/)

######################################################
# configure build options
add_definitions(-DLinux)
add_compile_options(-mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy -P)

######################################################
# Preprocessor options
set(CROCO_FORTAN_CPP_FLAGS -traditional -DLinux -P ${COMPILE_DEFINITIONS} ${CPP_INCLUDES})
set(CROCO_FORTAN_CPP cpp)

######################################################
# path to "mpc", the ocean internal source preprocessor
set(CROCO_MPC ${CMAKE_CURRENT_BINARY_DIR}/mpc)

######################################################
# Get list of files
file(GLOB OCEAN_CPP_H *.h)
file(GLOB OCEAN_SRC_F *.F)
file(GLOB OCEAN_SRC_F90 *.F90)

# we remove the unsed files or used for other binaries
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
)
list(REMOVE_ITEM OCEAN_SRC_F90 ${CMAKE_CURRENT_SOURCE_DIR}/module_checkmpi.F90 ${CMAKE_CURRENT_SOURCE_DIR}/module_qsort.F90 ${CMAKE_CURRENT_SOURCE_DIR}/module_tracetxt_out.F90)

######################################################
# Apply cpp + mpc on the sources to prepare them before build
# Inspired from https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
function(croco_cpp_and_mpc_preprocess oldfiles newfiles)
	set(_newfiles)
	foreach(oldfile IN LISTS oldfiles)
		message(STATUS ${oldfile})
		string(REGEX REPLACE "\\.F" ".preprocessed.F" newfile ${oldfile})
		get_filename_component(newfile ${newfile} NAME)
		set(newfile ${CMAKE_CURRENT_BINARY_DIR}/${newfile})
		add_custom_command(
			OUTPUT ${newfile}
			COMMAND ${CROCO_FORTAN_CPP} ${CROCO_FORTAN_CPP_FLAGS} ${oldfile} | ${CROCO_MPC} > ${newfile}
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
# build croco ocean simulation
croco_cpp_and_mpc_preprocess("${OCEAN_SRC_F}" OCEAN_SRC_F_PRE_PROCESSED)
set(CMAKE_Fortran_PREPROCESS ON)
add_executable(ocean ${OCEAN_SRC_F90} ${OCEAN_SRC_F_PRE_PROCESSED})
target_link_libraries(ocean /home/svalat/usr/lib/libnetcdff.so)

######################################################
# build ncjoin utility
set(OCEAN_NCJOIN_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncjoin.F)
croco_cpp_and_mpc_preprocess("${OCEAN_NCJOIN_SRC_F}" OCEAN_NCJOIN_SRC_F_PRE_PROCESSED)
add_executable(ncjoin ${OCEAN_NCJOIN_SRC_F_PRE_PROCESSED} lenstr.preprocessed.F)
target_link_libraries(ncjoin /home/svalat/usr/lib/libnetcdff.so)

######################################################
# build ncrename utility
set(OCEAN_NCRENAME_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncrename.F)
croco_cpp_and_mpc_preprocess("${OCEAN_NCRENAME_SRC_F}" OCEAN_NCRENAME_SRC_F_PRE_PROCESSED)
add_executable(ncrename ${OCEAN_NCRENAME_SRC_F_PRE_PROCESSED} lenstr.preprocessed.F)
target_link_libraries(ncrename /home/svalat/usr/lib/libnetcdff.so)

######################################################
# build mpc utility
add_executable(mpc mpc.F)

######################################################
# build checkkwds utility (used to produce setup_kwds.F)
add_executable(checkkwds checkkwds.F)

######################################################
# build checkkwds utility (used to produce setup_kwds.F)
add_executable(cppcheck cppcheck.F)

######################################################
# build checkkwds utility (used to produce check_srcs.F)
add_executable(srcscheck srcscheck.F)

######################################################
# build cross_matrix utility
add_executable(cross_matrix cross_matrix.F)

######################################################
# build partit utility
add_executable(partit partit.F insert_node.preprocessed.F lenstr.preprocessed.F)
target_link_libraries(partit /home/svalat/usr/lib/libnetcdff.so)
