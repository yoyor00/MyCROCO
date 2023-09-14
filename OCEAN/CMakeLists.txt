######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
# Add the source files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})

######################################################
# generate config.h
configure_file(config.h.in config.h)

######################################################
# -I directives
#include_directories(${CMAKE_CURRENT_BUILD_DIR})
include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${NETCDFF_INCLUDE_DIR})
set(CPP_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_BINARY_DIR} -I/home/svalat/usr/include/)

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
file(GLOB OCEAN_CPP_H *.h *.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
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
croco_cpp_and_mpc_preprocess(OCEAN_NCJOIN_SRC_F)
add_executable(ncjoin ${OCEAN_NCJOIN_SRC_F})
target_link_libraries(ncjoin ${NETCDFF_LIBRARY})
target_link_libraries(ncjoin croco_common)
set_property(TARGET ncjoin PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build ncrename utility
set(OCEAN_NCRENAME_SRC_F ${CMAKE_CURRENT_SOURCE_DIR}/ncrename.F)
croco_cpp_and_mpc_preprocess(OCEAN_NCRENAME_SRC_F)
add_executable(ncrename ${OCEAN_NCRENAME_SRC_F})
target_link_libraries(ncrename ${NETCDFF_LIBRARY})
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
target_link_libraries(partit ${NETCDFF_LIBRARY})
target_link_libraries(partit croco_common)
set_property(TARGET partit PROPERTY LINKER_LANGUAGE Fortran)

######################################################
# build croco ocean simulation

# TODO: remove when finished
if (ENABLE_PSYCLONE_OPENACC)
	croco_psyclone_pre_filter_acc(OCEAN_SRC_F)
endif (ENABLE_PSYCLONE_OPENACC)

# build
croco_cpp_and_mpc_preprocess(OCEAN_SRC_F)
croco_change_loop_preprocess(OCEAN_SRC_F)
set(CMAKE_Fortran_PREPROCESS ON)
add_executable(croco ${CROCO_SRC_F90} ${OCEAN_SRC_F})
target_link_libraries(croco ${NETCDFF_LIBRARY})
target_link_libraries(croco croco_common)
set_property(TARGET croco PROPERTY LINKER_LANGUAGE Fortran)

# Inject psyclone in the loop if enabled
if (ENABLE_PSYCLONE_OPENACC)
	set_target_properties(croco PROPERTIES Fortran_COMPILER_LAUNCHER ${CMAKE_SOURCE_DIR}/psyclone/psyclone.compiler.wrapper.py)
endif (ENABLE_PSYCLONE_OPENACC)

######################################################
# Make some of the executables landing directly in root dir
set_property(TARGET partit croco cross_matrix PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
