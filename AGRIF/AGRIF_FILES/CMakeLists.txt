######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
# Apply croco flags
add_definitions    (-DHAVE_CMAKE_CONFIG)
add_compile_options(${CROCO_FORTRAN_FLAGS})
add_link_options   (${CROCO_FORTRAN_FLAGS})

######################################################
# Build AGRIF library
file(GLOB AGRIF_SRC_F90 *.F90)
add_library(agrif STATIC ${AGRIF_SRC_F90})
