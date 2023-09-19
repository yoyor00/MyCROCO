######################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
######################################################

######################################################
# Build AGRIF source converter
file(GLOB AGRIF_SRC_C *.c)
add_executable(conv ${AGRIF_SRC_C})
