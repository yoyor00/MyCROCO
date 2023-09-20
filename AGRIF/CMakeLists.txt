###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

# if agrif enabled
if (ENABLE_AGRIF)
	add_subdirectory(LIB)
	add_subdirectory(AGRIF_FILES)
endif ()
