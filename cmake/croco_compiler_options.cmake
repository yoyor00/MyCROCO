###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
###########################################################

###########################################################
# Determine what are the default flags to use depending on the
# Compiler & OS
#
# Output variables
# ----------------
# CROCO_FORTRAN_CPP:
#     The `cpp` command to be used to pre-process the code.
# CROCO_FORTRAN_CPP_FLAGS:
#     The flags to be given to `cpp` in addition to -I... -D...
# CROCO_FORTRAN_FLAGS:
#     The tuned compile flags to be used to build fortran files.
macro(croco_tune_compile_flags)
	######################################################
	if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
		if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
			set(CROCO_FORTRAN_CPP "cpp")
			set(CROCO_FORTRAN_CPP_FLAGS  -traditional -P -DLinux -DIfort)
			set(CROCO_FORTRAN_FLAGS      )
			set(CROCO_OPTIMIZE_LEVEL     )
		elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
			set(CROCO_FORTRAN_CPP "cpp")
			set(CROCO_FORTRAN_CPP_FLAGS  -traditional -P -DLinux)
			set(CROCO_FORTRAN_FLAGS      )
			set(CROCO_OPTIMIZE_LEVEL )
		elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
			set(CROCO_FORTRAN_CPP "cpp")
			set(CROCO_FORTRAN_CPP_FLAGS -traditional -P -DLinux -DXLF)
			set(CROCO_FORTRAN_FLAGS     )
			set(CROCO_OPTIMIZE_LEVEL    )
		else()
			message(FATAL_ERROR "Unsupported compiler : ${CMAKE_Fortran_COMPILER}")
		endif()
	else()
		message(FATAL_ERROR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"Unsupported operating system (uname -s) : ${CMAKE_HOST_SYSTEM_NAME}\n"
			                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
	endif()

	#######################################################
	# Specific sites
	site_name(SITE_NAME)
	if (SITE_NAME MATCHES "jean-zay.*")
		list(APPEND CROCO_FORTRAN_CPP_FLAGS "-DJEANZAY")
	endif ()

	#######################################################
	# OpenMP is enabled
	if (OPENMP)
		list(APPEND CROCO_FORTRAN_FLAGS ${OpenMP_Fortran_FLAGS})
	endif()

	#######################################################
	# OpenACC is enabled
	if (OPENACC)
		if(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
			list(APPEND CROCO_FORTRAN_FLAGS -acc=gpu)
			#list(APPEND CROCO_FORTRAN_FLAGS -Minfo=accel)
		else()
			message(FATAL_ERROR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
			                    "Unsupported compiler (${CMAKE_Fortran_COMPILER}) to use OpenACC, please use FC=nvfortran !\n"
			                    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		endif()
	endif()

	#######################################################
	# Assign default not to let cmake using -O3 which is default
	list(JOIN CROCO_FORTRAN_FLAGS " " CROCO_FORTRAN_FLAGS_STR)
	list(JOIN CROCO_OPTIMIZE_LEVEL " " CROCO_OPTIMIZE_LEVEL_STR)
	if (CMAKE_Fortran_FLAGS STREQUAL "")
		set(CMAKE_Fortran_FLAGS "-DNDEBUG ${CROCO_OPTIMIZE_LEVEL_STR}")
	endif()
	set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG ${CROCO_OPTIMIZE_LEVEL_STR}")
	set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-DNDEBUG -g ${CROCO_OPTIMIZE_LEVEL_STR}")
endmacro()
