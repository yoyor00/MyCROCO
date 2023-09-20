###########################################################
#  CROCO cmake build system, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
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
			set(CROCO_FORTRAN_FLAGS      -mcmodel=medium -fno-alias -i4 -r8 -fp-model precise)
			                             #"-O0 -g -i4 -r8 -traceback -check all -check bounds
			                             #-check uninit -CA -CB -CS -ftrapuv -fpe1"
			set(CROCO_OPTIMIZE_LEVEL     -O2)
		elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
			set(CROCO_FORTRAN_CPP "cpp")
			set(CROCO_FORTRAN_CPP_FLAGS  -traditional -P -DLinux)
			set(CROCO_FORTRAN_FLAGS      -mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy)
			                             #FFLAGS1="-O0 -g -fdefault-real-8 -fdefault-double-8 -std=legacy -fbacktrace \
			                             #-fbounds-check -finit-real=nan -finit-integer=8888"
			set(CROCO_OPTIMIZE_LEVEL -O2)
		elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
			set(CROCO_FORTRAN_CPP "cpp")
			set(CROCO_FORTRAN_CPP_FLAGS -traditional -P -DLinux -DXLF)
			set(CROCO_FORTRAN_FLAGS     -g -fast -r8 -i4 -mcmodel=medium)
			set(CROCO_OPTIMIZE_LEVEL    -O2)
		else()
			message(FATAL_ERROR "Unsupported compiler : ${CMAKE_Fortran_COMPILER}")
		endif()
	elseif (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "CYGWIN_NT-10.0")
		set(CROCO_FORTRAN_CPP "cpp")
		set(CROCO_FORTRAN_CPP_FLAGS -traditional -P -DLinux -DIfort)
		set(CROCO_FORTRAN_FLAGS     -fdefault-real-8 -fdefault-double-8 -march=native -mtune=native)
		set(CROCO_OPTIMIZE_LEVEL    -O4)
	elseif (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "AIX")
		set(CROCO_FORTRAN_CPP "cpp")
		set(CROCO_FORTRAN_CPP_FLAGS  -traditional -P -DLinux -DIfort)
		set(CROCO_FORTRAN_FLAGS      -q64 -qwarn64 -qfixed -qrealsize=8 -qintsize=8 -qhot
		                             -qalias=noaryovrlp -qthreaded -qarch=pwr4 -qtune=pwr4
		                             -qunroll=yes)
		set(CROCO_OPTIMIZE_LEVEL     -O3)
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
			list(APPEND CROCO_FORTRAN_FLAGS -acc=gpu -Minfo=accel)
		else()
			message(FATAL_ERROR "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
			                    "Unsupported compiler (${CMAKE_Fortran_COMPILER}) to use OpenACC !\n"
			                    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		endif()
	endif()

	#######################################################
	# Assign default not to let cmake using -O3 which is default
	list(JOIN CROCO_FORTRAN_FLAGS " " CROCO_FORTRAN_FLAGS_STR)
	if (CMAKE_Fortran_FLAGS STREQUAL "")
		set(CMAKE_Fortran_FLAGS "-DNDEBUG ${CROCO_OPTIMIZE_LEVEL}")
	endif()
	set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG ${CROCO_OPTIMIZE_LEVEL}")
	set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-DNDEBUG ${CROCO_OPTIMIZE_LEVEL}")
endmacro()
