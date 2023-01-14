# Add the source files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BUILD_DIR})
include_directories(${CMAKE_CURRENT_BUILD_DIR})
#include_directories(${CMAKE_CURRENT_BUILD_DIR})
file(GLOB OCEAN_SRC *.F90 *.F)
list(REMOVE_ITEM OCEAN_SRC mpc.F)
file(GLOB OCEAN_CPP_H *.h)

add_executable(mpc mpc.F)

add_definitions(-DLinux)
add_compile_options(-mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy -P)
set(CROCO_FORTAN_CPP_FLAGS -traditional -DLinux -P ${COMPILE_DEFINITIONS})

#https://fortran.cat/2021/09/24/cmake-and-fypp-preprocessor/
function(croco_cpp_and_mpc_preprocess oldfiles newfiles)
	set(_newfiles)
	foreach(oldfile IN LISTS oldfiles)
		string(REGEX REPLACE "\\.F" ".preprocessed.F" newfile ${oldfile})
		add_custom_command(
			OUTPUT ${newfile}
			COMMAND ${CROCO_FORTAN_CPP} ${CROCO_FORTAN_CPP_FLAGS} ${oldfile} | ${CROCO_MPC} > ${newfile}
			MAIN_DEPENDENCY ${oldfile}
			DEPENDS mpc ${OCEAN_CPP_H}
		)
		list(APPEND _newfiles ${newfile})
	endforeach()
	set(${newfiles} ${_newfiles} PARENT_SCOPE)
endfunction()

croco_cpp_and_mpc_preprocess("${OCEAN_SRC}" OCEAN_SRC_PRE_PROCESSED)
set(CMAKE_Fortran_PREPROCESS OFF)
add_executable(ocean ${OCEAN_SRC_PRE_PROCESSED})
