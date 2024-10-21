###########################################################
# - Check is the compiler supports MPI
# Once done this will define
#  MPI_FOUND - System has MPI enabled

###########################################################
# define MPI fortran compiler
set(MPI_FORT OFF CACHE STRING "Define the MPI fortran compiler to be used.")

###########################################################
# gen test source file
set(tmp_check_mpi_path ${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.check.mpi.F90)
file(WRITE ${tmp_check_mpi_path} "
#include <mpif.h>
program test
end program tests
")

# call compiler with -M to dump #include paths
execute_process(COMMAND bash -c "${CMAKE_Fortran_COMPILER}\ -E\ -cpp\ ${tmp_check_mpi_path}\ |\ grep /mpif.h\ |\ cut\ -d '\"'\ -f\ 2"
                RESULT_VARIABLE tmp_check_mpi_status
                OUTPUT_VARIABLE tmp_check_mpi_output)

# extract patent dir
if (tmp_check_mpi_status EQUAL 0)
	get_filename_component(MPI_INCLUDE_DIR ${tmp_check_mpi_output} DIRECTORY)
endif()

###########################################################
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDFF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MPI DEFAULT_MSG MPI_INCLUDE_DIR)

###########################################################
mark_as_advanced(MPI_INCLUDE_DIR)
