######################################################
# - Check is the compiler supports MPI
# Once done this will define
#  MPI_FOUND - System has MPI enabled

######################################################
# Needs to import it
include(CheckFortranSourceCompiles)

######################################################
# define MPI fortran compiler
set(MPI_FORT OFF CACHE STRING "Define the MPI fortran compiler to be used.")

######################################################
check_fortran_source_compiles("
#include \"mpif.h\"
program test
end program"
	HAVE_MPIF_H
	SRC_EXT .F90
	FFLAGS ${CROCO_FORTRAN_CPP_FLAGS})

######################################################
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDFF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MPI DEFAULT_MSG
HAVE_MPIF_H)

######################################################
mark_as_advanced(HAVE_MPIF_H)