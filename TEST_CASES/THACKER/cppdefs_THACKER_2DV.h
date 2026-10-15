/*
!                       Thacker Example
!                       ======= =======
!
! Thacker, W., (1981), Some exact solutions to the nonlinear
! shallow-water wave equations. J. Fluid Mech., 107, 499-508.
*/
# undef  OPENMP
# undef  MPI

# define SOLVE3D
# define UV_COR
# define UV_ADV
# undef  UV_VIS2
# define WET_DRY
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE
# define PGF_FLAT_BOTTOM


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
