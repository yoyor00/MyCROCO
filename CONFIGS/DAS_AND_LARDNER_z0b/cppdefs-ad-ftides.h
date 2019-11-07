! -*- fortran -*-
! $Id: cppdefs.h 1628 2015-01-10 13:53:00Z marchesiello $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#define FRICTION_TIDES
!-------------------------------------------------
! PRE-SELECTED OPTIONS
!
! ADVANCED OPTIONS ARE IN CPPDEFS_DEV.H
!-------------------------------------------------
#if defined FRICTION_TIDES
!     ADJ   EXAMPLE
!     ===== =======
# define Z0B_VAR ZOBI_VAR
# define OUT_DOUBLE
# define TAPENADE      

# undef OPENMP
# define MPI
# define AMPI

# undef UV_ADV
# undef UV_COR
# undef UV_VIS2
# undef SOLVE3D
# undef TS_DIF2
# undef BODYFORCE

# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX

# undef OBC_EAST
# undef OBC_WEST
# define OBC_NORTH
# undef OBC_SOUTH

#endif

# ifdef MPI
#  undef  PARALLEL_FILES
# endif

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
