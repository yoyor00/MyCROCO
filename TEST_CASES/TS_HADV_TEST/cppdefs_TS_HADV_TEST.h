/*
!                       Horizontal TRACER ADVECTION EXAMPLE
!                       ========== ====== ========= =======
!
!  Velocity mode is selected at runtime via testcase_name in the nml:
!    'TS_HADV_TEST'      : body-forced (no velocity override),  dt=7
!    'TS_HADV_TEST_ROT'  : solid-body rotation (freeze u/v),    dt=4
!    'TS_HADV_TEST_DIAG' : constant diagonal advection,         dt=4
!    'TS_HADV_TEST_PER'  : cosine-modulated solid-body,         dt=4
*/
# undef  OPENMP
# undef  MPI
# undef  UV_ADV
# define NEW_S_COORD
# undef  UV_COR
# define SOLVE3D
# define M2FILTER_NONE
# define ANA_VMIX
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SSFLUX
# define NO_FRCFILE
# define SALINITY
# define EW_PERIODIC
# define NS_PERIODIC

/* Choose specific advection scheme */
# define TS_HADV_UP3
# undef  TS_HADV_C4
# undef  TS_HADV_UP5
# undef  TS_HADV_WENO5
# undef  TS_HADV_C6


# define PGF_FLAT_BOTTOM

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
