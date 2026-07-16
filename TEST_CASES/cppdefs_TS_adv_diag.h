# define TS_HADV_TEST
/*
!                       Horizontal TRACER ADVECTION EXAMPLE
!                       ========== ====== ========= =======
!
*/
/* Choose ONE velocity mode (activate via cppkeys in ts-hadv-test.jsonc):  */
/* Body-forced (default): all three undef                   dt=7           */
# undef  SOLID_BODY_ROT
/* Constant advection in the diagonal                       dt=4           */
# define  DIAGONAL_ADV
/* Space/time-varying velocity (solid-body periodic)        dt=4           */
# undef  SOLID_BODY_PER

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


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
