#define CANON2D
#define CANON2D_HWI
#undef  KH_INSTY
#undef  KH_INST3D
#undef  MPI
#undef  NBQ
#define KNBQ3
#define KNBQ
#define NBQ_GRAV
#define DIAG_CFL
#undef  NBQ_PRECISE
#undef  XIOS
#define SOLVE3D
#define NEW_S_COORD
#define UV_ADV
#define TS_HADV_WENO5
#define TS_VADV_WENO5
#define UV_HADV_WENO5
#define UV_VADV_WENO5
#define W_HADV_WENO5
#define W_VADV_WENO5
#undef  SALINITY
#undef  PASSIVE_TRACER
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#undef  ANA_SRFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define EW_PERIODIC
#define NO_FRCFILE
#undef  CVTK_DEBUG

#include "cppdefs_dev.h"

/*
! Surcharge pour stabiliser CANON2D et retrouver les resultats 2025
*/
# undef  UV_VADV_WENO5_INTC6
# undef  K3FAST_AB3
# undef  NBQ_MASS
# define K3FAST_AM4d
# undef  K3FAST_COUPLING_SCH1
# undef  K3FAST_COUPLINGW_SCH1
# define K3FAST_COUPLING_SCH0
# define K3FAST_COUPLINGW_SCH0
# define K3FAST_PG2

#include "set_global_definitions.h"
