#define AgAc
#undef  MPI
#define AGRIF
#define AGRIF_2WAY
#undef  NC4PAR
#undef  K3FAST_HIS
#define KXX
#define KNBQ3
#define KNBQ
#undef  NBQ_NUDGING_W
#define K3FAST_SACOUS
#define K3FAST_DIAGACOUS
#undef  K3FAST_SEDLAYERS
#define K3FAST_CSVISC2K
#undef  XIOS
#define M2FILTER_NONE
#define SOLVE3D
#define UV_ADV
#define NEW_S_COORD
#define ANA_GRID
#define ANA_INITIAL
#undef  TS_HADV_WENO5
#undef  TS_VADV_WENO5
#undef  UV_HADV_WENO5
#undef  UV_VADV_WENO5
#undef  W_HADV_WENO5
#undef  W_VADV_WENO5
#undef  SALINITY
#undef  PASSIVE_TRACER
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define NO_FRCFILE
#define DIAG_CFL
#define K3FAST_NOBPG
#define K3FAST_SOFAR

#include "cppdefs_dev.h"

/*
! Surcharge pour stabiliser AgAc et retrouver les resultats 2025
*/
# undef  UV_VADV_WENO5_INTC6
# undef  K3FAST_COUPLING_SCH1
# undef  K3FAST_COUPLINGW_SCH1
# define K3FAST_COUPLING_SCH0
# define K3FAST_COUPLINGW_SCH0
# define K3FAST_AVG_CLASSIC

#include "set_global_definitions.h"
