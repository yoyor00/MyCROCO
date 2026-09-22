#define ISOLITON
#define ISOLITON_DJL
#define MPI
#define KNBQ3
#define KNBQ
#define SOLVE3D
#define NEW_S_COORD
#define UV_ADV
#define TS_HADV_WENO5
#define TS_VADV_WENO5
#define UV_HADV_WENO5
#define UV_VADV_WENO5
#define W_HADV_WENO5
#define W_VADV_WENO5
#define SALINITY
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define EW_PERIODIC
#define NO_FRCFILE

#include "cppdefs_dev.h"

/*
! Retour au schema de couplage stable SCH0 et desactivation de UV_VADV_WENO5_INTC6
! specifiquement pour ISOLITON afin d'eviter l'explosion (blow-up) du modele.
*/
# undef  UV_VADV_WENO5_INTC6
# undef  K3FAST_COUPLING_SCH1
# undef  K3FAST_COUPLINGW_SCH1
# define K3FAST_COUPLING_SCH0
# define K3FAST_COUPLINGW_SCH0

#include "set_global_definitions.h"
