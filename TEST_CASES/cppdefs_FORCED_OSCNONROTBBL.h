# define SINGLE_COLUMN
# define  FORCED_OSCNONROTBBL
# undef  KATO_PHILLIPS
# define  SINGLE_COLUMN
/*
!                       Single Column Example
!                       ====== ====== =======
!
!                              Seven  sets up are encompassed :
*/
/* erosion of linear strat by constant wind stress */
# undef  KATO_PHILLIPS
/* erosion of linear strat by constant surf buoyancy loss */
# undef  WILLIS_DEARDORFF
/* erosion of linear strat by constant surf buoyancy loss */
# undef  DIURNAL_CYCLE
/* forced Ekman bottom boundary layer */
# undef  FORCED_EKBBL
/* forced Ekman bottom and surface boundary layers */
# undef  FORCED_DBLEEK
/* non rotating forced bottom boundary layer : Prandt layer */
# undef  FORCED_NONROTBBL
/* non rotating oscillatory forced bottom boundary layer */
# define  FORCED_OSCNONROTBBL
# undef  OPENMP
# undef  MPI
# define UV_ADV
# define NEW_S_COORD
# define UV_COR
# define SOLVE3D
# undef  LMD_MIXING
# define GLS_MIXING
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define EW_PERIODIC
# define NS_PERIODIC


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
