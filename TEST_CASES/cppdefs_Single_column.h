# define SINGLE_COLUMN
/*
!                       Single Column Example
!                       ====== ====== =======
!
!                              Seven sets up are encompassed:
*/
/* Choose ONE experiment (activate via cppkeys in single-column.jsonc): */
# define KATO_PHILLIPS        /* erosion of linear strat by constant wind stress         (5x5x100) */
# undef  WILLIS_DEARDORFF     /* erosion of linear strat by constant surf buoyancy loss  (5x5x50)  */
# undef  DIURNAL_CYCLE        /* erosion of linear strat by diurnal surface forcing      (5x5x150) */
# undef  FORCED_EKBBL         /* forced Ekman bottom boundary layer                      (5x5x40)  */
# undef  FORCED_DBLEEK        /* forced Ekman bottom and surface boundary layers         (5x5x25)  */
# undef  FORCED_NONROTBBL     /* non-rotating forced bottom boundary layer (Prandtl)     (5x5x50)  */
# undef  FORCED_OSCNONROTBBL  /* non-rotating oscillatory forced bottom boundary layer   (5x5x50)  */
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
