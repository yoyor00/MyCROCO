# define TIDAL_FLAT
/*
!                       TIDAL_FLAT  Example
!                       ==========  =======
*/
# undef  OPENMP
# undef  MPI
# undef  NONLIN_EOS
# define NEW_S_COORD
# define SALINITY
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define UV_COR
# define SOLVE3D
# define UV_VIS2
# define GLS_MIXING
# define ANA_INITIAL
# define WET_DRY
# define TS_DIF2
# define SPONGE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define OBC_WEST
# define FRC_BRY
# ifdef FRC_BRY
#  define ANA_BRY
#  define Z_FRC_BRY
#  define OBC_M2CHARACT
#  define OBC_REDUCED_PHYSICS
#  define M2_FRC_BRY
#  undef  M3_FRC_BRY
#  define T_FRC_BRY
# endif
# undef  SEDIMENT
# define MUSTANG
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD
# endif
# ifdef MUSTANG
#  undef  key_MUSTANG_V2
# endif
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
