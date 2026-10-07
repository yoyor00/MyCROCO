/*
!                       Rip Current Example
!                       === ======= =======
!
!   Weir, B., Uchiyama, Y. (2010):
!      A vortex force analysis of the interaction of rip
!      currents and surface gravity wave
!      JGR Vol. 116
!
!  Idealized Grand Popo Beach in Benin (GRANDPOPO).
!
!  WAVE_MAKER & NBQ : wave-resolving (#undef MRL_WCI)
*/
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# undef  NBQ
# ifdef NBQ
#  define NBQ_PRECISE
#  define WAVE_MAKER
#  define WAVE_MAKER_SPECTRUM
#  define WAVE_MAKER_DSPREAD
#  define UV_HADV_WENO5
#  define UV_VADV_WENO5
#  define W_HADV_WENO5
#  define W_VADV_WENO5
#  define GLS_MIXING_3D
#  define GLS_KOMEGA
#  undef  MRL_WCI
#  define OBC_SPECIFIED_WEST
#  define FRC_BRY
#  define ANA_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
#  define AVERAGES
#  define AVERAGES_K
# else
#  define UV_VIS2
#  define UV_VIS_SMAGO
#  define LMD_MIXING
#  define LMD_SKPP
#  define LMD_BKPP
#  define MRL_WCI
# endif
# ifdef MRL_WCI
#  define WKB_WWAVE
#  define WKB_OBC_WEST
#  define WAVE_ROLLER
#  define WAVE_FRICTION
#  define WAVE_STREAMING
#  define MRL_CEW
#  define ANA_BRY_WKB
# else
#  define WAVE_MAKER_JONSWAP
#  undef  WAVE_MAKER_GAUSSIAN
# endif
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define OBC_WEST
# define SPONGE
# define WET_DRY
# define NS_PERIODIC
# undef SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
#  undef  MORPHODYN
# endif
# undef  DIAGNOSTICS_UV
# define PGF_BASIC_JACOBIAN


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
