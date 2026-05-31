# define SANDBAR
# define SANDBAR_OFFSHORE
/* LIP-1C */
# undef  SANDBAR_ONSHORE
# undef  OPENMP
# undef  MPI
# undef  NBQ
# define SOLVE3D
# define UV_ADV
# define NEW_S_COORD
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
/* ! NBQ */
# ifndef NBQ 
#  define MRL_WCI
#  ifdef MRL_WCI
#   define WKB_WWAVE
#   define MRL_CEW
#   define WKB_OBC_WEST
#   define WAVE_ROLLER
#   define WAVE_FRICTION
#   define WAVE_BREAK_TG86
#   define WAVE_STREAMING
#   undef  WAVE_RAMP
#  endif
#  define GLS_MIXING
#  define GLS_KOMEGA
#  undef  LMD_MIXING
#  ifdef LMD_MIXING
#   define LMD_SKPP
#   define LMD_BKPP
#  endif
#  define BBL
#  define BBL_BREAKING_STIR
/* NBQ */
# else
#  define MPI
#  define NBQ_PRECISE
#  define WAVE_MAKER
#  define UV_ADV
#  define UV_HADV_WENO5
#  define UV_VADV_WENO5
#  define W_HADV_WENO5
#  define W_VADV_WENO5
#  define GLS_MIXING_3D
#  define GLS_KOMEGA
#  define ANA_BRY
#  define Z_FRC_BRY
#  define M2_FRC_BRY
#  define M3_FRC_BRY
#  define T_FRC_BRY
#  define AVERAGES
#  define AVERAGES_K
#  define DIAGNOSTICS_EDDY
/* NBQ */
# endif 
# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  ifndef NBQ
#   define BEDLOAD
#  endif
#  define MORPHODYN
#  define TCLIMATOLOGY
#  define TNUDGING
#  define ANA_TCLIMA
# endif
# undef  STATIONS
# ifdef STATIONS
#  define ALL_SIGMA
# endif
# undef  DIAGNOSTICS_TS
# ifdef DIAGNOSTICS_TS
#  define DIAGNOSTICS_TS_ADV
# endif
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
