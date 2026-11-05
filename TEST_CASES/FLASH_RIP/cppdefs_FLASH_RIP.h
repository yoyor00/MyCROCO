# undef  MPI
# undef  NC4PAR
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define NO_TRACER
# define NO_TEMPERATURE
# undef  PASSIVE_TRACER
# undef  ANA_TIDES
# define NBQ
# define NBQ_PRECISE
# define LIMIT_BSTRESS
# define WAVE_MAKER
# define WAVE_MAKER_SPECTRUM
# define WAVE_MAKER_DSPREAD
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# define GLS_MIXING_3D
# define GLS_KOMEGA
# define NS_PERIODIC
# define OBC_WEST
# define OBC_SPECIFIED_WEST
# define SPONGE
# define FRC_BRY
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define WET_DRY
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define AVERAGES
# define AVERAGES_K
# undef  DIAGNOSTICS_EDDY


# define PGF_BASIC_JACOBIAN

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
