# define OBSTRUCTION

# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_ADV
# define UV_COR
# define NONLIN_EOS
# define SALINITY
# define ANA_GRID
# define MASKING
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define GLS_MIXING
# define PSOURCE
# undef  PSOURCE_MASS
# define ANA_PSOURCE
# define NS_PERIODIC
# define NO_FRCFILE
# define USE_CALENDAR
# define NEW_S_COORD
# define OBC_EAST
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


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
