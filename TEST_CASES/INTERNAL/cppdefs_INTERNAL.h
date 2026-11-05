# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define UV_COR
# define UV_ADV
# define BODYTIDE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_VMIX
# define EW_PERIODIC
# define NS_PERIODIC
# ifdef INTERNALSHELF
#  undef   EW_PERIODIC
#  define OBC_EAST
#  define OBC_WEST
#  define SPONGE
#  define ANA_SSH
#  define ANA_M2CLIMA
#  define ANA_M3CLIMA
#  define ANA_TCLIMA
#  define ZCLIMATOLOGY
#  define M2CLIMATOLOGY
#  define M3CLIMATOLOGY
#  define TCLIMATOLOGY
#  define M2NUDGING
#  define M3NUDGING
#  define TNUDGING
# endif
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
