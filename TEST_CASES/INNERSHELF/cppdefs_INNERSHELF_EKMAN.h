# undef  OPENMP
# undef  MPI
# undef  NBQ
# define SOLVE3D
# define UV_COR
# define ANA_GRID
# define ANA_INITIAL
# define AVERAGES
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define NS_PERIODIC
# define OBC_WEST
# define SPONGE
# define NO_FRCFILE
# define NO_SPONGE_GRID
# define NO_RESET_RHO0
# define NO_LIMIT_BSTRESS


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
