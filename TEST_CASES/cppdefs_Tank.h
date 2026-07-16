# define TANK
# undef  TANKY   /* Y-oriented tank: swap LLm0=1,MMm0=50 in param_Tank.h */
# undef  MPI
# define NBQ
# define NBQ_PRECISE
# define SOLVE3D
# undef  UV_ADV
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_BTFLUX
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
