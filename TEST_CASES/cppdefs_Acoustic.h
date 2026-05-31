# define ACOUSTIC
# undef  MPI
# define NBQ
# ifdef NBQ
#  undef  NBQ_PRECISE
# endif
# undef  UV_VIS2
# define SOLVE3D
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SRFLUX
# define ANA_BTFLUX
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
