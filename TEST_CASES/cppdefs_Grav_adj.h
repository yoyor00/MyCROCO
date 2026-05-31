# define GRAV_ADJ
# undef  OPENMP
# undef  MPI
# undef  NBQ
# undef  XIOS
# define SOLVE3D
# define NEW_S_COORD
# define UV_ADV
# define TS_HADV_WENO5
# define TS_VADV_WENO5
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# ifdef NBQ
#  define W_HADV_WENO5
#  define W_VADV_WENO5
# endif
# undef  UV_VIS2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# undef  PASSIVE_TRACER
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
