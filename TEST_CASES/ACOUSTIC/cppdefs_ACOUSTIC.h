# define ACOUSTIC
# undef  MPI
# undef  NBQ
# define KNBQ3
# define KNBQ
# ifdef NBQ
#  undef  NBQ_PRECISE
# endif
# ifdef KNBQ
#  define K3FAST_SACOUS
#  define K3FAST_DIAGACOUS
#  undef  K3FAST_SEDLAYERS
#  define K3FAST_CSVISC2K
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
# undef  RVTK_DEBUG


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
