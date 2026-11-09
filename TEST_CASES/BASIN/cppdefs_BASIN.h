# undef  OPENMP
# undef  MPI
# define UV_ADV
# define UV_COR
# define UV_VIS2
# define SOLVE3D
# define TS_DIF2
# define BODYFORCE
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE
# undef  STOGEN
# define PGF_FLAT_BOTTOM

/* XIOS output frequency must be an integer multiple of dt (9600 s). */
# define CROCO_XIOS_OUTPUT_FREQ "9600s"

#include "cppdefs_dev.h"
#include "set_global_definitions.h"
