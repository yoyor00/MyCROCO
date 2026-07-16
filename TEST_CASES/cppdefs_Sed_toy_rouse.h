# define SED_TOY
# define SED_TOY_ROUSE
/* SED TOY — Rouse profile (5 x 5 x 100) */
# undef  OPENMP
# undef  MPI
# define NEW_S_COORD
# define SOLVE3D
# undef  NONLIN_EOS
# define SALINITY
# undef  UV_VIS2
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define EW_PERIODIC
# define NS_PERIODIC
# define ANA_VMIX
# define BODYFORCE
# define SEDIMENT
# undef  MUSTANG
# define SUSPLOAD
# undef  BEDLOAD
# define SED_TAU_CD_CONST
# undef  MORPHODYN
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
