# define SED_TOY
# define SED_TOY_FLOC_0D
/* SED TOY — Flocculation 0D (5 x 5 x 50) */
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

# undef  SEDIMENT
# define MUSTANG

# ifdef MUSTANG
#   define key_MUSTANG_flocmod
#   define GLS_MIXING
#   define GLS_KOMEGA
# endif

# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD
#  define FLOC_TURB_DISS
#  undef FLOC_BBL_DISS
#  define SED_FLOCS
#  undef SED_DEFLOC
# endif

# undef  MORPHODYN
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
