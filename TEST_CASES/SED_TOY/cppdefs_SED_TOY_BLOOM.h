/* SED TOY — BLOOM biology + MUSTANG sediment (4 x 3 x 20 : 240) */
# undef  OPENMP
# undef  MPI
# define NEW_S_COORD
# define SOLVE3D
# undef  NONLIN_EOS
# define SALINITY
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
# define BODYFORCE
# define USE_CALENDAR

# define GLS_MIXING
# define GLS_KOMEGA

# undef  SEDIMENT
# define SUBSTANCE
# define MUSTANG

# define BIOLink
# define BIOLink_UPDATE_CONCBIO
# define BIOLink_PAR_eval

# define BLOOM
# define key_oxygen
# define key_BLOOM_insed
# define DIAGNOSTICS_BIO

# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
