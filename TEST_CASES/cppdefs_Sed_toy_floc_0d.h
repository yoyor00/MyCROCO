# define SED_TOY
/*
!                       SED TOY (1D Single Column example)
!                       === === === ====== ====== ========
!
*/ 
/* Choose an experiment :               */
/*   Rouse                              */
# undef  SED_TOY_ROUSE
/*   Consolidation                      */ 
# undef  SED_TOY_CONSOLID
/*   Erosion and sediment resuspension  */
# undef  SED_TOY_RESUSP
/*   Flocculation                       */
# define  SED_TOY_FLOC_0D
/*   Flocculation                       */
# undef  SED_TOY_FLOC_1D

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

# ifdef SED_TOY_ROUSE
#  define ANA_VMIX
#  define BODYFORCE
# endif

# ifdef SED_TOY_FLOC_1D
#  define ANA_VMIX
#  define BODYFORCE
# endif

# ifdef SED_TOY_FLOC_0D
#  define ANA_VMIX
#  define BODYFORCE
# endif

# undef  SEDIMENT
# define  MUSTANG

# ifdef MUSTANG
#  if defined SED_TOY_FLOC_0D || defined SED_TOY_FLOC_1D
#   define key_MUSTANG_flocmod
#   define GLS_MIXING
#   define GLS_KOMEGA
#  endif
# endif


# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD

#  ifdef SED_TOY_ROUSE
#   define SED_TAU_CD_CONST
#  endif

#  if defined SED_TOY_FLOC_1D || defined SED_TOY_CONSOLID || \
	defined SED_TOY_RESUSP
#   undef  BBL
#   define GLS_MIXING
#   define GLS_KOMEGA
#   define MIXED_BED
#   undef  COHESIVE_BED
#  endif

#  if defined SED_TOY_FLOC_0D || defined SED_TOY_FLOC_1D
#   define FLOC_TURB_DISS
#   undef FLOC_BBL_DISS
#   define SED_FLOCS
#   undef SED_DEFLOC
#  endif

# endif
# undef  MORPHODYN
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
