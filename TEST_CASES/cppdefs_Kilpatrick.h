# define KILPATRICK
# define MPI
# define AVERAGES
# define NONLIN_EOS
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define NO_FRCFILE
# define ABL1D
# ifdef  ABL1D
#  define BULK_FLUX
#  undef  BULK_ECUMEV0
#  undef  BULK_ECUMEV6
#  define BULK_GUSTINESS
#  define ANA_ABL_LSDATA
#  define ANA_ABL_VGRID
#  define STRESS_AT_RHO_POINTS
#  undef  ABL_NUDGING
#  undef  ABL_NUDGING_DYN
#  undef  ABL_NUDGING_TRA
#  undef  ABL_DYN_RESTORE_EQ
#  undef  SFLUX_CFB
# else
#  undef BULK_FLUX
# endif


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
