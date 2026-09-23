/* Dune test case example */
# undef  OPENMP
# undef  MPI
# define M2FILTER_NONE
# define UV_ADV
# define NEW_S_COORD
# undef  UV_COR
# define SOLVE3D
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_STFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_SMFLUX
# define OBC_WEST
# define OBC_EAST
# define ANA_SSH
# define ZCLIMATOLOGY
# define ANA_M2CLIMA
# define M2CLIMATOLOGY
# define SEDIMENT
# undef  MUSTANG
# define MORPHODYN
# ifdef SEDIMENT
#  undef  SUSPLOAD
#  define BEDLOAD
#  undef  BEDLOAD_WENO5
#  define BEDLOAD_WULIN
#  define TAU_CRIT_WULIN
# endif
# ifdef MUSTANG
#  define key_MUSTANG_V2
#  define key_MUSTANG_bedload
#  define key_tauskin_c_upwind
#  define key_noTSdiss_insed
# endif
# define GLS_MIXING
# define NO_FRCFILE
# define OBC_COM_ZCHAPMAN
# define OBC_COM_M2SPECIFIED_WEST
# define OBC_COM_M2CHARACT_EAST


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
