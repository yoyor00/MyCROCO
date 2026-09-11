/*  SWASH_GLOBEX_B3 : bichromatic (default): wa1=0.07, wa2=0.03 */
/*  SWASH_GLOBEX_B2 : alt bichromatic: wa1=0.09, wa2=0.01       */
/*  SWASH_GLOBEX_A3 : alt JONSWAP: wp=2.25, wa=0.0354, gamma=20 */
/*  SWASH_GLOBEX_A3 : you will have to  define WAVE_MAKER_JONSWAP else define WAVE_MAKER_BICHROMATIC*/

# define WAVE_MAKER_BICHROMATIC
# undef  OPENMP
# undef  MPI
# define SOLVE3D
# define AVERAGES
# define NBQ
# define NBQ_PRECISE
# define WAVE_MAKER
# define UV_ADV
# define UV_HADV_WENO5
# define UV_VADV_WENO5
# define W_HADV_WENO5
# define W_VADV_WENO5
# define GLS_MIXING_3D
# define GLS_KOMEGA
# define NEW_S_COORD
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_SRFLUX
# define ANA_SST
# define ANA_BTFLUX
# define OBC_WEST
# undef  OBC_SPECIFIED_WEST
# define ANA_BRY
# define Z_FRC_BRY
# define M2_FRC_BRY
# define M3_FRC_BRY
# define T_FRC_BRY
# define WET_DRY
# define NO_FRCFILE


#include "cppdefs_dev.h"
#include "set_global_definitions.h"
