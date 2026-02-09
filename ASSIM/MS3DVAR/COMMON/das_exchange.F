#include "cppdefs.h"
#if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI

# define das_exchange_2d_tile das_exchange_Ir2d_tile
# define JSTART JstrR
# define ISTART IstrR
# include "das_exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef das_exchange_2d_tile

# define das_exchange_2d_tile das_exchange_Iu2d_tile
# define JSTART JstrR
# define ISTART Istr
# include "das_exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef exchange_2d_tile

# define das_exchange_2d_tile das_exchange_Iv2d_tile
# define JSTART Jstr
# define ISTART IstrR
# include "das_exchange_2d_tile.h"
# undef ISTART
# undef JSTART
# undef das_exchange_2d_tile

# ifdef SOLVE3D
#  define KSTART 1
#  define das_exchange_3d_tile das_exchange_r3d_tile
#  define JSTART JstrR
#  define ISTART IstrR
#  include "das_exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef das_exchange_3d_tile

#  define das_exchange_3d_tile das_exchange_u3d_tile
#  define JSTART JstrR
#  define ISTART Istr
#  include "das_exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef das_exchange_3d_tile

#  define das_exchange_3d_tile das_exchange_v3d_tile
#  define JSTART Jstr
#  define ISTART IstrR
#  include "das_exchange_3d_tile.h"
#  undef ISTART
#  undef JSTART
#  undef das_exchange_3d_tile
# endif
#else
      subroutine das_exchange_empty
      end
#endif




