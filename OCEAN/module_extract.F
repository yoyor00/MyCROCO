#include "cppdefs.h"

#if defined ONLINE_EXTRACT

      module module_extract

   ! pass from fix (F) to free format (F90)
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ncscrum.h"
# include "ocean2d.h"
# ifdef SOLVE3D
#  include "ocean3d.h"
# endif
# ifdef NC4PAR_OE
# include "mpi_cpl.h"
# endif
      end module module_extract
#else
      module module_extract_empty
      end module module_extract_empty
      
#endif /* ONLINE_EXTRACT */
