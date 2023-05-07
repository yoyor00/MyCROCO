#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to OBSTRUCTIONS

#if defined OBSTRUCTION

      module module_OBSTRUCTIONS

   ! pass from fix (F) to free format (F90)

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"
# include "ocean3d.h"
# include "mixing.h"

      end module module_OBSTRUCTIONS
#else
      module module_OBSTRUCTIONS_empty
      end module module_OBSTRUCTIONS_empty
      
#endif /* OBSTRUCTIONS */