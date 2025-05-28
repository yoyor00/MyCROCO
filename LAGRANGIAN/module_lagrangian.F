#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to lagrangian

#if defined LAGRANGIAN || defined DEB_IBM

      MODULE module_lagrangian

   ! pass from fix (F) to free format (F90)

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"
#include "ocean3d.h"
#include "ocean2d.h"
#include "mixing.h"
#include "coupling.h"
#include "set_global_definitions.h"

      END MODULE module_lagrangian
#else
      MODULE module_lagrangian_empty
      END MODULE module_lagrangian_empty
      
#endif /* LAGRANGIAN */