#include "cppdefs.h"

!!! give access to dimensions, variables, parameters from hydro model to ibm

#if defined DEB_IBM 

      MODULE module_ibm

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

      END MODULE module_ibm
#else
      MODULE module_ibm_empty
      END MODULE module_ibm_empty
      
#endif /* DEB_IBM */