      MODULE module_HYBIOSED

#include "cppdefs.h"

#if defined HYBIOSED
      !!========================================================================
      !!                   ***  MODULE  module_HYBIOSED  ***
      !!        
      !!  give access to dimensions, variables, parameters from 
      !!  hydro model to HYBIOSED
      !! 
      !!  pass variables from fix (F) to free format (F90) using a module
      !!========================================================================

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "coupling.h"

#endif /* HYBIOSED */

      END MODULE module_HYBIOSED