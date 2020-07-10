! $Id: module_substance  2018-08- bthouvenin 
#include "cppdefs.h"

#if defined SUBSTANCE
      module module_substance

#include "param.h"
#include "ncscrum.h"
#include "scalars.h"
#include "forces.h"
#include "grid.h"

      end module module_substance
#else
      module module_substance_empty
      end module module_substance_empty
      
#endif /* SUBSTANCE */

        !!==============================================================================

