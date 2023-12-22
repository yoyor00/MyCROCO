#include "cppdefs.h"
#if defined NBQ || defined K3FAST
      module module_ocean3d
# include "param.h"
# include "ocean3d.h"
      end module module_ocean3d
#else
      module module_ocean3d_empty
      end module module_ocean3d_empty
#endif
