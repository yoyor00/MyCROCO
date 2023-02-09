#include "cppdefs.h"
#ifdef M3FAST
      module module_param
# include "param.h"
      end module module_param
#else
      module module_param_empty
      end module module_param_empty
#endif
