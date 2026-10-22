#include "cppdefs.h"
      module param
! The KXX common block is AGRIF-managed in the model sources (see param.h):
! do not declare it inside module param.
#define PARAM_MODULE_INCLUDE
#include "param.h"
      end module param
