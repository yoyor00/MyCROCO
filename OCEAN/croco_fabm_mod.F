#include "cppdefs.h"
! Modified by SINTEF 2024
! ---------------------
#ifdef CFABM
      module croco_fabm
# include "param.h"
      ! == Added by SINTEF 2024
      USE fabm
      USE fabm_config

      
      ! ==
      end module croco_fabm
#else
      module croco_fabm_empty
      end module croco_fabm_empty
#endif/* CFABM */