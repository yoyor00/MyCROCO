#include "cppdefs.h"

        module module_parameter_oa

        USE param
        integer :: imax,jmax,kmax



#ifdef MPI
# define imax Lmmpi
# define jmax Mmmpi
#else
# define imax Lm
# define jmax Mm
#endif   


        end module module_parameter_oa
 
