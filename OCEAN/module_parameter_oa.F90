#include "cppdefs.h"

        module module_parameter_oa

        integer :: imax,jmax,kmax

        USE param


#ifdef MPI
# define imax Lmmpi
# define jmax Mmmpi
#else
# define imax Lm
# define jmax Mm
#endif   


        end module module_parameter_oa
 
