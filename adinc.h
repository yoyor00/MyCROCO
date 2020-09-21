c -*- fortran -*-

c common includes for AD tools

# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# ifdef MPI
# include "mpi_cpl.h"
# ifdef AMPI
# include "ampi/ampif.h"
# else
      include 'mpif.h'
# endif
# else
      integer mynode, nnodes
      parameter (mynode=0)
      parameter (nnodes=1)
#endif

#ifndef MASKING
#define rmask(X,Y) 1
#endif

#include "adcppdefs.h"
