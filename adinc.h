c -*- fortran -*-

c common includes for AD tools

# include "cppdefs.h"
# include "param.h"
# include "grid.h"
# include "scalars.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "mixing.h"
# include "coupling.h"
# include "ncscrum.h"
# include "private_scratch.h"
      real*QUAD buff(6)
      common /xyz/ buff
# include "climat.h"
# include "forces.h"
# ifdef MPI
# include "mpi_cpl.h"
# include "ampi/ampif.h"
# else
      integer nnodes
      parameter (nnodes=1)
#endif

#ifndef MASKING
#define rmask(X,Y) 1
#endif
