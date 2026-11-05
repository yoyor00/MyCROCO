!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
! /* Buffer array to allow reshaping of input/output data fot
!  the purpose of compatibility with the SCRUM plotting package.
!  This implies that variables defined at RHO-, VORTICITY-, U-
!  and V-points are written into netCDF files in such a way as
!  if they would be dimensioned as follows:

!  Location     name                  dimensions

!    RHO-        r2dvar   zeta-type   (0:Lm+1,0:Mm+1)
!    VORT-       p2dvar   vort-type   (1:Lm+1,1:Mm+1)
!    U-          u2dvar   ubar-type   (1:Lm+1,0:Mm+1)
!    V-          v2dvar   vbar-type   (0:Lm+1,1:Mm+1)

!    RHO-,RHO-   r3dvar   RHO-type    (0:Lm+1,0:Mm+1,    N)
!    VORT-,RHO-  p3dvar               (1:Lm+1,1:Mm+1,    N)
!    U-,RHO-     u3dvar   U-type      (1:Lm+1,0:Mm+1,    N)
!    V-,RHO-     v3dvar   V-type      (0:Lm+1,1:Mm+1,    N)
!    RHO-,W-     w3dvar   W-type      (0:Lm+1,0:Mm+1,  0:N)
!    RHO-,BED-   b3dvar   BED-type    (0:Lm+1,0:Mm+1, NLAY)
!    ABL-        abl3dvar ABL-type    (0:Lm+1,0:Mm+1,N_abl)
! */
#include "cppdefs.h"

MODULE buffer

   USE param, only: N, Lm, Mm
#ifdef MUSTANG
   USE param, only: ksdmax
#endif
#ifdef SEDIMENT
   USE param, only: NLAY
#endif
#ifdef ABL1D
   USE param, only: N_abl
#endif

   IMPLICIT NONE

   ! default
   PRIVATE
   PUBLIC buff, max_buffer_N_size, init_buffer

   INTEGER :: max_buffer_N_size
   REAL, DIMENSION(:), ALLOCATABLE :: buff

CONTAINS

   SUBROUTINE init_buffer
      INTEGER :: tmp_size

      tmp_size = N + 1
#ifdef MUSTANG
      tmp_size = max(tmp_size, ksdmax)
#endif
#ifdef SEDIMENT
      tmp_size = max(tmp_size, NLAY)
#endif
#ifdef ABL1D
      tmp_size = max(tmp_size, N_abl + 1)
#endif
      max_buffer_N_size = tmp_size

      ALLOCATE (buff((Lm + 5)*(Mm + 5)*max_buffer_N_size))
   END SUBROUTINE init_buffer

END MODULE buffer
